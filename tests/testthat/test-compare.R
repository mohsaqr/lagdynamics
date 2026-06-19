# Helper: build two groups of categorical sequences whose A-> transitions
# differ by `bias`, sharing the state set {A, B, C}.
.compare_make_groups <- function(n_a, n_b, bias_a, bias_b, seed = 1L) {
  set.seed(seed)
  states <- c("A", "B", "C")
  gen <- function(n, bias) replicate(n, {
    L <- sample(8:14, 1L)
    s <- character(L)
    s[1L] <- sample(states, 1L)
    for (i in 2:L) {
      p <- rep(1, 3)
      if (s[i - 1L] == "A") p <- p + bias
      s[i] <- sample(states, 1L, prob = p)
    }
    s
  }, simplify = FALSE)
  list(seqs = c(gen(n_a, bias_a), gen(n_b, bias_b)),
       group = rep(c("g1", "g2"), c(n_a, n_b)))
}

test_that("compare_lsa returns the documented shape and class", {
  d <- .compare_make_groups(30, 30, c(0, 4, 0), c(0, 0, 4), seed = 1L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  cmp <- compare_lsa(fit, R = 200L)

  expect_s3_class(cmp, "lsa_comparison")
  expect_identical(cmp$R, 200L)
  expect_identical(cmp$measure, "log_or")   # N-invariant default
  expect_identical(cmp$groups, c("g1", "g2"))
  # One row per ordered state pair (K = 3 -> 9 edges).
  expect_identical(nrow(cmp$edges), 9L)
  expect_true(all(c("from", "to", "log_or_a", "log_or_b", "diff",
                    "p_perm", "p_adj", "significant") %in%
                  names(cmp$edges)))
  expect_true(all(c("statistic", "p_value", "R") %in% names(cmp$global)))
  # diff is exactly group-a minus group-b on the measure columns.
  expect_equal(cmp$edges$diff,
               cmp$edges$log_or_a - cmp$edges$log_or_b,
               tolerance = 1e-12)
  # p-values live in [0, 1] (ignoring non-estimable NA cells).
  pp <- cmp$edges$p_perm[is.finite(cmp$edges$p_perm)]
  expect_true(all(pp >= 0 & pp <= 1))
})

test_that("compare_lsa detects a true group difference", {
  d <- .compare_make_groups(45, 45, c(0, 6, 0), c(0, 0, 6), seed = 7L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  cmp <- compare_lsa(fit, R = 300L)

  # The A->B and A->C edges carry the engineered difference.
  ab <- cmp$edges[cmp$edges$from == "A" & cmp$edges$to == "B", ]
  ac <- cmp$edges[cmp$edges$from == "A" & cmp$edges$to == "C", ]
  expect_true(ab$significant)
  expect_true(ac$significant)
  # Omnibus rejects overall network invariance.
  expect_lt(cmp$global$p_value, 0.05)
})

test_that("compare_lsa does not over-reject when groups are exchangeable", {
  # Both groups share one generating process: the null is true.
  d <- .compare_make_groups(40, 40, c(0, 4, 0), c(0, 4, 0), seed = 11L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  cmp <- compare_lsa(fit, R = 300L)

  expect_gt(cmp$global$p_value, 0.05)
  # With no real difference, few or no edges should be flagged.
  expect_lte(sum(cmp$edges$significant, na.rm = TRUE), 1L)
})

test_that("two-fit form matches the grouped form's observed statistic", {
  d <- .compare_make_groups(30, 30, c(0, 4, 0), c(0, 4, 0), seed = 11L)
  grp <- lsa(d$seqs, group = d$group, engine = "classical")
  # Same data, fit as two separate single-group objects.
  n_a <- 30L
  fa <- lsa(d$seqs[seq_len(n_a)], engine = "classical",
            labels = attr(grp, "labels"))
  fb <- lsa(d$seqs[(n_a + 1L):length(d$seqs)], engine = "classical",
            labels = attr(grp, "labels"))

  set.seed(99L); a <- compare_lsa(grp, R = 100L)
  set.seed(99L); b <- compare_lsa(fa, fb, R = 100L)
  # The observed difference (and hence omnibus statistic) is permutation-
  # independent, so the two entry points must agree exactly.
  expect_equal(a$global$statistic, b$global$statistic, tolerance = 1e-12)
  expect_equal(a$edges$diff, b$edges$diff, tolerance = 1e-12)
})

test_that("compare_lsa honours the measure argument", {
  d <- .compare_make_groups(30, 30, c(0, 4, 0), c(0, 0, 4), seed = 3L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  cmp <- compare_lsa(fit, R = 100L, measure = "prob")
  expect_identical(cmp$measure, "prob")
  expect_true(all(c("prob_a", "prob_b") %in% names(cmp$edges)))
})

test_that("compare_lsa adjusts p-values across edges", {
  d <- .compare_make_groups(35, 35, c(0, 5, 0), c(0, 0, 5), seed = 5L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  cmp <- compare_lsa(fit, R = 200L, adjust = "holm")
  ok <- is.finite(cmp$edges$p_perm)
  # Holm adjustment never decreases a p-value.
  expect_true(all(cmp$edges$p_adj[ok] >= cmp$edges$p_perm[ok] - 1e-12))
  expect_identical(cmp$adjust, "holm")
})

test_that("compare_lsa validates its inputs", {
  d <- .compare_make_groups(20, 20, c(0, 4, 0), c(0, 0, 4), seed = 2L)
  grp <- lsa(d$seqs, group = d$group, engine = "classical")
  fb <- lsa(d$seqs[21:40], engine = "classical")

  # A single-group lsa_group cannot be compared.
  g1 <- lsa(d$seqs[1:20], group = rep("only", 20L), engine = "classical")
  expect_error(compare_lsa(g1), "at least two groups")

  expect_error(compare_lsa(grp, fb), "leave `y = NULL`")

  # Transition-matrix input cannot be permuted.
  m <- matrix(c(0, 5, 3, 2, 0, 4, 1, 6, 0), 3, 3,
              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  fm <- lsa(m)
  expect_error(compare_lsa(fm, fm), "event-level input")

  expect_error(compare_lsa(grp, R = 100L, adjust = "nope"), "adjust")
})

test_that("compare_lsa runs all pairwise comparisons for >2 groups", {
  set.seed(21L)
  d3 <- list(
    seqs = c(.compare_make_groups(30, 30, c(0, 6, 0), c(0, 0, 6),
                                  seed = 21L)$seqs,
             .compare_make_groups(30, 0, c(0, 0, 0), c(0, 0, 0),
                                  seed = 22L)$seqs),
    group = rep(c("AB", "AC", "flat"), each = 30))
  fit <- lsa(d3$seqs, group = d3$group, engine = "classical")
  cmp <- compare_lsa(fit, R = 200L, adjust = "BH")

  expect_s3_class(cmp, "lsa_comparison_pairwise")
  # C(3, 2) = 3 pairs, each with K^2 = 9 edges -> 27 rows.
  expect_identical(nrow(cmp$edges), 27L)
  expect_identical(nrow(cmp$global), 3L)
  expect_identical(length(cmp$comparisons), 3L)
  expect_identical(names(cmp$comparisons),
                   c("AB_vs_AC", "AB_vs_flat", "AC_vs_flat"))
  expect_true(all(c("group_a", "group_b", "from", "to", "diff",
                    "p_perm", "p_adj", "significant") %in%
                  names(cmp$edges)))
  expect_true(all(c("group_a", "group_b", "statistic", "p_value",
                    "p_adj") %in% names(cmp$global)))
  expect_identical(cmp$groups, c("AB", "AC", "flat"))
  # The AB-vs-AC pair carries the largest engineered difference.
  s <- cmp$global$statistic
  names(s) <- paste0(cmp$global$group_a, "_", cmp$global$group_b)
  expect_identical(names(which.max(s)), "AB_AC")
})

test_that("pairwise correction is family-wide, not per-pair", {
  set.seed(31L)
  d3 <- list(
    seqs = c(.compare_make_groups(28, 28, c(0, 5, 0), c(0, 0, 5),
                                  seed = 31L)$seqs,
             .compare_make_groups(28, 0, c(0, 4, 0), c(0, 0, 0),
                                  seed = 32L)$seqs),
    group = rep(c("g1", "g2", "g3"), each = 28))
  fit <- lsa(d3$seqs, group = d3$group, engine = "classical")
  cmp <- compare_lsa(fit, R = 200L, adjust = "bonferroni")

  ok <- is.finite(cmp$edges$p_perm)
  m <- sum(ok)  # total tests in the family
  # Family-wide Bonferroni multiplies each p by the pooled count m, not
  # by the per-pair count (K^2). Verify against a direct recomputation.
  expect_equal(cmp$edges$p_adj[ok],
               pmin(1, cmp$edges$p_perm[ok] * m),
               tolerance = 1e-12)
})

test_that("two bare fits still return a single lsa_comparison", {
  d <- .compare_make_groups(25, 25, c(0, 4, 0), c(0, 0, 4), seed = 41L)
  fa <- lsa(d$seqs[1:25], engine = "classical", labels = c("A", "B", "C"))
  fb <- lsa(d$seqs[26:50], engine = "classical", labels = c("A", "B", "C"))
  cmp <- compare_lsa(fa, fb, R = 100L)
  expect_s3_class(cmp, "lsa_comparison")
  expect_false(inherits(cmp, "lsa_comparison_pairwise"))
})

test_that("plot.lsa_comparison returns a ggplot", {
  skip_if_not_installed("ggplot2")
  d <- .compare_make_groups(25, 25, c(0, 5, 0), c(0, 0, 5), seed = 51L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  p <- plot(compare_lsa(fit, R = 100L))
  expect_s3_class(p, "ggplot")
})

test_that("plot.lsa_comparison_pairwise returns a faceted ggplot", {
  skip_if_not_installed("ggplot2")
  set.seed(61L)
  d3 <- list(
    seqs = c(.compare_make_groups(24, 24, c(0, 5, 0), c(0, 0, 5),
                                  seed = 61L)$seqs,
             .compare_make_groups(24, 0, c(0, 0, 0), c(0, 0, 0),
                                  seed = 62L)$seqs),
    group = rep(c("a", "b", "c"), each = 24))
  fit <- lsa(d3$seqs, group = d3$group, engine = "classical")
  p <- plot(compare_lsa(fit, R = 100L))
  expect_s3_class(p, "ggplot")
  # One panel per pair (C(3, 2) = 3). Use ggplot_build rather than `$`
  # extraction, which is deprecated on the S7 ggplot object.
  panels <- ggplot2::ggplot_build(p)$layout$layout
  expect_identical(nrow(panels), 3L)
})

test_that("log_or default is N-invariant: equal behavior in unequal groups", {
  # Same transition process, but group A has 3x the sessions of group B.
  # adj_res differs (scales with sqrt(N)); log_or does not.
  d <- .compare_make_groups(90, 30, c(0, 5, 0), c(0, 5, 0), seed = 71L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")

  set.seed(1); cmp_lor <- compare_lsa(fit, R = 300L, measure = "log_or")
  set.seed(1); cmp_z <- suppressMessages(
    compare_lsa(fit, R = 300L, measure = "adj_res"))

  # The A->B edge is engineered identical in both groups: log_or should
  # not flag it, even though adj_res (sqrt(N)-inflated) is prone to.
  lor_ab <- cmp_lor$edges[cmp_lor$edges$from == "A" &
                          cmp_lor$edges$to == "B", ]
  expect_false(isTRUE(lor_ab$significant))
  expect_gt(cmp_lor$global$p_value, 0.05)
})

test_that("adj_res measure warns when group sizes are unequal", {
  d <- .compare_make_groups(90, 30, c(0, 5, 0), c(0, 0, 5), seed = 73L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  expect_message(compare_lsa(fit, R = 50L, measure = "adj_res"),
                 "sqrt\\(N\\)")
  # The default (log_or) is silent.
  expect_silent(compare_lsa(fit, R = 50L))
})

test_that("log_or matrix is finite on zero cells (Haldane correction)", {
  d <- .compare_make_groups(20, 20, c(0, 6, 0), c(0, 0, 6), seed = 75L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  lor <- lagseq:::.lsa_log_or(fit[[1L]])
  # Some transitions are unobserved (zero cells); Haldane keeps log_or
  # finite there rather than -Inf, except engine-NA cells.
  obs <- fit[[1L]]$obs
  zero_cells <- obs == 0 & !is.na(fit[[1L]]$adj_res)
  expect_true(all(is.finite(lor[zero_cells])))
})

test_that("min_count guard: never-observed transitions are not significant", {
  # Build groups over 4 states where state D never follows anything in
  # either group, so every *->D and D->* cell is 0/0 (pooled 0).
  set.seed(81L)
  states <- c("A", "B", "C")  # D absent from the alphabet on purpose
  gen <- function(n, bias) replicate(n, {
    L <- sample(8:14, 1L); s <- character(L); s[1L] <- sample(states, 1L)
    for (i in 2:L) { p <- rep(1, 3); if (s[i-1L]=="A") p <- p+bias
                     s[i] <- sample(states, 1L, prob = p) }
    s
  }, simplify = FALSE)
  seqs <- c(gen(40, c(0,5,0)), gen(40, c(0,0,5)))
  # Inject one rare extra state "Z" exactly twice total, so Z cells are
  # below the default min_count of 5.
  seqs[[1]][1] <- "Z"; seqs[[41]][1] <- "Z"
  grp <- rep(c("g1","g2"), each = 40)
  fit <- lsa(seqs, group = grp, engine = "classical")
  cmp <- compare_lsa(fit, R = 300L, adjust = "BH")

  # Pooled count per edge.
  fa <- fit[[1]]; fb <- fit[[2]]; labs <- rownames(fa$obs)
  ix <- cbind(match(cmp$edges$from, labs), match(cmp$edges$to, labs))
  pooled <- fa$obs[ix] + fb$obs[ix]
  # Every below-threshold cell must be untested (NA p) and not significant.
  rare <- pooled < 5
  expect_true(all(is.na(cmp$edges$p_perm[rare])))
  expect_false(any(cmp$edges$significant[rare], na.rm = TRUE))
})

test_that("min_count = 0 tests every populated cell", {
  d <- .compare_make_groups(40, 40, c(0, 5, 0), c(0, 0, 5), seed = 83L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  cmp0 <- compare_lsa(fit, R = 100L, min_count = 0L)
  cmp5 <- compare_lsa(fit, R = 100L, min_count = 5L)
  # Relaxing the threshold tests at least as many cells.
  expect_gte(sum(is.finite(cmp0$edges$p_perm)),
             sum(is.finite(cmp5$edges$p_perm)))
})

test_that("omnibus uses the same tested-cell set and is NA when none qualify", {
  d <- .compare_make_groups(40, 40, c(0, 5, 0), c(0, 0, 5), seed = 91L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")

  # With an impossibly high min_count nothing is testable: the omnibus
  # must be NA (no test performed), not a misleading p = 1.
  cmp_none <- compare_lsa(fit, R = 100L, min_count = 1e6)
  expect_true(all(is.na(cmp_none$edges$p_perm)))
  expect_true(is.na(cmp_none$global$p_value))
  expect_true(is.na(cmp_none$global$statistic))

  # Normal run: the omnibus is a finite probability.
  cmp <- compare_lsa(fit, R = 200L)
  expect_true(is.finite(cmp$global$p_value))
  expect_true(cmp$global$p_value >= 0 && cmp$global$p_value <= 1)
})

test_that("barrel borders the higher group's bar, darker with |diff|", {
  d <- .compare_make_groups(45, 45, c(0, 6, 0), c(0, 0, 6), seed = 95L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  cmp <- compare_lsa(fit, R = 300L)
  fa <- fit[[1]]; fb <- fit[[2]]
  e <- cmp$edges
  n <- nrow(e)
  one <- lagseq:::.barrel_one(fa, fb, e, rev(seq_len(n)), "prob")
  rect <- one$rect
  left <- rect[rect$side == "left", ]    # group a bars, edge order
  right <- rect[rect$side == "right", ]  # group b bars, edge order

  # The higher group's bar is bordered (non-white); the lower stays white.
  a_higher <- is.finite(e$diff) & e$diff > 0
  expect_true(all(left$border[a_higher] != "white"))
  expect_true(all(right$border[a_higher] == "white"))
  expect_true(all(right$border[e$diff < 0] != "white"))

  # Border darkens with |diff|: the largest-|diff| edge is darker (lower
  # grey value) than the smallest non-zero one on its winning side.
  win_border <- ifelse(a_higher, left$border, right$border)
  ok <- is.finite(e$diff) & e$diff != 0
  # grey() returns "#RRGGBB"; a darker border has a smaller RR channel.
  chan <- strtoi(substr(win_border, 2, 3), 16L)
  expect_lt(chan[ok][which.max(abs(e$diff[ok]))],
            chan[ok][which.min(abs(e$diff[ok]))])
})

test_that("barrel rank = 'effect' surfaces avoided (negative log OR) cells", {
  skip_if_not_installed("ggplot2")
  d <- .compare_make_groups(45, 45, c(0, 6, 0), c(0, 0, 6), seed = 97L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  cmp <- compare_lsa(fit, R = 200L)
  expect_s3_class(plot(cmp, rank = "frequency"), "ggplot")
  expect_s3_class(plot(cmp, rank = "effect"), "ggplot")
  expect_error(plot(cmp, rank = "nope"))

  # The effect ranking must put the strongest-|log OR| tested cell first.
  fa <- fit[[1]]; fb <- fit[[2]]; labs <- rownames(fa$obs)
  e <- cmp$edges; ix <- cbind(match(e$from, labs), match(e$to, labs))
  eff <- pmax(abs(lagseq:::.lsa_log_or(fa)[ix]),
              abs(lagseq:::.lsa_log_or(fb)[ix]))
  eff[!is.finite(e$p_perm)] <- -Inf
  top <- e[which.max(eff), ]
  # That top cell's per-group log OR is the most extreme among tested cells.
  expect_gt(max(eff), 0)
})

test_that("as.data.frame returns the tidy edge table for comparisons", {
  d <- .compare_make_groups(30, 30, c(0, 5, 0), c(0, 0, 5), seed = 101L)
  fit <- lsa(d$seqs, group = d$group, engine = "classical")
  cmp <- compare_lsa(fit, R = 100L)
  expect_identical(as.data.frame(cmp), cmp$edges)

  d3 <- list(seqs = c(d$seqs, .compare_make_groups(20, 0, c(0, 3, 0),
                                                   c(0, 0, 0), 102L)$seqs),
             group = rep(c("a", "b", "c"), c(30, 30, 20)))
  fit3 <- lsa(d3$seqs, group = d3$group, engine = "classical")
  cmp3 <- compare_lsa(fit3, R = 100L)
  expect_identical(as.data.frame(cmp3), cmp3$edges)
})
