# Non-parallel-dominance analysis. Implements §12 of
# inst/REFERENCES.md following Sackett (1979) and Wampold (1984)
# Table 3.
#
# Same dominance numerator as the parallel engine, but standard error
# is based on observed counts rather than expected:
#
#   D[i, j]    = O[i, j] - O[j, i]
#   SE_np[i,j] = sqrt( O[i, j] + O[j, i] )
#   z[i, j]    = D[i, j] / SE_np[i, j]
#
# Under H0 of within-pair symmetry, O[i,j] | (O[i,j] + O[j,i]) is
# Binomial(O[i,j] + O[j,i], 1/2). The z statistic is the normal
# approximation to that binomial; an exact binomial p-value is also
# computed and returned in $binomial_p.

.engine_nonparallel_dominance <- function(transitions,
                                          structural_zeros = NULL,
                                          alternative = c("two.sided",
                                                          "greater",
                                                          "less"),
                                          n_events = NULL,
                                          ...) {
  stopifnot(inherits(transitions, "lsa_transitions"))
  alternative <- match.arg(alternative)

  obs <- transitions$obs
  K   <- nrow(obs)
  R   <- rowSums(obs)
  C   <- colSums(obs)
  N   <- transitions$n_transitions
  if (N == 0) stop("No transitions in input.", call. = FALSE)

  exp_mat <- outer(R, C) / N

  D    <- obs - t(obs)
  W    <- obs + t(obs)              # pair totals
  se   <- sqrt(W)
  z    <- matrix(NA_real_, K, K)
  ok   <- is.finite(se) & se > 0
  z[ok] <- D[ok] / se[ok]
  diag(z) <- 0

  p_val <- .normal_p(z, alternative)

  # Exact binomial test on each cell (i, j) where pair total > 0.
  # Convention: alternative "two.sided" tests P(i->j) != 1/2 given the
  # pair total. The test is symmetric so binomial_p[i,j] == binomial_p[j,i].
  binomial_p <- matrix(NA_real_, K, K)
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      n_ij <- W[i, j]
      if (i == j || n_ij == 0) next
      x <- obs[i, j]
      alt <- switch(alternative,
                    "two.sided" = "two.sided",
                    "greater"   = "greater",
                    "less"      = "less")
      bt <- stats::binom.test(x, n_ij, p = 0.5, alternative = alt)
      binomial_p[i, j] <- bt$p.value
    }
  }
  dimnames(binomial_p) <- dimnames(obs)

  prob <- .row_normalize(obs, R)
  yulesq <- .yules_q(obs, R, C, N)

  na_K <- matrix(NA_real_, K, K)
  labs <- rownames(obs)
  dimnames(na_K) <- list(labs, labs)

  list(
    obs              = obs,
    exp              = exp_mat,
    prob             = prob,
    adj_res          = z,
    p                = p_val,
    yules_q          = yulesq,
    kappa            = na_K,
    kappa_z          = na_K,
    kappa_p          = na_K,
    lrx2             = NULL,
    structural_zeros = structural_zeros,
    alternative      = alternative,
    n_events_used    = if (is.null(n_events)) NA_integer_ else n_events,
    dominance        = D,
    pair_totals      = W,
    binomial_p       = binomial_p
  )
}
