# Split-half network reliability for an lsa fit. Implements the
# random-half resampling described in Epskamp et al. (2018) for
# psychometric networks, adapted to LSA edge matrices.
#
# Each replicate: partition sequences uniformly at random into halves
# H1 and H2, refit lsa on each half, correlate the two edge-weight
# vectors. The distribution of these correlations across R replicates
# is the reliability summary.

#' Estimate the split-half reliability of a lag transition network
#'
#' @description
#' `reliability_lsa()` estimates the split-half reliability of a lag transition
#' network from one data set. A lag transition network is reliable when two
#' independent samples of the same population give similar networks. Split-half
#' reliability follows the resampling approach that Epskamp, Borsboom and Fried
#' (2018) describe for the accuracy of psychological networks: the sequences are
#' divided at random into two halves, a network is estimated on each half, and
#' the two networks are compared transition by transition. A high correlation
#' between the halves means that the structure of the network does not depend on
#' which sequences happened to be sampled.
#'
#' The random split is repeated `R` times. Each split assigns half of the
#' sequences to one half and the rest to the other, refits the model on each
#' with the same engine and settings, and correlates the two halves' transition
#' probabilities, counts or adjusted residuals over all transitions.
#' `reliability_lsa()` reports the mean, standard deviation and 95% interval
#' (2.5% and 97.5% quantiles) of the split-half correlations, and the number of
#' failed replicates, in which a half had no transitions or a half's values did
#' not vary. A grouped fit gives one result per group.
#'
#' Reliability describes the network as a whole. For the precision of
#' each transition use [bootstrap_lsa()] or [certainty_lsa()]; for
#' whether each transition stays significant when sequences are dropped,
#' use [stability_lsa()].
#'
#' @param fit An `lsa` fit from [lsa()], or a grouped fit from
#'   `lsa(..., group = )`, built from at least two sequences (not from a
#'   transition matrix).
#' @param R Number of random splits. Default `100`.
#' @param weights Transition values correlated between the halves:
#'   `"prob"` (default, transition probabilities), `"count"` or
#'   `"adj_res"` (adjusted residuals).
#' @param method Correlation: `"pearson"` (default) or `"spearman"`.
#' @param parallel Run the splits on several cores. Default `FALSE`.
#' @param n_cores Number of cores when `parallel = TRUE`. Default `NULL`,
#'   all cores but one.
#' @param verbose Report progress every 100 splits. Default `FALSE`.
#' @param ... For a grouped fit, the arguments above, passed on to each
#'   group; otherwise not used.
#' @param x An `lsa_reliability` or `lsa_reliability_group` object, for
#'   the `print()` methods.
#'
#' @return An object of class `c("lsa_reliability", "list")`. Read it with
#'   [as.data.frame()][as.data.frame.lsa_reliability]:
#'   `as.data.frame(x)` gives one row per split with its correlation
#'   (`NA` for a failed replicate), and `as.data.frame(x, what =
#'   "summary")` one row with the settings, the number of failed
#'   replicates `n_failed`, and the `mean`, `sd`, `ci_low` and `ci_high`
#'   of the correlations. `print()` shows the same summary and returns
#'   `x` invisibly. For a grouped fit, an object of class
#'   `c("lsa_reliability_group", "list")` with one result per group, whose
#'   tables gain a leading `group` column.
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and one of:
#' * `lagdynamics_needs_events`: the fit was built from a transition
#'   matrix, so there are no sequences to split.
#' * `lagdynamics_bad_input`: the fit has fewer than two sequences.
#'
#' @examples
#' fit <- lsa(engagement)
#' rel <- reliability_lsa(fit, R = 50)
#' rel
#' as.data.frame(rel, what = "summary")
#'
#' # One result per group
#' g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
#'          time = "Time", group = "Achiever")
#' reliability_lsa(g, R = 20) |>
#'   as.data.frame(what = "summary")
#' \donttest{
#' # A full analysis uses more splits
#' reliability_lsa(fit, R = 500)
#' }
#'
#' @references
#' Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
#' psychological networks and their accuracy: A tutorial paper.
#' \emph{Behavior Research Methods, 50}(1), 195-212.
#'
#' @seealso [as.data.frame.lsa_reliability()] to read the result;
#'   [bootstrap_lsa()], [certainty_lsa()], [stability_lsa()] and
#'   [permute_lsa()] for transition-level checks.
#'
#' @export
reliability_lsa <- function(fit, ...) UseMethod("reliability_lsa")

#' @rdname reliability_lsa
#' @export
reliability_lsa.lsa <- function(fit,
                                R = 100L,
                                weights = c("prob", "count", "adj_res"),
                                method = c("pearson", "spearman"),
                                parallel = FALSE,
                                n_cores = NULL,
                                verbose = FALSE,
                                ...) {
  weights <- match.arg(weights)
  method <- match.arg(method)
  stopifnot(
    "`R` must be a single whole number >= 1" =
      is.numeric(R) && length(R) == 1L && is.finite(R) && R >= 1L &&
      R == floor(R)
  )
  R <- as.integer(R)

  recipe <- fit$params
  d <- fit$data
  if (identical(d$source, "transitions")) {
    .lsa_stop(class = "needs_events", "reliability_lsa() requires event-level input. The fit was ",
         "built from a pre-computed transition matrix.", call. = FALSE)
  }
  per_seq <- split(d$events, d$seq_id)
  S <- length(per_seq)
  if (S < 2L) {
    .lsa_stop("reliability_lsa() needs at least 2 sequences; got ", S, ".",
         call. = FALSE)
  }

  pull_weights <- function(refitted) {
    switch(weights,
      prob    = refitted$prob,
      count   = refitted$obs,
      adj_res = refitted$adj_res
    )
  }

  worker <- function(b) {
    half <- max(1L, S %/% 2L)
    idx_1 <- sample.int(S, half, replace = FALSE)
    idx_2 <- setdiff(seq_len(S), idx_1)
    if (length(idx_2) == 0L) return(NA_real_)

    refit_half <- function(idx) {
      hs <- per_seq[idx]
      ev <- unlist(hs, use.names = FALSE)
      sid <- rep.int(seq_along(idx),
                     times = vapply(hs, length, integer(1)))
      # A half made up only of singleton sequences contributes zero
      # transitions, and the engine signals `lagdynamics_no_transitions`.
      # That is a property of the random split, not of the fit, so the
      # replicate is recorded as failed (NA) and counted in `n_failed`.
      # Any other error is a real problem and propagates.
      tryCatch(
        .refit_from_events(events = ev, seq_id = sid,
                           labels = d$labels, recipe = recipe),
        lagdynamics_no_transitions = function(e) NULL
      )
    }

    f1 <- refit_half(idx_1)
    f2 <- refit_half(idx_2)
    if (is.null(f1) || is.null(f2)) return(NA_real_)
    w1 <- as.vector(pull_weights(f1))
    w2 <- as.vector(pull_weights(f2))
    ok <- is.finite(w1) & is.finite(w2)
    # A correlation needs two finite pairs and variation in both halves;
    # otherwise the replicate is undefined (NA), not an error.
    if (sum(ok) < 2L || !(stats::sd(w1[ok]) > 0) ||
        !(stats::sd(w2[ok]) > 0)) {
      return(NA_real_)
    }
    stats::cor(w1[ok], w2[ok], method = method)
  }

  if (verbose) message("Running ", R, " split-half replicates ...")
  results <- .run_parallel(worker, R = R, parallel = parallel,
                           n_cores = n_cores, verbose = verbose)
  cors <- vapply(results, function(v) {
    if (is.numeric(v) && length(v) == 1L) v else NA_real_
  }, numeric(1))

  # Replicates whose split left a half with no transitions are failed
  # replicates; replicates with an undefined correlation are NA as well.
  n_failed <- sum(!is.finite(cors))
  finite <- is.finite(cors)
  ci <- if (sum(finite) >= 2L) {
    stats::quantile(cors[finite], c(0.025, 0.975), names = FALSE)
  } else {
    c(NA_real_, NA_real_)
  }
  structure(
    list(
      correlations = cors,
      mean         = if (any(finite)) mean(cors[finite]) else NA_real_,
      sd           = if (sum(finite) >= 2L) stats::sd(cors[finite])
                     else NA_real_,
      ci_low       = ci[1L],
      ci_high      = ci[2L],
      R            = R,
      weights      = weights,
      method       = method,
      n_sequences  = S,
      n_failed     = n_failed,
      fit          = fit
    ),
    class = c("lsa_reliability", "list")
  )
}

#' @rdname reliability_lsa
#' @export
reliability_lsa.lsa_group <- function(fit, ...) {
  rels <- lapply(fit, function(f) reliability_lsa(f, ...))
  names(rels) <- names(fit)
  structure(rels, levels = names(fit),
            class = c("lsa_reliability_group", "list"))
}

#' @rdname reliability_lsa
#' @export
print.lsa_reliability_group <- function(x, ...) {
  cat("<lsa_reliability_group>\n")
  cat(sprintf("  groups: %d\n\n", length(x)))
  for (i in seq_along(x)) {
    cat(sprintf("[%s]\n", names(x)[i]))
    print(x[[i]])
    if (i < length(x)) cat("\n")
  }
  invisible(x)
}

#' Read the tables of a split-half reliability result
#'
#' @description
#' `as.data.frame()` returns the tables of a [reliability_lsa()] result, which
#' holds the split-half correlation of every random split: these correlations,
#' one row per split, or a one-row summary of their distribution, chosen with
#' `what`. For a grouped result the tables of the groups are stacked.
#'
#' @param x An `lsa_reliability` or `lsa_reliability_group` object from
#'   [reliability_lsa()].
#' @param what `"replicates"` (default): one row per split. `"summary"`:
#'   one row summarizing all splits.
#' @param row.names,optional Not used; present for consistency with the
#'   generic.
#' @param ... For a grouped result, passed on to each group (use it for
#'   `what`); otherwise not used.
#'
#' @return A data frame.
#'   * `what = "replicates"`: one row per split, with `replicate` and
#'     `correlation`; a failed replicate (a half with no transitions, or
#'     no variation to correlate) is kept as `NA`.
#'   * `what = "summary"`: one row with `weights`, `method`, `R`,
#'     `n_failed`, `mean`, `sd`, `ci_low` and `ci_high` (2.5% and 97.5%
#'     quantiles of the finite correlations).
#'
#'   A grouped result gains a leading `group` column.
#'
#' @examples
#' rel <- reliability_lsa(lsa(engagement), R = 20)
#' as.data.frame(rel) |> head()
#' as.data.frame(rel, what = "summary")
#'
#' g <- lsa(group_regulation_long, actor = "Actor", action = "Action",
#'          time = "Time", group = "Achiever")
#' grel <- reliability_lsa(g, R = 10)
#' as.data.frame(grel, what = "summary")
#'
#' @seealso [reliability_lsa()].
#'
#' @export
as.data.frame.lsa_reliability <- function(x, row.names = NULL,
                                          optional = FALSE,
                                          what = c("replicates", "summary"),
                                          ...) {
  what <- match.arg(what)
  if (identical(what, "summary")) {
    return(data.frame(
      weights = x$weights, method = x$method, R = x$R,
      n_failed = x$n_failed %||% sum(!is.finite(x$correlations)),
      mean = x$mean, sd = x$sd, ci_low = x$ci_low, ci_high = x$ci_high,
      stringsAsFactors = FALSE, row.names = row.names
    ))
  }
  data.frame(replicate = seq_len(x$R),
             correlation = as.numeric(x$correlations),
             stringsAsFactors = FALSE, row.names = row.names)
}

#' @rdname as.data.frame.lsa_reliability
#' @export
as.data.frame.lsa_reliability_group <- function(x, row.names = NULL,
                                                optional = FALSE, ...) {
  .grouped_df(x, ...)
}

#' @rdname reliability_lsa
#' @export
print.lsa_reliability <- function(x, ...) {
  cat("<lsa_reliability>\n")
  cat(sprintf("  engine:        %s\n", x$fit$method))
  cat(sprintf("  replicates:    %d\n", x$R))
  cat(sprintf("  weights:       %s\n", x$weights))
  cat(sprintf("  method:        %s\n", x$method))
  cat(sprintf("  n sequences:   %d\n", x$n_sequences))
  cat(sprintf("  split-half r:  %.3f  (sd = %.3f)\n", x$mean, x$sd))
  cat(sprintf("  95%% CI:        [%.3f, %.3f]\n", x$ci_low, x$ci_high))
  n_failed <- x$n_failed %||% sum(!is.finite(x$correlations))
  cat(sprintf("  failed:        %d of %d replicates\n", n_failed, x$R))
  invisible(x)
}
