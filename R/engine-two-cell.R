# Two-cell test (2x2 cell collapse). Implements §9 of
# inst/REFERENCES.md following Bakeman & Gottman (1997) ch. 7.
#
# For each cell (i, j), collapse the K x K table to a 2x2:
#
#                 j         not j
#             +---------+-------------+
#       i     |   a     |     b       |
#             +---------+-------------+
#      !i     |   c     |     d       |
#             +---------+-------------+
#
#  a = O[i, j]
#  b = R[i] - O[i, j]
#  c = C[j] - O[i, j]
#  d = N - R[i] - C[j] + O[i, j]
#
# Computes per-cell odds ratio, log-OR with Wald standard error,
# Wald z, and corresponding p-value. Yule's Q is the standard 2x2
# Q statistic computed inline. When any of a, b, c, d is zero, a
# 0.5 continuity correction is applied for the log-OR and its SE
# (Haldane-Anscombe correction).

.engine_two_cell <- function(transitions,
                             structural_zeros = NULL,
                             alternative = c("two.sided", "greater", "less"),
                             n_events = NULL,
                             continuity = 0.5,
                             ...) {
  stopifnot(inherits(transitions, "lsa_transitions"))
  alternative <- match.arg(alternative)
  stopifnot(is.numeric(continuity), length(continuity) == 1L,
            continuity >= 0)

  obs <- transitions$obs
  K   <- nrow(obs)
  R   <- rowSums(obs)
  C   <- colSums(obs)
  N   <- transitions$n_transitions
  if (N == 0) stop("No transitions in input.", call. = FALSE)

  exp_mat <- outer(R, C) / N
  prob    <- .row_normalize(obs, R)
  yulesq  <- .yules_q(obs, R, C, N)

  # 2x2 collapse cell-by-cell (vectorized in matrix form).
  a <- obs
  b <- outer(R, rep(1, K)) - obs
  c <- outer(rep(1, K), C) - obs
  d <- N - outer(R, rep(1, K)) - outer(rep(1, K), C) + obs

  # Haldane-Anscombe continuity adjustment: add `continuity` (default 0.5)
  # wherever any of a, b, c, d is exactly zero.
  any_zero <- (a == 0) | (b == 0) | (c == 0) | (d == 0)
  a_adj <- a + continuity * any_zero
  b_adj <- b + continuity * any_zero
  c_adj <- c + continuity * any_zero
  d_adj <- d + continuity * any_zero

  odds_ratio <- (a_adj * d_adj) / (b_adj * c_adj)
  log_or     <- log(odds_ratio)
  log_or_se  <- sqrt(1 / a_adj + 1 / b_adj + 1 / c_adj + 1 / d_adj)
  log_or_z   <- log_or / log_or_se
  log_or_z[!is.finite(log_or_z)] <- NA_real_

  p_val <- .normal_p(log_or_z, alternative)

  labs <- rownames(obs)
  dimnames(odds_ratio) <- list(labs, labs)
  dimnames(log_or)     <- list(labs, labs)
  dimnames(log_or_se)  <- list(labs, labs)
  dimnames(log_or_z)   <- list(labs, labs)

  na_K <- matrix(NA_real_, K, K, dimnames = list(labs, labs))

  list(
    obs              = obs,
    exp              = exp_mat,
    prob             = prob,
    adj_res          = log_or_z,
    p                = p_val,
    yules_q          = yulesq,
    kappa            = na_K,
    kappa_z          = na_K,
    kappa_p          = na_K,
    lrx2             = NULL,
    structural_zeros = structural_zeros,
    alternative      = alternative,
    n_events_used    = if (is.null(n_events)) NA_integer_ else n_events,
    odds_ratio       = odds_ratio,
    log_or           = log_or,
    log_or_se        = log_or_se,
    continuity       = continuity
  )
}
