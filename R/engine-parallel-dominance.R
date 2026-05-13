# Parallel-dominance analysis. Implements §11 of inst/REFERENCES.md
# following Sackett (1979) and tested in Wampold (1984), Table 2.
#
# For an ordered pair (i, j), tests whether i -> j is more likely
# than j -> i under the null of symmetric independence:
#
#   D[i, j]  = O[i, j] - O[j, i]
#   SE[i, j] = sqrt( E[i, j] + E[j, i] )   (independence expectation)
#   z[i, j]  = D[i, j] / SE[i, j]
#
# Antisymmetric: z[i, j] = -z[j, i]. Diagonal is zero by construction.

.engine_parallel_dominance <- function(transitions,
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

  # Independence expected (same as classical, used here for SE).
  exp_mat <- outer(R, C) / N

  # Dominance statistic: directional difference and its expected-SE form.
  D  <- obs - t(obs)
  S2 <- exp_mat + t(exp_mat)         # symmetric, equals expected-pair sum
  se <- sqrt(S2)
  z  <- matrix(NA_real_, K, K)
  ok <- is.finite(se) & se > 0
  z[ok] <- D[ok] / se[ok]
  diag(z) <- 0                       # D[i,i] = 0, SE > 0 -> z = 0

  p_val <- .normal_p(z, alternative)

  # Transitional probabilities for reference.
  prob <- .row_normalize(obs, R)

  # Yule's Q on the original asymmetric table (interpretation as
  # "association" still meaningful, but the engine's main test is z).
  yulesq <- .yules_q(obs, R, C, N)

  # No kappa for dominance.
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
    dominance_se     = se
  )
}
