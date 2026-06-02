# Bidirectional / matched-pair test. Implements §10 of
# inst/REFERENCES.md following Sackett (1979) and Wampold (1984).
#
# The test asks whether the unordered pair {i, j} co-occurs at lag 1
# more (or less) than expected under independence, ignoring direction.
# Achieved by symmetrizing the transition matrix: W = O + t(O), then
# computing adjusted residuals on the symmetric table.

.engine_bidirectional <- function(transitions,
                                  structural_zeros = NULL,
                                  alternative = c("two.sided", "greater", "less"),
                                  n_events = NULL,
                                  ...) {
  stopifnot(inherits(transitions, "lsa_transitions"))
  alternative <- match.arg(alternative)
  if (!is.null(structural_zeros)) {
    stop("Engine 'bidirectional' does not support structural_zeros. ",
         "Use engine = 'classical' for structural-zero handling.",
         call. = FALSE)
  }

  obs <- transitions$obs
  K   <- nrow(obs)
  N   <- transitions$n_transitions
  if (N == 0) stop("No transitions in input.", call. = FALSE)

  # Symmetrize. W is symmetric by construction; rowSums(W) = colSums(W) =
  # rowSums(O) + colSums(O) at each code (the marginal "incidence"
  # under either direction).
  W   <- obs + t(obs)
  R_W <- rowSums(W)
  N_W <- sum(W)

  # Expected under symmetric independence: rowSums of the symmetrized
  # table are equal, so E_W[i, j] = R_W[i] * R_W[j] / N_W.
  exp_mat <- outer(R_W, R_W) / N_W

  # Haberman's standardized residual on the symmetric table:
  # z[i, j] = (W[i, j] - E_W[i, j]) / sqrt(E_W[i, j] * (1 - p[i]) * (1 - p[j]))
  # where p[i] = R_W[i] / N_W. Symmetric in (i, j).
  p_marg <- R_W / N_W
  denom  <- exp_mat * outer(1 - p_marg, 1 - p_marg)
  z      <- matrix(NA_real_, K, K)
  ok     <- is.finite(denom) & denom > 0
  z[ok]  <- (W[ok] - exp_mat[ok]) / sqrt(denom[ok])

  p_val <- .normal_p(z, alternative)

  # Transitional probabilities for reference (computed on original obs,
  # not the symmetrized table) — same as classical.
  prob <- .row_normalize(obs, rowSums(obs))

  # Yule's Q on the symmetrized table.
  yulesq <- .yules_q(W, R_W, R_W, N_W)

  # No kappa for bidirectional (different statistical question).
  kappa_na <- matrix(NA_real_, K, K)
  labs <- rownames(obs)
  dimnames(kappa_na) <- list(labs, labs)

  # No tablewise LR test for this engine — the symmetric model is a
  # different fit. Provide NULL so $build_lsa_object handles it.
  list(
    obs              = obs,
    exp              = exp_mat,
    prob             = prob,
    adj_res          = z,
    p                = p_val,
    yules_q          = yulesq,
    kappa            = kappa_na,
    kappa_z          = kappa_na,
    kappa_p          = kappa_na,
    lrx2             = NULL,
    structural_zeros = structural_zeros,
    alternative      = alternative,
    n_events_used    = if (is.null(n_events)) NA_integer_ else n_events,
    symmetric_obs    = W,
    symmetric_exp    = exp_mat,
    marginals        = R_W
  )
}
