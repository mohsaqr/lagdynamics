# Classical lag sequential analysis engine. Implements:
#   §2  Expected frequencies (independence, optionally with IPF)
#   §3  Transitional probabilities
#   §4  Adjusted residuals (Haberman; Christensen for structural zeros)
#   §5  p-values for adjusted residuals
#   §6  Yule's Q
#   §7  Unidirectional kappa
#   §8  Tablewise likelihood-ratio chi-square
# from inst/REFERENCES.md.
#
# Pure matrix-algebra implementation; no for-loops over cells.

.engine_classical <- function(transitions,
                              structural_zeros = NULL,
                              alternative = c("two.sided", "greater", "less"),
                              n_events = NULL,
                              ...) {
  stopifnot(inherits(transitions, "lsa_transitions"))
  alternative <- match.arg(alternative)

  obs <- transitions$obs
  K <- nrow(obs)
  rt <- transitions$row_totals
  ct <- transitions$col_totals
  N  <- transitions$n_transitions

  if (N == 0) {
    stop("No transitions in input.", call. = FALSE)
  }

  has_struct <- !is.null(structural_zeros)
  if (has_struct) {
    stopifnot(
      is.matrix(structural_zeros),
      nrow(structural_zeros) == K, ncol(structural_zeros) == K,
      all(structural_zeros %in% c(0, 1))
    )
  }

  # §2: expected frequencies under independence.
  if (!has_struct) {
    exp_mat <- outer(rt, ct) / N
    ipf_info <- NULL
  } else {
    ipf <- lsa_ipf(obs, structure = structural_zeros)
    if (!isTRUE(ipf$converged)) {
      warning(sprintf(
        paste0("IPF did not converge after %d iterations (max margin ",
               "diff = %.3g). Residuals and p-values for the ",
               "structural-zero fit may be unreliable."),
        ipf$iterations, ipf$max_margin_diff
      ), call. = FALSE)
    }
    exp_mat <- ipf$fit
    ipf_info <- ipf[c("iterations", "converged", "max_margin_diff")]
  }

  # §3: transitional probabilities.
  prob <- .row_normalize(obs, rt)

  # §4 + §5: adjusted residuals and their p-values.
  resid <- .adjusted_residuals(
    obs = obs, exp_mat = exp_mat, rt = rt, ct = ct, N = N,
    structural_zeros = structural_zeros
  )
  z <- resid$z
  p <- .normal_p(z, alternative)

  # §6: Yule's Q (per-cell 2x2 collapse).
  yulesq <- .yules_q(obs, rt, ct, N)

  # §7: unidirectional kappa (event-level expectations).
  # Wampold's transformed kappa uses *event* totals, not transition
  # totals. For a single sequence of length T, event count for code i
  # equals row_total(obs)[i] + 1{i is last event} and equivalently
  # equals col_total(obs)[i] + 1{i is first event}. When the caller
  # supplies the event vector we use exact event totals; when only a
  # transition matrix is available we approximate via the larger of
  # the row/col total (a safe upper bound that recovers the exact
  # event count for the cell containing the boundary code).
  event_totals <- attr(transitions, "event_totals_col")
  if (is.null(event_totals)) {
    event_totals <- pmax(rt, ct)
  }
  if (is.null(n_events)) {
    n_events <- sum(event_totals)
  }
  kappa_out <- .unidirectional_kappa(
    obs = obs, n_i = event_totals, n_j = event_totals, n = n_events
  )

  # §8: tablewise likelihood-ratio chi-square.
  lr <- .lrx2(obs, exp_mat, K = K, structural_zeros = structural_zeros)

  list(
    obs = obs,
    exp = exp_mat,
    prob = prob,
    adj_res = z,
    p = p,
    yules_q = yulesq,
    kappa = kappa_out$kappa,
    kappa_z = kappa_out$z,
    kappa_p = .normal_p(kappa_out$z, alternative),
    lrx2 = lr,
    ipf = ipf_info,
    structural_zeros = structural_zeros,
    alternative = alternative,
    n_events_used = n_events
  )
}

# --- helpers ------------------------------------------------------------

# §3: prob[i, j] = obs[i, j] / rowtot[i]; NA when row total is zero.
.row_normalize <- function(obs, rt) {
  prob <- obs / rt  # rt recycles by column
  prob[rt == 0, ] <- NA_real_
  prob
}

# §4: adjusted residuals.
.adjusted_residuals <- function(obs, exp_mat, rt, ct, N,
                                structural_zeros = NULL) {
  K <- nrow(obs)
  if (is.null(structural_zeros)) {
    # §4.1 Haberman: matrix form
    p_row <- rt / N
    p_col <- ct / N
    denom <- exp_mat * outer(1 - p_row, 1 - p_col)
    z <- matrix(NA_real_, K, K)
    ok <- is.finite(denom) & denom > 0
    z[ok] <- (obs[ok] - exp_mat[ok]) / sqrt(denom[ok])
    return(list(z = z, hat = NULL))
  }
  # §4.2 Christensen with structural zeros via design-matrix hat.
  # Non-estimable cells (structural zeros and cells with degenerate
  # denominators) are returned as NA, not 0: zero would falsely
  # encode "exactly expected" and produce p = 1 on cells where no
  # test is defined.
  S <- structural_zeros
  keep <- which(S == 1, arr.ind = TRUE)
  if (nrow(keep) == 0L) {
    return(list(z = matrix(NA_real_, K, K), hat = NULL))
  }
  rows_kept <- keep[, 1L]
  cols_kept <- keep[, 2L]
  X <- stats::model.matrix(
    ~ factor(rows_kept, levels = seq_len(K)) +
      factor(cols_kept, levels = seq_len(K))
  )
  ev <- exp_mat[cbind(rows_kept, cols_kept)]
  W <- diag(ev, nrow = length(ev))
  XtWX <- crossprod(X, W %*% X)
  # Rank-based singularity check: det() loses precision on larger
  # tables and can falsely accept near-singular matrices.
  if (qr(XtWX)$rank < ncol(XtWX)) {
    warning("Design matrix is singular under the supplied ",
            "structural-zero pattern; falling back to the ",
            "no-structural-zero residual formula. This is an ",
            "approximation only.", call. = FALSE)
    return(.adjusted_residuals(obs, exp_mat, rt, ct, N,
                               structural_zeros = NULL))
  }
  H <- X %*% solve(XtWX) %*% crossprod(X, W)
  h_diag <- diag(H)
  z <- matrix(NA_real_, K, K)
  denom <- ev * (1 - h_diag)
  ok <- is.finite(denom) & denom > 0
  z_vals <- rep(NA_real_, length(ev))
  z_vals[ok] <- (obs[cbind(rows_kept, cols_kept)][ok] - ev[ok]) /
    sqrt(denom[ok])
  z[cbind(rows_kept, cols_kept)] <- z_vals
  list(z = z, hat = h_diag)
}

# §5: normal-CDF p-values.
.normal_p <- function(z, alternative) {
  switch(alternative,
    "two.sided" = 2 * stats::pnorm(-abs(z)),
    "greater"   = stats::pnorm(z, lower.tail = FALSE),
    "less"      = stats::pnorm(z)
  )
}

# §6: Yule's Q from the cell-wise 2x2 collapse.
.yules_q <- function(obs, rt, ct, N) {
  a <- obs
  b <- outer(rt, rep(1, ncol(obs))) - obs
  c <- outer(rep(1, nrow(obs)), ct) - obs
  d <- N - outer(rt, rep(1, ncol(obs))) -
        outer(rep(1, nrow(obs)), ct) + obs
  num <- a * d - b * c
  den <- a * d + b * c
  out <- matrix(NA_real_, nrow(obs), ncol(obs))
  ok <- is.finite(den) & den > 0
  out[ok] <- num[ok] / den[ok]
  out
}

# §7: unidirectional kappa, event-level expectation.
# Matches the canonical convention used by O'Connor's SEQUENTIAL
# (1999, p. 720) and Bakeman & Quera's GSEQ: et[i,j] = n_i * n_j / n
# for ALL cells including the diagonal. The Wampold (1989)
# theoretical formula with et[i,i] = n_i * (n_i - 1) / n is the
# alternative but is not what published LSA tools actually compute.
.unidirectional_kappa <- function(obs, n_i, n_j, n) {
  K <- nrow(obs)
  et <- outer(n_i, n_j) / n
  var <- outer(n_i, n_j) *
         outer(n - n_i, n - n_j) /
         (n^2 * (n - 1))
  mmax <- pmin(matrix(n_i, K, K), matrix(n_j, K, K, byrow = TRUE))
  num  <- obs - et
  pos  <- num >= 0
  k    <- matrix(NA_real_, K, K)
  denom_pos <- mmax - et
  denom_neg <- et
  ok_pos <- pos & is.finite(denom_pos) & denom_pos != 0
  ok_neg <- !pos & is.finite(denom_neg) & denom_neg != 0
  k[ok_pos] <- num[ok_pos] / denom_pos[ok_pos]
  k[ok_neg] <- num[ok_neg] / denom_neg[ok_neg]
  z <- matrix(NA_real_, K, K)
  ok_z <- is.finite(var) & var > 0
  z[ok_z] <- num[ok_z] / sqrt(var[ok_z])
  list(kappa = k, z = z)
}

# §8: likelihood-ratio chi-square and its p-value.
.lrx2 <- function(obs, exp_mat, K, structural_zeros = NULL) {
  ok <- obs > 0 & is.finite(exp_mat) & exp_mat > 0
  g2 <- 2 * sum(obs[ok] * log(obs[ok] / exp_mat[ok]))
  if (is.null(structural_zeros)) {
    df <- (K - 1L)^2
  } else {
    # Quasi-independence df: (# estimable cells) - rank(row + col
    # design under the structural-zero pattern). For non-degenerate
    # patterns this equals (K-1)^2 - s_zeros, but if a pattern leaves
    # an entire row or column with no estimable cell, the
    # corresponding effect is non-identifiable and rank(X) drops,
    # giving the correctly larger residual df. Forcing df >= 1 (as
    # the previous formula did) is a silent misuse of pchisq().
    keep <- which(structural_zeros == 1, arr.ind = TRUE)
    n_estimable <- nrow(keep)
    if (n_estimable == 0L) {
      return(list(statistic = NA_real_, df = NA_integer_, p = NA_real_))
    }
    X <- stats::model.matrix(
      ~ factor(keep[, 1L], levels = seq_len(K)) +
        factor(keep[, 2L], levels = seq_len(K))
    )
    df <- n_estimable - qr(X)$rank
  }
  if (!is.finite(df) || df <= 0L) {
    return(list(statistic = g2, df = df, p = NA_real_))
  }
  p <- stats::pchisq(g2, df = df, lower.tail = FALSE)
  list(statistic = g2, df = df, p = p)
}
