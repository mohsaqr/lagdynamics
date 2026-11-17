# Iterative proportional fitting (IPF) for expected frequencies with
# structural zeros. Implements §2.2 of inst/REFERENCES.md following
# Wickens (1989) ch. 4, pp. 107-112.
#
# The function returns the fitted expected-frequency matrix E with the
# property that rowSums(E) == rowSums(O), colSums(E) == colSums(O),
# and E[i, j] == 0 wherever S[i, j] == 0.
#
# Validated against stats::loglin() — see tests/testthat/test-ipf.R.

#' Fit expected counts with structural zeros by iterative proportional fitting
#'
#' @description
#' `lsa_ipf()` fits the expected counts of a count table with structural zeros
#' by iterative proportional fitting. When some transitions cannot occur by
#' construction (structural zeros, such as self-transitions in a coding scheme
#' that merges repeated codes), the quasi-independence model fits expected
#' counts that are zero in the forbidden cells and reproduce the observed row
#' and column totals; it has no closed form and is found by iterative
#' proportional fitting (Wickens, 1989), which rescales the rows and columns of
#' the table in turn until both sets of totals match. Under ordinary
#' independence, the expected count of a transition is the product of its row
#' and column totals divided by the grand total, a formula that assigns expected
#' counts to impossible cells.
#'
#' [lsa()] calls `lsa_ipf()` when `loops = FALSE` or `structural_zeros` is
#' given, so it is needed directly only to inspect the expected table or to use
#' it outside [lsa()].
#'
#' @param obs Square numeric matrix of observed transition counts.
#' @param structure Square 0/1 matrix of the same size: `0` marks a
#'   structural zero, `1` an estimable cell. Default `NULL`: every cell is
#'   estimable, which gives the ordinary independence table.
#' @param tol Convergence tolerance on the largest difference between
#'   fitted and observed totals (default `1e-8`).
#' @param max_iter Maximum number of row-and-column passes (default `200`).
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and one of:
#' * `lagdynamics_bad_input`: `obs` contains missing, infinite or negative
#'   values.
#' * `lagdynamics_ipf_infeasible`: a row or column has observed
#'   transitions but every one of its cells is a structural zero.
#'
#' @return A list with the fitted `K x K` expected-count matrix (`fit`), the
#'   number of passes used (`iterations`), whether the totals converged
#'   within `max_iter` (`converged`), and the largest remaining difference
#'   between fitted and observed totals (`max_margin_diff`).
#'
#' @examples
#' obs <- matrix(c(0, 4, 6,
#'                 3, 0, 5,
#'                 7, 2, 0), nrow = 3, byrow = TRUE)
#' # Every cell estimable: the ordinary independence table
#' lsa_ipf(obs)
#'
#' # Self-transitions forbidden: zero expected counts on the diagonal,
#' # observed row and column totals kept
#' lsa_ipf(obs, structure = 1 - diag(3))
#'
#' @references
#' Wickens, T. D. (1989). \emph{Multiway contingency tables analysis for
#' the social sciences} (pp. 107-112). Lawrence Erlbaum.
#'
#' @seealso [lsa()] and its `loops` and `structural_zeros` arguments.
#'
#' @export
lsa_ipf <- function(obs, structure = NULL, tol = 1e-8, max_iter = 200L) {
  stopifnot(
    "`obs` must be a square numeric matrix with at least 2 rows" =
      is.matrix(obs) && is.numeric(obs) && nrow(obs) == ncol(obs) &&
      nrow(obs) >= 2L
  )
  # IPF is defined for a non-negative finite count table. NA/Inf would
  # propagate into the scaling passes and return a converged-looking but
  # meaningless fit; negative "counts" are outside the model entirely.
  if (anyNA(obs) || !all(is.finite(obs))) {
    .lsa_stop("`obs` must contain only finite counts (no NA, NaN, or Inf).",
         call. = FALSE)
  }
  if (any(obs < 0)) {
    .lsa_stop("`obs` must be a non-negative count matrix.", call. = FALSE)
  }
  K <- nrow(obs)
  # No structural zeros by default. Self-transitions are kept; the
  # diagonal is part of the model. Users who want to forbid
  # self-transitions must pass structure = 1 - diag(K) explicitly.
  if (is.null(structure)) structure <- matrix(1L, K, K)
  stopifnot(
    "`structure` must be a K x K matrix of 0/1 values" =
      is.matrix(structure) && nrow(structure) == K && ncol(structure) == K &&
      all(structure %in% c(0, 1))
  )
  rt <- rowSums(obs)
  ct <- colSums(obs)
  # Sanity: if a row has total > 0 but no estimable cell in it, IPF cannot
  # converge. Same for columns.
  if (any(rt > 0 & rowSums(structure) == 0)) {
    .lsa_stop(class = "ipf_infeasible", "IPF infeasible: a row has positive marginal but no estimable cells.",
         call. = FALSE)
  }
  if (any(ct > 0 & colSums(structure) == 0)) {
    .lsa_stop(class = "ipf_infeasible", "IPF infeasible: a column has positive marginal but no estimable cells.",
         call. = FALSE)
  }

  E <- structure * 1.0
  converged <- FALSE
  iter <- 0L
  max_diff <- Inf
  while (iter < max_iter) {
    iter <- iter + 1L

    # Row scaling
    rsE <- rowSums(E)
    scale_r <- ifelse(rsE > 0, rt / rsE, 0)
    E <- E * scale_r

    # Column scaling
    csE <- colSums(E)
    scale_c <- ifelse(csE > 0, ct / csE, 0)
    E <- sweep(E, 2L, scale_c, `*`)

    max_diff <- max(
      max(abs(rowSums(E) - rt)),
      max(abs(colSums(E) - ct))
    )
    if (max_diff < tol) {
      converged <- TRUE
      break
    }
  }
  list(
    fit = E,
    iterations = iter,
    converged = converged,
    max_margin_diff = max_diff
  )
}
