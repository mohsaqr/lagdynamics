# Iterative proportional fitting (IPF) for expected frequencies with
# structural zeros. Implements §2.2 of inst/REFERENCES.md following
# Wickens (1989) ch. 4, pp. 107-112.
#
# The function returns the fitted expected-frequency matrix E with the
# property that rowSums(E) == rowSums(O), colSums(E) == colSums(O),
# and E[i, j] == 0 wherever S[i, j] == 0.
#
# Validated against stats::loglin() — see tests/testthat/test-ipf.R.

#' Iterative Proportional Fitting for Two-Way Tables with Structural Zeros
#'
#' Fits expected cell frequencies under the row + column independence
#' model when some cells are constrained to be exactly zero
#' (structural zeros). The fitted table has the same row and column
#' marginals as `obs`, satisfies `E[i, j] = 0` wherever
#' `structure[i, j] = 0`, and is the maximum likelihood estimate under
#' the independence model restricted to the non-zero pattern.
#'
#' Implementation follows Wickens (1989), pp. 107-112: alternately
#' scale rows then columns of an initialized expected table until row
#' and column marginals converge to those of `obs` within `tol`.
#'
#' @param obs Numeric `K x K` matrix of observed counts.
#' @param structure Numeric `K x K` 0/1 matrix. A `1` means the cell is
#'   estimable; a `0` means the cell is a structural zero. Defaults to
#'   `1 - diag(K)` (the "no self-transitions" pattern).
#' @param tol Numeric. Convergence tolerance on marginal differences.
#'   Default `1e-8`.
#' @param max_iter Integer. Maximum number of row+column scaling
#'   passes. Default `200L`.
#'
#' @return A list with elements:
#' \describe{
#'   \item{fit}{The fitted `K x K` expected-frequency matrix.}
#'   \item{iterations}{Number of row+column passes used.}
#'   \item{converged}{Logical scalar: whether convergence was reached
#'     within `max_iter`.}
#'   \item{max_margin_diff}{Maximum absolute difference between observed
#'     and fitted marginals at termination.}
#' }
#'
#' @examples
#' obs <- matrix(c(0, 4, 6,
#'                 3, 0, 5,
#'                 7, 2, 0), nrow = 3, byrow = TRUE)
#' fit <- lsa_ipf(obs)
#' fit$fit                                       # fitted expected counts
#' all.equal(rowSums(fit$fit), rowSums(obs))     # TRUE
#' all.equal(colSums(fit$fit), colSums(obs))     # TRUE
#'
#' @references
#' Wickens, T. D. (1989). \emph{Multiway contingency tables analysis
#' for the social sciences}, pp. 107-112. Lawrence Erlbaum.
#'
#' @export
lsa_ipf <- function(obs, structure = NULL, tol = 1e-8, max_iter = 200L) {
  stopifnot(
    is.matrix(obs), is.numeric(obs),
    nrow(obs) == ncol(obs), nrow(obs) >= 2L
  )
  K <- nrow(obs)
  if (is.null(structure)) structure <- 1 - diag(K)
  stopifnot(
    is.matrix(structure), nrow(structure) == K, ncol(structure) == K,
    all(structure %in% c(0, 1))
  )
  rt <- rowSums(obs)
  ct <- colSums(obs)
  # Sanity: if a row has total > 0 but no estimable cell in it, IPF cannot
  # converge. Same for columns.
  if (any(rt > 0 & rowSums(structure) == 0)) {
    stop("IPF infeasible: a row has positive marginal but no estimable cells.",
         call. = FALSE)
  }
  if (any(ct > 0 & colSums(structure) == 0)) {
    stop("IPF infeasible: a column has positive marginal but no estimable cells.",
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
