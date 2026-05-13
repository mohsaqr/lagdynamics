#' Canonical LSA Worked Example (O'Connor 1999)
#'
#' The complete published input/output pair from the canonical
#' lag-sequential-analysis methods paper: O'Connor, B. P. (1999),
#' "Simple and flexible SAS and SPSS programs for analyzing
#' lag-sequential categorical data", *Behavior Research Methods,
#' Instruments, & Computers*, 31(4), 718-726.
#' \doi{10.3758/BF03200753}.
#'
#' The paper publishes a 393-event input sequence (couple-interaction
#' data from Gottman & Roy 1990, p. 78; Appendix A) plus the full
#' numerical output of the SEQUENTIAL program for that exact input
#' (Appendix B). Six behavior codes; 392 transitions at lag 1.
#'
#' This is the gold-standard validation oracle for any lag-sequential
#' analysis implementation. lagseq's classical engine reproduces every
#' cell of every output matrix to better than the paper's own printed
#' precision (3-4 decimal places); see
#' `tests/testthat/test-published-oconnor.R`.
#'
#' @format A named list with elements:
#' \describe{
#'   \item{sequence}{Integer vector of length 393. Codes 1-6.}
#'   \item{obs}{6 x 6 integer matrix: published transition frequency
#'     matrix (`sum(obs) = 392`).}
#'   \item{expected}{6 x 6 numeric matrix: published expected
#'     frequencies under row x column independence.}
#'   \item{prob}{6 x 6 numeric matrix: published transitional
#'     probabilities.}
#'   \item{lrx2}{List `(statistic, df, p)` for the tablewise
#'     likelihood-ratio chi-square: statistic = 202.5009, df = 25,
#'     p approximately 0.}
#'   \item{adj_res}{6 x 6 numeric matrix of published adjusted
#'     residuals.}
#'   \item{adj_p}{6 x 6 matrix of published p-values for the
#'     adjusted residuals.}
#'   \item{yules_q}{6 x 6 matrix of published Yule's Q values.}
#'   \item{kappa}{6 x 6 matrix of published Wampold-style
#'     unidirectional kappas.}
#'   \item{kappa_z}{6 x 6 matrix of published kappa z-scores.}
#'   \item{kappa_p}{6 x 6 matrix of published kappa p-values.}
#'   \item{permutation}{List with the permutation test outputs
#'     (`n_blocks = 10`, `n_perm_block = 1000`, plus `p_mean`,
#'     `p_high`, `p_low` matrices). Used by Step 5
#'     `permute_lsa()` validation.}
#'   \item{source, notes}{Citation and metadata.}
#' }
#'
#' @source O'Connor, B. P. (1999). \emph{Behavior Research Methods,
#'   Instruments, & Computers}, \emph{31}(4), 718-726.
#'   \doi{10.3758/BF03200753}. The underlying interaction data are
#'   from Gottman, J. M., & Roy, A. K. (1990), p. 78.
#'
#' @examples
#' fit <- lsa(oconnor_couple$sequence, engine = "classical")
#' # Reproduce every output matrix in Appendix B:
#' all.equal(unname(fit$obs),     unname(oconnor_couple$obs))
#' all.equal(unname(fit$exp),     unname(oconnor_couple$expected),
#'           tolerance = 1e-3)
#' all.equal(unname(fit$adj_res), unname(oconnor_couple$adj_res),
#'           tolerance = 1e-3)
#' all.equal(unname(fit$yules_q), unname(oconnor_couple$yules_q),
#'           tolerance = 1e-3)
#' all.equal(unname(fit$kappa),   unname(oconnor_couple$kappa),
#'           tolerance = 1e-3)
#'
"oconnor_couple"
