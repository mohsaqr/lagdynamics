#' Grandmother Behaviour Transitions, Qi An et al. (2026)
#'
#' The published lag-1 transition matrix and adjusted-residual /
#' Yule's Q output matrices for grandmother behaviour across two
#' Beijing dual-income households, as printed in Tables 4 and 5 of
#' Qi An, W. Xing, Y. Wang, X. Li (2026), *Sustainability*, 18(5),
#' 2326, \doi{10.3390/su18052326}.
#'
#' The 10 behaviour codes are:
#' `CO` (Cooking), `DH` (Doing housework), `WO` (Working),
#' `CK` (Caring for kid), `SM` (Self-management), `EA` (Eating),
#' `CM` (Communicating), `ED` (Education), `RE` (Resting),
#' `UO` (Using object).
#'
#' This object is the cleanest mathematical oracle shipped with
#' lagseq: the published input is a complete transition matrix and
#' the published output is a complete residual matrix, so feeding
#' `$obs` into [lsa()] and comparing the result to `$adj_res` is a
#' direct correctness check. Cross-validation against
#' [stats::chisq.test()] on `$obs` is exact at floating-point
#' precision (< 1e-12).
#'
#' Four cells in the paper's published Table 5 print values that do
#' not match the math computed from the paper's own Table 4 input.
#' These are documented in `$known_typos` together with the values
#' the math actually yields, so that lagseq's tests can distinguish
#' the paper's transcription errors from any genuine engine
#' regression.
#'
#' @format A named list with elements:
#' \describe{
#'   \item{obs}{10 x 10 integer matrix of transition frequencies
#'     (the paper's Table 4). `sum(obs) = 1531`.}
#'   \item{adj_res}{10 x 10 numeric matrix of published Z-scores
#'     (Table 5, upper of each cell pair).}
#'   \item{yules_q}{10 x 10 numeric matrix of published Yule's Q
#'     values (Table 5, lower of each cell pair).}
#'   \item{known_typos}{Data frame of 4 cells where Table 5 disagrees
#'     with the math computed from Table 4. Columns: `from`, `to`,
#'     `paper_printed`, `math_computed`, `category`.}
#'   \item{code_descriptions}{Named character vector mapping each
#'     code to its plain-English description.}
#'   \item{source}{Citation string.}
#'   \item{license}{Licensing note. The numerical tables are facts
#'     and are not copyrightable; the paper itself is published Open
#'     Access.}
#'   \item{n_transitions, k_states, notes}{Summary metadata.}
#' }
#'
#' @source Qi An, Wanli Xing, Yuzhe Wang, & Xiuyu Li. (2026).
#'   Behavioural Trajectories and Spatial Responses: A Study on Lag
#'   Sequential Analysis and Design Framework for Elderly Caregivers
#'   in Chinese Dual-Earner Households. \emph{Sustainability},
#'   \emph{18}(5), 2326. \doi{10.3390/su18052326}.
#'
#' @examples
#' # Reproduce the paper's residual analysis from its own input matrix
#' fit <- lsa(qi2026_grandmother$obs, engine = "classical")
#' # Compare to the published Z-scores (modulo the documented typos)
#' fit$adj_res - qi2026_grandmother$adj_res
#' # The 4 paper typos:
#' qi2026_grandmother$known_typos
#'
"qi2026_grandmother"
