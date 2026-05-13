#' Knowledge-Graph Learning Logs (Du Jun 2026, 29 learners)
#'
#' Event sequences from 29 undergraduate learners interacting with a
#' knowledge-graph learning environment over an exam-preparation
#' session. Each sequence is a chronologically ordered character
#' vector of action codes drawn from
#' `c("L01", "L02", ..., "L12", "E")` where `L01..L12` are knowledge
#' nodes the learner visited (`L01` = "Civil disputes", `L02` =
#' "Civil litigation", ...) and `E` = exercise attempt.
#'
#' Learners are labeled by 10-digit numeric IDs (used as `names()`).
#' The attribute `"group"` is a named character vector classifying
#' each learner into the post-test performance group used by the
#' source paper: `"低"` (low), `"中"` (medium), or `"高"` (high).
#'
#' This data set is intended as a real-world third-party validation
#' input for [lsa()] and is shipped alongside [kg_lsa_oracle], which
#' contains the source paper's own published LSA results computed on
#' the same data.
#'
#' @format A named list of 29 character vectors with a `"group"`
#'   attribute (also a named character vector).
#'
#' @source Du, J. (2026). The dataset of "Sequential Behavioral
#'   Mechanisms Linking Learning Paths to Academic Performance in
#'   Knowledge Graph Environments". \emph{Mendeley Data} V1.
#'   \doi{10.17632/bdwcj7vw94.1}. License: CC BY 4.0.
#'   Re-distributed without modification.
#'
#' @examples
#' length(kg_logs)
#' head(kg_logs[[1]], 10)
#' table(attr(kg_logs, "group"))
#'
#' fit <- lsa(kg_logs, engine = "classical")
#' fit
#'
#' @seealso [kg_lsa_oracle] for the source paper's published LSA
#'   matrices used as a validation oracle.
"kg_logs"


#' Published LSA Results for the Knowledge-Graph Dataset
#'
#' The lag-sequential-analysis outputs published by Du Jun (2026) for
#' the same 29-learner dataset shipped in [kg_logs]. Used as an
#' independent third-party validation oracle for [lsa()].
#'
#' The published values were extracted verbatim from the spreadsheet
#' deposited at \doi{10.17632/bdwcj7vw94.1} (sheets
#' `"整体分析"`, `"低结果"`, `"中结果"`, `"高结果"`). The author used
#' GSEQ-style output (matrices labeled JNTF for joint transition
#' frequencies and ADJR for adjusted residuals) computed from the
#' wide-format sequence sheet `"整体数据"`.
#'
#' Note that the published total `sum(obs)` is 870 transitions; running
#' lagseq on the same wide-format sheet yields 871, an
#' off-by-one difference attributable to a minor undocumented
#' preprocessing step in the source paper. Cell-level agreement is
#' typically within 1-5 events out of 870, and adjusted-residual
#' agreement is within 0.5 in roughly 90% of cells. See the test file
#' `tests/testthat/test-published-kg.R` for the precise agreement
#' thresholds.
#'
#' @format A named list with four elements (`overall`, `low`, `mid`,
#'   `high`), each itself a list with:
#' \describe{
#'   \item{obs}{The published 13 x 13 transition-frequency matrix
#'     (JNTF), with rows = "given" codes and columns = "target"
#'     codes.}
#'   \item{adj_res}{The published 13 x 13 adjusted-residual matrix
#'     (ADJR), printed to two decimal places in the source.}
#' }
#'
#' @source Du, J. (2026). The dataset of "Sequential Behavioral
#'   Mechanisms Linking Learning Paths to Academic Performance in
#'   Knowledge Graph Environments". \emph{Mendeley Data} V1.
#'   \doi{10.17632/bdwcj7vw94.1}. License: CC BY 4.0.
#'
#' @seealso [kg_logs] for the raw event sequences.
"kg_lsa_oracle"
