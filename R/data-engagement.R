#' Student Engagement Trajectories
#'
#' Wide-format categorical sequence data: 138 students observed over 15
#' weekly time points. Each row is one student; each column is a
#' week. Entries are the student's engagement state for that week,
#' one of `"Active"`, `"Average"`, `"Disengaged"`, or `NA` for missing
#' weeks.
#'
#' This is a standard small-K, multi-sequence example for lag
#' sequential analysis: K = 3 states, S = 138 sequences, mean sequence
#' length about 15. It exercises the wide-matrix input path of
#' [lsa_data()] and produces a stable transition pattern with clear
#' adjusted-residual signals.
#'
#' @format A character matrix with 138 rows and 15 columns.
#'
#' @source Derived without modification from the `trajectories` matrix
#'   in the `Nestimate` package
#'   (\url{https://github.com/mohsaqr/Nestimate}), which is MIT-licensed
#'   and produced by Saqr and collaborators as a synthetic engagement
#'   trajectory example. Re-shipped here for convenience and offline
#'   testing; both attribution and license are preserved.
#'
#' @examples
#' fit <- lsa(engagement, engine = "classical")
#' fit
#' fit$adj_res
#'
"engagement"
