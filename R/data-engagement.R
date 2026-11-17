#' Engagement states of students over time
#'
#' @description
#' `engagement` is a coded dataset, in wide format, of engagement: 138 students,
#' each observed at 15 ordered time points and coded at each as `Active`,
#' `Average` or `Disengaged`. Missing time points are `NA`. Two students have no
#' observation at all, so a fit reports 136 sequences; the remaining sequences
#' have on average about 14 observed states.
#'
#' With three states and many short sequences, the data show the typical
#' pattern of a lag sequential analysis: each state tends to persist, and
#' moves between `Active` and `Disengaged` are rare. They are the main
#' example of the help pages.
#'
#' @format A data frame with 138 rows (students) and 15 character columns
#'   named `1` to `15` (time points, in order). Values are `"Active"`,
#'   `"Average"`, `"Disengaged"` or `NA`.
#'
#' @source Derived from the synthetic `trajectories` example data of the
#'   \pkg{Nestimate} package (<https://github.com/mohsaqr/Nestimate>, MIT
#'   license), with the states kept and the matrix stored as a data frame.
#'
#' @seealso [group_regulation] for a larger set of sequences; [lsa()].
#'
#' @examples
#' fit <- lsa(engagement)
#' fit
#' transitions(fit, significant = TRUE)
#'
#' @keywords datasets
"engagement"
