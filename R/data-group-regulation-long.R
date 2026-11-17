#' Regulation actions in collaborative learning, as an event log
#'
#' @description
#' `group_regulation_long` is a coded dataset, in long format, of regulation
#' during collaborative learning: each row is one coded action of one student,
#' with the student's collaboration group, course and achievement level and the
#' time of the action. It contains the nine actions of [group_regulation]
#' (`adapt`, `cohesion`, `consensus`, `coregulate`, `discuss`, `emotion`,
#' `monitor`, `plan`, `synthesis`) and serves to show how [lsa()] reads an event
#' log and fits and compares groups.
#'
#' The log holds 27,533 actions of 2,000 students in 200 collaboration
#' groups of 10 students each. Achievement (`High` or `Low`, 1,000
#' students each) is constant within a student, and the course (`A`, `B`
#' or `C`) is constant within a group, so either can be passed to the
#' `group` argument of [lsa()] for a between-group comparison.
#'
#' @format A data frame with 27,533 rows (actions) and 6 columns:
#' \describe{
#'   \item{`Actor`}{Integer. Student identifier (2,000 students).}
#'   \item{`Achiever`}{Character. Achievement level, `High` or `Low`.}
#'   \item{`Group`}{Numeric. Collaboration group (200 groups).}
#'   \item{`Course`}{Character. Course, `A`, `B` or `C`.}
#'   \item{`Time`}{POSIXct. Time of the action.}
#'   \item{`Action`}{Character. Regulation action (9 values).}
#' }
#'
#' @source Taken without modification from `tna::group_regulation_long`
#'   (\pkg{tna} package, MIT license).
#'
#' @seealso [group_regulation] for the wide sequence table of the same
#'   actions; [lsa()]; [compare_lsa()].
#'
#' @examples
#' fit <- lsa(group_regulation_long, actor = "Actor", action = "Action",
#'            time = "Time")
#' transitions(fit, direction = "over", sort = "strength") |> head()
#'
#' # One model per achievement level
#' gfit <- lsa(group_regulation_long, actor = "Actor", action = "Action",
#'             time = "Time", group = "Achiever")
#' summary(gfit)
#'
#' @keywords datasets
"group_regulation_long"
