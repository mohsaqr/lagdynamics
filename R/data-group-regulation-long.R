#' Group Regulation Long Event Log
#'
#' A long-format event log of collaborative-learning regulation actions,
#' one row per coded event. It is the long companion to the wide
#' [group_regulation] sequence matrix: the same nine regulation actions,
#' recorded here with actor, session, timing, and grouping columns so the
#' package's long-format import, grouping, and between-group comparison
#' paths can be demonstrated on a single realistic data set.
#'
#' The nine actions are `adapt`, `cohesion`, `consensus`, `coregulate`,
#' `discuss`, `emotion`, `monitor`, `plan`, and `synthesis`. The `Achiever`
#' column (`High` / `Low`) is a recorded grouping variable suitable for
#' between-group analysis.
#'
#' @format A `data.frame` with 27,533 rows and 6 columns:
#' \describe{
#'   \item{Actor}{Integer student identifier (2,000 students).}
#'   \item{Achiever}{Achievement group, `High` or `Low`.}
#'   \item{Group}{Numeric collaboration-group identifier.}
#'   \item{Course}{Course identifier (`A`, `B`, or `C`).}
#'   \item{Time}{Event timestamp (`POSIXct`).}
#'   \item{Action}{The regulation action code.}
#' }
#'
#' @source Derived without modification from `tna::group_regulation_long`
#'   (`tna` package, MIT license).
#'
#' @seealso [group_regulation] (the wide sequence matrix of the same
#'   actions)
#'
#' @examples
#' fit <- lsa(group_regulation_long, actor = "Actor",
#'            action = "Action", time = "Time")
#' transitions(fit, significant = TRUE)
#'
#' @keywords datasets
"group_regulation_long"
