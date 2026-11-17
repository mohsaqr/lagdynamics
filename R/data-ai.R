#' Coded AI behaviours in human-AI programming sessions
#'
#' @description
#' `ai_long` is a coded dataset, in long format, of the behaviours of an AI
#' coding assistant in "vibe-coding" sessions, in which a person builds software
#' by instructing the assistant in natural language. Each row is one coded AI
#' behaviour, with the project and session it belongs to and its position in the
#' session. A single assistant message can carry several codes, which are
#' ordered within the message.
#'
#' The log contains 8,551 coded events in 428 sessions from 34 projects.
#' Eight behaviour codes are grouped into three clusters: `Action`
#' (`Delegate`, `Execute`, `Investigate`), `Communication` (`Ask`,
#' `Explain`, `Report`) and `Repair` (`Plan`, `Repair`). With `session =
#' "session_id"`, `action = "code"` and `order = "order_in_session"`,
#' [lsa()] reads each session as one sequence and tests which behaviour
#' tends to follow which.
#'
#' @format A data frame with 8,551 rows (coded events) and 9 columns:
#' \describe{
#'   \item{`message_id`}{Integer. Identifier of the assistant message.}
#'   \item{`project`}{Character. Project identifier (34 projects).}
#'   \item{`session_id`}{Character. Session identifier (428 sessions,
#'     each within one project).}
#'   \item{`timestamp`}{Integer. Unix time in seconds.}
#'   \item{`session_date`}{Character. Session date as `YYYY-MM-DD`.}
#'   \item{`code`}{Character. AI behaviour code (8 values).}
#'   \item{`cluster`}{Character. Behaviour cluster (3 values).}
#'   \item{`code_order`}{Integer. Position of the code within its
#'     message.}
#'   \item{`order_in_session`}{Integer. Position of the event within its
#'     session.}
#' }
#'
#' @source Taken without modification from `Nestimate::ai_long`
#'   (\pkg{Nestimate} 0.7.7, MIT license).
#'
#' @seealso [lsa()] and its long-format arguments; [group_regulation_long]
#'   for another event log.
#'
#' @examples
#' fit <- lsa(ai_long, session = "session_id", action = "code",
#'            order = "order_in_session")
#' transitions(fit, direction = "over", sort = "strength")
#'
#' @keywords datasets
"ai_long"
