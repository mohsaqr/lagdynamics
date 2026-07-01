#' Human-AI Vibe Coding Interaction Events
#'
#' A long-format event log of coded AI-side behaviours in human-AI vibe
#' coding sessions. Each row is one coded AI event, with project/session
#' identifiers, ordering variables, a fine-grained AI behaviour code, and a
#' broader behaviour cluster. It is useful for demonstrating [lsa()]'s
#' long-format event-log import path.
#'
#' The eight AI behaviour codes are `Ask`, `Delegate`, `Execute`, `Explain`,
#' `Investigate`, `Plan`, `Repair`, and `Report`. The three clusters are
#' `Action`, `Communication`, and `Repair`.
#'
#' @format A `data.frame` with 8,551 rows and 9 columns:
#' \describe{
#'   \item{message_id}{Integer message identifier.}
#'   \item{project}{Project identifier.}
#'   \item{session_id}{Session identifier.}
#'   \item{timestamp}{Unix timestamp in seconds.}
#'   \item{session_date}{Session date as `YYYY-MM-DD`.}
#'   \item{code}{Fine-grained AI behaviour code.}
#'   \item{cluster}{Broader AI behaviour cluster.}
#'   \item{code_order}{Order of multiple codes within the same message.}
#'   \item{order_in_session}{Event order within session.}
#' }
#'
#' @source Derived without modification from `Nestimate::ai_long`
#'   (`Nestimate` version 0.7.7, MIT license).
#'
#' @examples
#' fit <- lsa(ai_long, actor = "project", session = "session_id",
#'            action = "code", order = "order_in_session")
#' transitions(fit, significant = TRUE)
#'
#' @keywords datasets
"ai_long"

