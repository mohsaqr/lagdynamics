#' Regulation actions in collaborative learning, as sequences
#'
#' @description
#' `group_regulation` is a coded dataset, in wide format, of regulation during
#' collaborative learning: each row is one sequence of coded actions and each
#' column an ordered position. The nine actions are `adapt`, `cohesion`,
#' `consensus`, `coregulate`, `discuss`, `emotion`, `monitor`, `plan` and
#' `synthesis`. Sequences are 2 to 26 actions long; shorter sequences are padded
#' with `NA`.
#'
#' With nine states and 2,000 sequences, the data give a dense transition
#' table in which most transitions depart from independence, and they are
#' used to show lag transition networks and group comparisons. The same
#' actions are available as an event log in [group_regulation_long].
#'
#' @format A data frame with 2,000 rows (sequences) and 26 character
#'   columns `T1` to `T26` (positions, in order).
#'
#' @source Identical to `tna::group_regulation` (\pkg{tna} package, MIT
#'   license).
#'
#' @seealso [group_regulation_long] for the event log; [engagement];
#'   [lsa()].
#'
#' @examples
#' fit <- lsa(group_regulation)
#' transitions(fit, direction = "over", sort = "strength") |> head()
#'
#' @keywords datasets
"group_regulation"
