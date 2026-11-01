#' Collaborative Learning Self-Regulation Sequences
#'
#' A wide-format dataset of group regulation during collaborative
#' learning: each row is one group's session, each column an ordered
#' time point, and each cell a coded regulation action. Shorter sessions
#' are padded with `NA`. It is the flagship example for the package
#' vignette: 9 states over 2,000 sequences gives a realistically rich
#' transition network.
#'
#' The nine coded actions are `adapt`, `cohesion`, `consensus`,
#' `coregulate`, `discuss`, `emotion`, `monitor`, `plan`, and
#' `synthesis`.
#'
#' @format A `data.frame` with 2,000 rows (sequences) and 26 columns
#'   (ordered time points `T1` to `T26`), character-coded.
#'
#' @source Identical in content to `tna::group_regulation` (`tna`
#'   package, MIT license).
#'
#' @seealso [group_regulation_long] (the long event log of the same
#'   actions)
#'
#' @examples
#' fit <- lsa(group_regulation)
#' transitions(fit, significant = TRUE)
#'
#' @keywords datasets
"group_regulation"
