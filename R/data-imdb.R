#' IMDB Primary-Genre Sequence (1970-2024)
#'
#' Chronological sequence of primary genres for 1,000 highly-rated
#' IMDB films (`averageRating >= 7.0`, `numVotes >= 1000`, release
#' years 1970-2024). Each event is one film's first-listed genre; the
#' sequence is sorted by `startYear` ascending, breaking ties by
#' descending rating.
#'
#' This dataset serves as a medium-K, medium-N validation input
#' outside the learning-analytics domain that dominates the rest of
#' the lagseq battery. It exercises K = 16 with N = 1,000 events and
#' has no published LSA result to validate against; instead, the
#' test suite cross-validates lagseq's classical engine against
#' [stats::chisq.test()] on this exact sequence.
#'
#' @format A named list with elements:
#' \describe{
#'   \item{sequence}{Character vector of length 1,000: the
#'     chronological primary-genre sequence.}
#'   \item{year}{Integer vector of corresponding release years.}
#'   \item{decade}{Character vector of decade labels
#'     ("1970s", "1980s", ...).}
#'   \item{rating}{Numeric vector of IMDB average ratings.}
#'   \item{title}{Character vector of primary movie titles.}
#'   \item{alphabet}{Sorted character vector of the 16 distinct
#'     genres in `sequence`.}
#'   \item{source}{Citation string.}
#'   \item{license}{MIT (cooccure package); IMDB raw data is Open
#'     Data.}
#'   \item{n_events, k_states, description}{Summary metadata.}
#' }
#'
#' @source Derived from `cooccure::movies` (MIT,
#'   \url{https://github.com/mohsaqr/cooccure}); IMDB raw data
#'   (\url{https://www.imdb.com/interfaces/}) is Open Data.
#'
#' @examples
#' fit <- lsa(imdb_genres$sequence, engine = "classical")
#' fit
#' # Top over-represented genre transitions
#' edges <- fit$edges[order(fit$edges$adj_res, decreasing = TRUE), ]
#' head(edges[, c("from", "to", "count", "adj_res")], 5)
#'
"imdb_genres"
