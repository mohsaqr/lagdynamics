# Build the imdb_genres sequence dataset for lagseq from the 1,000-movie
# IMDB subset shipped in the MIT-licensed cooccure package.
#
# The original cooccure::movies is tabular (one row per movie, with a
# comma-delimited `genres` field). For lag sequential analysis we need
# a temporal sequence; here we sort movies chronologically and take
# the primary (first-listed) genre of each. The result is a single
# 1,000-event sequence over a 22-genre alphabet.
#
# This serves as a medium-K (~20), medium-N (~1,000) stress test
# outside the learning-analytics domain that dominates the rest of
# the verification battery. It is validated against
# stats::chisq.test()$stdres (the base-R Haberman residual oracle),
# not against any published LSA result, since no LSA paper has been
# published on this slice.
#
# Run with: Rscript data-raw/imdb_genres.R

stopifnot(requireNamespace("cooccure", quietly = TRUE))

m <- cooccure::movies
stopifnot(all(c("startYear", "genres", "averageRating") %in% names(m)))

# Sort chronologically; within a year, sort by averageRating descending
# so ties resolve deterministically.
m <- m[order(m$startYear, -m$averageRating, m$genres), ]

# Primary genre: first token before the first comma.
primary <- vapply(strsplit(m$genres, ",", fixed = TRUE),
                  function(x) trimws(x[1]),
                  character(1))

imdb_genres <- list(
  sequence    = unname(primary),
  year        = m$startYear,
  decade      = m$decade,
  rating      = m$averageRating,
  title       = m$primaryTitle,
  alphabet    = sort(unique(primary)),
  source      = "cooccure::movies (MIT, https://github.com/mohsaqr/cooccure); IMDB extract 1970-2024, rating >= 7.0",
  license     = "MIT (cooccure package); IMDB raw data is Open Data",
  n_events    = length(primary),
  k_states    = length(unique(primary)),
  description = paste(
    "Chronological sequence of primary genres for 1,000 highly-rated",
    "IMDB movies (rating >= 7.0, votes >= 1000, 1970-2024), used as a",
    "medium-K medium-N validation input for lagseq. No published LSA",
    "exists for this slice; the test suite cross-validates lagseq's",
    "classical engine on it against stats::chisq.test()$stdres."
  )
)

save(imdb_genres, file = "data/imdb_genres.rda", compress = "bzip2")
cat("data/imdb_genres.rda saved (", file.size("data/imdb_genres.rda"),
    "bytes ); K =", imdb_genres$k_states, "  N =", imdb_genres$n_events, "\n")
