# Build lagdynamics's `engagement` data set from the published Nestimate
# `trajectories` matrix. This script regenerates `data/engagement.rda`.
# Run with: Rscript data-raw/engagement.R
#
# The source matrix is shipped in the MIT-licensed Nestimate package
# (https://github.com/mohsaqr/Nestimate). We store it as a data frame so
# tna::tna(engagement) uses tna's sequence-data method directly.

stopifnot(requireNamespace("Nestimate", quietly = TRUE))

engagement <- as.data.frame(Nestimate::trajectories,
                            stringsAsFactors = FALSE)

stopifnot(
  is.data.frame(engagement),
  nrow(engagement) == 138L,
  ncol(engagement) == 15L,
  all(unique(unlist(engagement, use.names = FALSE)) %in%
        c("Active", "Average", "Disengaged", NA))
)

usethis::use_data(engagement, overwrite = TRUE)
