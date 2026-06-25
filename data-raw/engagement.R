# Build lagdynamics's `engagement` data set from the published Nestimate
# `trajectories` matrix. This script regenerates `data/engagement.rda`.
# Run with: Rscript data-raw/engagement.R
#
# The source matrix is shipped in the MIT-licensed Nestimate package
# (https://github.com/mohsaqr/Nestimate). It is a 138 x 15 character
# matrix of student engagement states over 15 weekly observations.

stopifnot(requireNamespace("Nestimate", quietly = TRUE))

engagement <- Nestimate::trajectories

stopifnot(
  is.matrix(engagement),
  nrow(engagement) == 138L,
  ncol(engagement) == 15L,
  all(unique(as.vector(engagement)) %in%
        c("Active", "Average", "Disengaged", NA))
)

usethis::use_data(engagement, overwrite = TRUE)
