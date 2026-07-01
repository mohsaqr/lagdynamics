# Provenance for data/group_regulation_long.rda
#
# Derived without modification from tna::group_regulation_long (the `tna`
# package, MIT license), stored as a plain data.frame so the package
# carries no tibble/tidyverse coupling. This script is kept for
# reproducibility only; `tna` is NOT a dependency of lagdynamics.
#
# Run from the package root with tna installed:
#   source("data-raw/group_regulation_long.R")

data("group_regulation_long", package = "tna", envir = environment())
group_regulation_long <- as.data.frame(group_regulation_long,
                                       stringsAsFactors = FALSE)
usethis::use_data(group_regulation_long, overwrite = TRUE, compress = "xz")
