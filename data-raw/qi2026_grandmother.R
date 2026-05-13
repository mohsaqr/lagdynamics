# Build the qi2026_grandmother dataset: the published transition
# frequency matrix (Table 4), Z-score matrix and Yule's Q matrix
# (Table 5) for the grandmother in Qi An et al. (2026),
# "Behavioural Trajectories and Spatial Responses",
# Sustainability 18(5), 2326. doi:10.3390/su18052326. Open Access.
#
# All values transcribed verbatim from the published tables. Run with:
#   Rscript data-raw/qi2026_grandmother.R

LABS <- c("CO","DH","WO","CK","SM","EA","CM","ED","RE","UO")

# --- Table 4: transition frequencies (priming rows x lag columns) ---
obs <- matrix(c(
  344,  44,  10,  56,  48,  39,  70,  0,  17,  0,
   62,  55,   4,   8,  18,   2,  27,  0,   9,  0,
   14,   1,   8,   4,   2,   0,   2,  0,   1,  0,
   34,   8,   4,  37,   1,   4,  11,  0,   1,  1,
   68,  16,   4,  13,  12,   3,  20,  0,   0,  0,
   32,   5,   1,   4,   1,   1,   6,  0,   4,  1,
  117,  27,  18,  31,  10,   2,  47,  0,  18,  0,
   14,   0,   8,   3,   0,   0,   0,  0,   0,  0,
   31,   8,  11,   9,   2,   2,  16,  0,  11,  0,
    3,   2,   1,   0,   2,   0,   1,  0,   0,  0
), 10, 10, byrow = TRUE, dimnames = list(LABS, LABS))

stopifnot(sum(obs) == 1531L)

# --- Table 5: Z-scores (published) ---
adj_res <- matrix(c(
   5.09,   4.033, -4.588, -1.965,  1.842,  4.901, -1.864, 0, -2.135, -0.835,
  -3.918,  8.807, -1.641, -3.021,  2.067, -1.89,   0.655, 0,  0.651, -0.371,
  -0.372, -1.42,   5.645,  0.316, -0.006, -1.083, -1.157, 0, -0.252, -0.146,
  -2.777, -0.979, -0.275,  8.666, -2.266,  0.282, -0.673, 0, -1.593,  3.763,
   0.736,  0.359, -0.924, -0.483,  1.284, -0.841,  0.592, 0, -2.49,  -0.312,
   1.839, -0.383, -0.958, -0.815, -1.364, -0.66,  -0.435, 0,  1.308, -0.191,
  -1.328, -0.495,  1.882,  0.407, -1.92,  -2.697,  2.329, 0,  2.48,  -0.463,
   0.91,  -1.759,  1.678,  0.198, -1.304, -0.955, -1.955, 0, -1.027, -0.129,
  -2.459, -0.617,  3.634, -0.247, -1.634, -0.664,  1.365, 0,  4.116, -0.25,
  -0.824,  1.1,    0.957, -1.046,  1.979, -0.57,  -0.175, 0, -0.613, -0.077
), 10, 10, byrow = TRUE, dimnames = list(LABS, LABS))

# --- Table 5: Yule's Q values (published) ---
yules_q <- matrix(c(
   0.26,  -0.35,  -0.624, -0.168,  0.191,  0.615, -0.146, 0, -0.297, -1,
  -0.309,  0.649, -0.394, -0.49,   0.273, -0.566,  0.073, 0,  0.12,  -1,
  -0.067, -0.587,  0.774,  0.085, -0.002, -1,     -0.391, 0, -0.128, -1,
  -0.289, -0.182, -0.072,  0.709, -0.754,  0.075, -0.11,  0, -0.628,  1,
   0.066,  0.05,  -0.235, -0.073,  0.203, -0.245,  0.075, 0, -1,     -1,
   0.251, -0.091, -0.438, -0.21,  -0.57,  -0.319, -0.095, 0,  0.331, -1,
  -0.089, -0.055,  0.257,  0.043, -0.311, -0.699,  0.208, 0,  0.338, -1,
   0.182, -1,      0.235,  0.061, -1,     -1,     -1,     0, -1,     -1,
  -0.27,  -0.116,  0.537, -0.045, -0.509, -0.235,  0.192, 0,  0.589, -1,
  -0.28,   0.405,  0.455, -1,     0.625,  -1,     -0.093, 0, -1,     -1
), 10, 10, byrow = TRUE, dimnames = list(LABS, LABS))

# --- Catalogued paper typos -----------------------------------------------
# Cells where Table 5 prints a value that does not match the math
# computed from Table 4. The "math_computed" column is what
# stats::chisq.test(obs, correct=FALSE)$stdres yields.
known_typos <- data.frame(
  from = c("CO",    "ED",    "EA",    "CK"),
  to   = c("DH",    "WO",    "UO",    "UO"),
  paper_printed = c( 4.033,   1.678,  -0.191,   3.763),
  math_computed = c(-4.026,   6.681,   3.529,   2.474),
  category = c("sign error", "transcription error",
               "transcription error", "transcription error"),
  stringsAsFactors = FALSE
)

# --- Code labels ----------------------------------------------------------
code_descriptions <- c(
  CO  = "Cooking",
  DH  = "Doing housework",
  WO  = "Working",
  CK  = "Caring for kid",
  SM  = "Self-management",
  EA  = "Eating",
  CM  = "Communicating",
  ED  = "Education",
  RE  = "Resting",
  UO  = "Using object"
)

# --- Final object ---------------------------------------------------------
qi2026_grandmother <- list(
  obs               = obs,
  adj_res           = adj_res,
  yules_q           = yules_q,
  known_typos       = known_typos,
  code_descriptions = code_descriptions,
  source            = "Qi An, W. Xing, Y. Wang, X. Li (2026). Sustainability, 18(5), 2326. doi:10.3390/su18052326",
  license           = "Open Access (paper); numerical tables are facts and not copyrightable",
  n_transitions     = 1531L,
  k_states          = 10L,
  notes             = paste(
    "Grandmother behaviour transitions at lag = 1 across 14 days in two",
    "Beijing dual-income households. Published Z-scores in $adj_res",
    "contain 4 documented typos (see $known_typos) where the printed",
    "values do not match the math computed from $obs."
  )
)

save(qi2026_grandmother,
     file = "data/qi2026_grandmother.rda", compress = "bzip2")
cat("data/qi2026_grandmother.rda saved (",
    file.size("data/qi2026_grandmother.rda"), "bytes )\n")
