# Build the oconnor_couple dataset: the canonical LSA worked example
# from O'Connor (1999), "Simple and flexible SAS and SPSS programs
# for analyzing lag-sequential categorical data", Behavior Research
# Methods, Instruments, & Computers, 31(4), 718-726.
# doi:10.3758/BF03200753
#
# Appendix A publishes a 393-event input sequence (couple-interaction
# data from Gottman & Roy 1990, p. 78). Appendix B publishes the full
# numerical output of the SEQUENTIAL program for this exact input.
# We ship both, so lagseq's clean-room reimplementation can be tested
# against the canonical published outputs cell by cell.
#
# Run with: Rscript data-raw/oconnor_couple.R

# --- Appendix A: input sequence (verbatim from p. 723) ----------------------
# 12 rows printed; concatenated into one length-393 integer sequence.
sequence <- as.integer(c(
  2,5,2,4,2,5,2,4,2,4,2,2,5,2,4,2,5,3,4,3,4,1,4,4,1,2,5,2,1,5,4,2,5,
  1,4,4,2,3,6,6,2,2,6,5,5,5,6,3,6,3,3,6,3,2,5,3,5,2,6,2,2,2,5,2,4,4,
  5,2,5,2,5,5,4,2,5,5,5,2,5,2,5,2,5,2,2,5,2,5,2,2,5,5,2,2,2,5,2,5,2,
  5,2,5,2,2,5,2,5,3,5,2,5,5,2,2,5,2,2,5,2,6,2,5,4,2,5,2,4,4,2,1,2,5,
  2,5,2,5,5,5,3,5,2,5,5,2,5,5,2,5,2,5,2,4,2,4,2,5,4,2,2,5,2,5,5,5,2,
  5,2,5,5,2,5,2,5,5,2,1,2,5,1,5,1,5,5,4,2,2,2,3,6,3,6,3,6,3,6,3,6,3,
  3,6,6,6,5,1,2,5,2,5,5,2,4,5,5,2,5,5,2,5,2,5,2,5,1,5,4,2,2,5,2,5,2,
  5,4,1,4,4,2,4,4,2,4,2,3,5,4,2,2,5,2,6,1,4,1,4,5,5,5,4,5,2,4,5,5,2,
  5,2,2,2,5,2,2,5,2,2,4,2,2,2,6,2,4,3,3,3,3,2,3,6,3,5,3,3,5,3,2,2,2,
  2,6,3,6,3,6,2,6,1,6,3,3,3,3,5,3,5,2,5,3,3,6,3,1,4,1,5,1,5,6,3,3,6,
  3,6,3,6,3,6,3,6,3,3,6,3,6,3,2,5,2,5,2,5,2,5,2,2,2,3,5,1,6,1,6,6,6,
  1,6,6,1,6,1,6,2,5,1,6,6,6,2,2,6,6,6,2,6,2,6,2,5,5,5,5,5,5,6
))
stopifnot(length(sequence) == 393L,
          all(sequence %in% 1:6))

LABS <- paste("Code", 1:6)

# --- Appendix B: published output matrices ----------------------------------

# Cell frequencies (Table on p. 724)
obs <- matrix(c(
  0,  4,  0,  6,  6,  7,
  3, 28,  5, 14, 63, 10,
  1,  4, 12,  2,  9, 20,
  5, 19,  2,  6,  5,  0,
  8, 57,  8,  9, 27,  3,
  6, 10, 21,  0,  2, 10
), 6, 6, byrow = TRUE, dimnames = list(LABS, LABS))
stopifnot(sum(obs) == 392L)

# Expected frequencies
expected <- matrix(c(
  1.3495,  7.1582,  2.8163,  2.1709,  6.5714,  2.9337,
  7.2168, 38.281,  15.061,  11.610,  35.143,  15.689,
  2.8163, 14.939,   5.8776,  4.5306, 13.714,   6.1224,
  2.1709, 11.515,   4.5306,  3.4923, 10.571,   4.7194,
  6.5714, 34.857,  13.714,  10.571,  32.000,  14.286,
  2.8750, 15.250,   6.0000,  4.6250, 14.000,   6.2500
), 6, 6, byrow = TRUE, dimnames = list(LABS, LABS))

# Transitional probabilities
prob <- matrix(c(
  0.0000, 0.1739, 0.0000, 0.2609, 0.2609, 0.3043,
  0.0244, 0.2276, 0.0407, 0.1138, 0.5122, 0.0813,
  0.0208, 0.0833, 0.2500, 0.0417, 0.1875, 0.4167,
  0.1351, 0.5135, 0.0541, 0.1622, 0.1351, 0.0000,
  0.0714, 0.5089, 0.0714, 0.0804, 0.2411, 0.0268,
  0.1224, 0.2041, 0.4286, 0.0000, 0.0408, 0.2041
), 6, 6, byrow = TRUE, dimnames = list(LABS, LABS))

# Tablewise likelihood-ratio chi-square test
lrx2 <- list(statistic = 202.5009, df = 25L, p = 0.0000)

# Adjusted residuals
adj_res <- matrix(c(
  -1.234, -1.466, -1.846,  2.815, -0.272,  2.620,
  -1.953, -2.417, -3.341,  0.890,  6.712, -1.856,
  -1.191, -3.640,  2.878, -1.334, -1.608,  6.410,
   2.080,  2.793, -1.334,  1.482, -2.131, -2.444,
   0.680,  5.347, -1.949, -0.601, -1.237, -3.782,
   2.031, -1.732,  6.988, -2.416, -4.057,  1.717
), 6, 6, byrow = TRUE, dimnames = list(LABS, LABS))

# P-values for adjusted residuals
adj_p <- matrix(c(
  0.2172, 0.1427, 0.0648, 0.0049, 0.7857, 0.0088,
  0.0508, 0.0157, 0.0008, 0.3735, 0.0000, 0.0634,
  0.2337, 0.0003, 0.0040, 0.1823, 0.1079, 0.0000,
  0.0376, 0.0052, 0.1823, 0.1384, 0.0331, 0.0145,
  0.4967, 0.0000, 0.0513, 0.5479, 0.2159, 0.0002,
  0.0423, 0.0833, 0.0000, 0.0157, 0.0000, 0.0860
), 6, 6, byrow = TRUE, dimnames = list(LABS, LABS))

# Yule's Q values
yules_q <- matrix(c(
  -1.000, -0.381, -1.000,  0.587, -0.066,  0.537,
  -0.525, -0.291, -0.636,  0.157,  0.650, -0.327,
  -0.525, -0.703,  0.481, -0.445, -0.299,  0.764,
   0.490,  0.442, -0.445,  0.338, -0.468, -1.000,
   0.152,  0.548, -0.368, -0.120, -0.157, -0.760,
   0.456, -0.308,  0.795, -1.000, -0.835,  0.320
), 6, 6, byrow = TRUE, dimnames = list(LABS, LABS))

# Unidirectional kappas (Wampold's transformed kappa)
kappa <- matrix(c(
  -1.0000, -0.4443, -1.0000,  0.1840, -0.0846,  0.2029,
  -0.5832, -0.2727, -0.6672,  0.0952,  0.3632, -0.3610,
  -0.6440, -0.7337,  0.1457, -0.5574, -0.3421,  0.3316,
   0.1361,  0.2919, -0.5574,  0.0751, -0.5258, -1.0000,
   0.0879,  0.2852, -0.4152, -0.1465, -0.1541, -0.7895,
   0.1531, -0.3610,  0.3555, -1.0000, -0.8596,  0.0834
), 6, 6, byrow = TRUE, dimnames = list(LABS, LABS))

# z values for kappas
kappa_z <- matrix(c(
  -1.231, -1.480, -1.841,  2.818, -0.264,  2.624,
  -1.943, -2.459, -3.326,  0.900,  6.726, -1.842,
  -1.186, -3.657,  2.884, -1.327, -1.595,  6.415,
   2.083,  2.760, -1.327,  1.487, -2.119, -2.437,
   0.687,  5.282, -1.936, -0.590, -1.216, -3.768,
   1.980, -1.842,  6.876, -2.437, -4.103,  1.651
), 6, 6, byrow = TRUE, dimnames = list(LABS, LABS))

# P-values for kappas
kappa_p <- matrix(c(
  0.2184, 0.1388, 0.0656, 0.0048, 0.7920, 0.0087,
  0.0520, 0.0139, 0.0009, 0.3680, 0.0000, 0.0655,
  0.2357, 0.0003, 0.0039, 0.1845, 0.1107, 0.0000,
  0.0372, 0.0058, 0.1845, 0.1371, 0.0341, 0.0148,
  0.4920, 0.0000, 0.0529, 0.5550, 0.2240, 0.0002,
  0.0477, 0.0655, 0.0000, 0.0148, 0.0000, 0.0988
), 6, 6, byrow = TRUE, dimnames = list(LABS, LABS))

# Permutation test results (saved for Step 5 validation of permute_lsa)
# Mean p-values from 10 blocks x 1000 permutations
perm_p_mean <- matrix(c(
  0.37890, 0.16610, 0.09050, 0.01400, 0.81660, 0.01410,
  0.05970, 0.01850, 0.00060, 0.46340, 0.00000, 0.07570,
  0.33300, 0.00090, 0.00580, 0.20660, 0.12750, 0.00000,
  0.05570, 0.00790, 0.20750, 0.13100, 0.03490, 0.01670,
  0.63200, 0.00000, 0.06550, 0.57750, 0.26500, 0.00010,
  0.05810, 0.10330, 0.00000, 0.01360, 0.00020, 0.10380
), 6, 6, byrow = TRUE, dimnames = list(LABS, LABS))
perm_p_high <- matrix(c(
  0.38868, 0.17553, 0.09550, 0.01644, 0.82396, 0.01566,
  0.06533, 0.02028, 0.00112, 0.47211, 0.00000, 0.08132,
  0.34134, 0.00144, 0.00747, 0.21125, 0.13420, 0.00000,
  0.05923, 0.00989, 0.21426, 0.13874, 0.03992, 0.01843,
  0.64930, 0.00000, 0.07108, 0.58828, 0.27164, 0.00030,
  0.06310, 0.10766, 0.00000, 0.01585, 0.00046, 0.11023
), 6, 6, byrow = TRUE, dimnames = list(LABS, LABS))
perm_p_low <- matrix(c(
  0.36912, 0.15667, 0.08550, 0.01156, 0.80924, 0.01254,
  0.05407, 0.01672, 0.00008, 0.45469, 0.00000, 0.07008,
  0.32466, 0.00036, 0.00413, 0.20195, 0.12080, 0.00000,
  0.05217, 0.00591, 0.20074, 0.12326, 0.02988, 0.01497,
  0.61470, 0.00000, 0.05992, 0.56672, 0.25836, -0.00010,
  0.05310, 0.09894, 0.00000, 0.01135, -0.00006, 0.09737
), 6, 6, byrow = TRUE, dimnames = list(LABS, LABS))

# --- Final object ----------------------------------------------------------
oconnor_couple <- list(
  sequence    = sequence,
  obs         = obs,
  expected    = expected,
  prob        = prob,
  lrx2        = lrx2,
  adj_res     = adj_res,
  adj_p       = adj_p,
  yules_q     = yules_q,
  kappa       = kappa,
  kappa_z     = kappa_z,
  kappa_p     = kappa_p,
  permutation = list(
    n_blocks     = 10L,
    n_perm_block = 1000L,
    p_mean       = perm_p_mean,
    p_high       = perm_p_high,
    p_low        = perm_p_low
  ),
  source = paste0(
    "O'Connor, B. P. (1999). Simple and flexible SAS and SPSS ",
    "programs for analyzing lag-sequential categorical data. ",
    "Behavior Research Methods, Instruments, & Computers, 31(4), ",
    "718-726. doi:10.3758/BF03200753. Input data originally from ",
    "Gottman, J. M., & Roy, A. K. (1990), p. 78."
  ),
  notes = paste(
    "Canonical LSA worked example. Input sequence is 393 events long",
    "(starting with code 2, ending with code 6); 392 transitions at",
    "lag 1. Output matrices transcribed verbatim from Appendix B."
  )
)

save(oconnor_couple, file = "data/oconnor_couple.rda",
     compress = "bzip2")
cat("data/oconnor_couple.rda saved (",
    file.size("data/oconnor_couple.rda"), "bytes )\n")
cat("Per-code event counts:\n")
print(table(sequence))
