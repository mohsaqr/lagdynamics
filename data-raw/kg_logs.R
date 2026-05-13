# Build lagseq's `kg_logs` and `kg_lsa_oracle` data objects from the
# Du Jun (2026) Mendeley dataset
# "Sequential Behavioral Mechanisms Linking Learning Paths to Academic
# Performance in Knowledge Graph Environments"
# DOI: 10.17632/bdwcj7vw94.1  -- License: CC BY 4.0
#
# Run with: Rscript data-raw/kg_logs.R
#
# Two objects are produced and saved to data/:
#   kg_logs        — list of 29 character vectors, one event sequence
#                    per learner, taken from the wide-format sheet
#                    ("整体数据"), which is the sheet the authors used
#                    to compute the published LSA matrices.
#   kg_lsa_oracle  — list of 4 published LSA outputs (overall, low,
#                    mid, high) each containing $obs (the published
#                    transition frequency matrix JNTF) and $adj_res
#                    (the published adjusted-residual matrix ADJR).
#                    These are the *authors' own* LSA results, used as
#                    an independent third-party oracle.

stopifnot(requireNamespace("readxl", quietly = TRUE))
suppressMessages(library(readxl))

xlsx_path <- "data-raw/kg_logs.xlsx"
LABELS <- c(paste0("L0", 1:9), paste0("L", 10:12), "E")

# --- 1. Event sequences from the wide-format sheet -----------------------

suppressWarnings({
  wide <- read_excel(xlsx_path, sheet = "整体数据",
                     col_types = "text", col_names = FALSE)
})
wide <- wide[-1, , drop = FALSE]                # drop header row
ids    <- wide[[1]]
groups <- wide[[2]]                              # 低/中/高
event_cols <- wide[, -c(1L, 2L), drop = FALSE]

kg_logs <- lapply(seq_len(nrow(event_cols)), function(i) {
  v <- unlist(event_cols[i, ])
  v <- v[!is.na(v) & nzchar(v)]
  unname(v)
})
names(kg_logs) <- ids

# Group attribute for downstream group_lsa() testing.
attr(kg_logs, "group") <- setNames(groups, ids)

# --- 2. Published LSA oracle matrices ------------------------------------

read_jntf_adjr <- function(sheet, label_set = LABELS) {
  suppressWarnings({
    s <- read_excel(xlsx_path, sheet = sheet,
                    col_types = "text", col_names = FALSE)
  })
  K <- length(label_set)
  obs     <- matrix(as.numeric(unlist(s[3:(2 + K),       2:(1 + K)])),
                    K, K, dimnames = list(label_set, label_set))
  adj_res <- matrix(as.numeric(unlist(s[20:(19 + K),     2:(1 + K)])),
                    K, K, dimnames = list(label_set, label_set))
  list(obs = obs, adj_res = adj_res)
}

kg_lsa_oracle <- list(
  overall = read_jntf_adjr("整体分析"),
  low     = read_jntf_adjr("低结果"),
  mid     = read_jntf_adjr("中结果"),
  high    = read_jntf_adjr("高结果")
)

# --- 3. Sanity ------------------------------------------------------------

stopifnot(
  length(kg_logs) == 29L,
  all(unique(unlist(kg_logs)) %in% LABELS),
  identical(rownames(kg_lsa_oracle$overall$obs), LABELS),
  sum(kg_lsa_oracle$overall$obs) == 870
)

# --- 4. Save --------------------------------------------------------------

save(kg_logs,        file = "data/kg_logs.rda",       compress = "bzip2")
save(kg_lsa_oracle,  file = "data/kg_lsa_oracle.rda", compress = "bzip2")

cat("data/kg_logs.rda:        ", length(kg_logs), "sequences,",
    sum(vapply(kg_logs, length, integer(1))), "events\n")
cat("data/kg_lsa_oracle.rda:  4 oracle matrices",
    "(overall + low + mid + high), each 13 x 13\n")
