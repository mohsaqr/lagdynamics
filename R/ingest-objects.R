# Sequence recovery from external objects. lagseq is the sequence
# *analyser*; whenever an object that already carries event sequences is
# passed in, we extract those sequences and feed them through the normal
# lsa_data() pipeline. This is an input convenience and does not couple
# lagseq to the source package: every extractor reads plain
# list/matrix/attribute structure, never the source package's functions.
#
# Supported sources:
#   tna / group_tna          -> $data (a tna_seq_data code matrix)
#   tna_seq_data             -> the code matrix directly
#   tna_data                 -> $sequence_data (wide label df)
#   nestimate_data           -> $sequence_data (wide label df)
#   stslist                  -> factor state matrix + alphabet
#
# A model carrying the `tna` class is handled by the `tna` branch with
# no extra code. The `tna_data` and `nestimate_data` containers are
# structurally identical (both carry $sequence_data), so they share one
# branch.

# Does `x` carry recoverable event sequences from a known source?
.is_seq_object <- function(x) {
  inherits(x, c("tna", "group_tna", "tna_seq_data",
                "tna_data", "nestimate_data", "stslist"))
}

# Dispatch on the object class and return a list of event sequences
# (character or integer vectors, one per sequence).
.sequences_from_object <- function(x) {
  if (inherits(x, "group_tna")) {
    # Pool the per-group sequences into one set; grouping is dropped
    # because lsa() builds a single fit. (Use lsa(..., group=) on the
    # recovered data to re-group.)
    parts <- lapply(unclass(x), .seqs_from_tna)
    return(unlist(parts, recursive = FALSE, use.names = FALSE))
  }
  if (inherits(x, "tna"))           return(.seqs_from_tna(x))
  if (inherits(x, "tna_seq_data"))  return(.seqs_from_seqdata_matrix(x))
  if (inherits(x, c("tna_data", "nestimate_data"))) {
    # Both carry a wide one-row-per-session label table in
    # $sequence_data.
    sd <- x$sequence_data
    if (is.null(sd)) {
      stop("This ", class(x)[1L], " object has no $sequence_data.",
           call. = FALSE)
    }
    return(.as_sequence_list(sd))
  }
  if (inherits(x, "stslist"))       return(.seqs_from_stslist(x))
  stop("Unsupported sequence object: ",
       paste(class(x), collapse = "/"), call. = FALSE)
}

# A tna model carries its source sequences in $data only when it was
# built from sequence data. Models built from a bare weight/transition
# matrix have $data = NULL and weights that are row-normalised
# probabilities, from which transition *counts* cannot be faithfully
# recovered -- so we refuse rather than fabricate.
.seqs_from_tna <- function(x) {
  d <- x$data
  if (is.null(d)) {
    stop("This tna object carries no sequence data (it was built from ",
         "a matrix, so its weights are probabilities, not counts). ",
         "lsa() needs the original sequences: rebuild the tna from ",
         "sequence data, or pass the sequences to lsa() directly.",
         call. = FALSE)
  }
  .seqs_from_seqdata_matrix(d)
}

# Decode a tna_seq_data code matrix (rows = sequences, cols = time,
# integer codes, NA-padded) into a list of label sequences using its
# `alphabet` attribute. Falls back to integer codes when no alphabet
# is attached (lsa_data() then derives labels itself).
.seqs_from_seqdata_matrix <- function(d) {
  alphabet <- attr(d, "alphabet")
  codes <- matrix(as.integer(d), nrow = nrow(d))
  rows <- lapply(seq_len(nrow(codes)), function(i) {
    v <- codes[i, ]
    v <- v[!is.na(v)]
    if (!is.null(alphabet)) alphabet[v] else v
  })
  rows[vapply(rows, length, integer(1L)) > 0L]
}

# Decode an stslist into a list of label sequences. Cells are
# factor-coded states; the void/NA sentinels mark gaps, which we drop
# (consistent with lagseq's NA-as-missingness convention).
.seqs_from_stslist <- function(x) {
  void <- attr(x, "void")
  m <- as.matrix(as.data.frame(x, stringsAsFactors = FALSE))
  rows <- lapply(seq_len(nrow(m)), function(i) {
    v <- as.character(m[i, ])
    keep <- !is.na(v) & nzchar(v)
    if (!is.null(void)) keep <- keep & v != as.character(void)
    v[keep]
  })
  rows[vapply(rows, length, integer(1L)) > 0L]
}
