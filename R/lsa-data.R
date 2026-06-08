# Canonical input form. Every entry point in lagseq accepts diverse
# user input shapes, but everything downstream reads only an `lsa_data`
# object. Centralizing input handling here means engines, bootstrap,
# permutation, and plotting all see the same simple structure.

#' Canonicalize Sequence Input for Lag Sequential Analysis
#'
#' Coerces a wide variety of user input shapes into a single canonical
#' representation used by every downstream lagseq function (engines,
#' bootstrap, permutation, grouping, plotting).
#'
#' Accepted input forms:
#'
#' - An atomic vector of integer or character codes — treated as a
#'   single sequence.
#' - A list of atomic vectors — treated as multiple independent
#'   sequences; transitions are not counted across sequence boundaries.
#' - A wide matrix or data.frame with rows = sequences, columns =
#'   ordered time points. Missing values (`NA`) and empty strings are
#'   treated as missingness, not as a state: they are dropped wherever
#'   they occur in a row and the surrounding events close up, so no
#'   transition is counted into or out of a gap. To model missingness
#'   as its own state, recode it (e.g. `NA -> "missing"`) before
#'   calling [lsa()].
#' - A square numeric matrix of pre-computed transition counts. Row
#'   `i`, column `j` is the count of `i -> j` transitions. In this case
#'   `events` and `seq_id` are not available and downstream resampling
#'   tools that need event-level data will error.
#' - A sequence-bearing object: a `tna` or
#'   `group_tna` (sequences read from its `$data` slot), a `tna_data`
#'   or `nestimate_data` (`$sequence_data`), or an `stslist`. The stored event
#'   sequences are recovered and analysed.
#'   A `tna` built from a bare matrix (no retained sequences) errors,
#'   because transition *counts* cannot be recovered from probability
#'   weights.
#'
#' @param x Sequence input. See Details.
#' @param labels Optional character vector of label names for the
#'   states. When `NULL`, labels are extracted from the data: unique
#'   sorted values of character input, or `"Code 1", "Code 2", ...` for
#'   integer input.
#'
#' @return An object of class `c("lsa_data", "list")` with elements:
#' \describe{
#'   \item{events}{Integer vector of event codes (1-indexed), or
#'     `NULL` if input was a transition matrix.}
#'   \item{seq_id}{Integer vector of sequence membership, same length
#'     as `events`, or `NULL` if input was a transition matrix.}
#'   \item{labels}{Character vector of state labels.}
#'   \item{n_states}{Number of distinct states (`K`).}
#'   \item{n_sequences}{Integer count of independent sequences.}
#'   \item{n_events}{Total number of events across all sequences.}
#'   \item{transitions_per_seq}{Integer vector: number of transitions
#'     each sequence contributes at lag 1.}
#'   \item{source}{One of `"events"`, `"transitions"` — flags whether
#'     event-level data is available.}
#'   \item{obs_input}{If `source = "transitions"`, the original `K x K`
#'     count matrix. Otherwise `NULL`.}
#' }
#'
#' @examples
#' # Single character sequence
#' d1 <- lsa_data(c("a", "b", "a", "c", "b"))
#' d1$n_events
#'
#' # Multiple sequences
#' d2 <- lsa_data(list(c("a", "b", "a"), c("b", "c", "a", "b")))
#' d2$n_sequences
#'
#' # Pre-computed transition matrix
#' tm <- matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3,
#'              dimnames = list(c("a","b","c"), c("a","b","c")))
#' d3 <- lsa_data(tm)
#' d3$source
#'
#' @seealso [lsa_transitions()], [lsa()]
#'
#' @export
lsa_data <- function(x, labels = NULL) {
  if (inherits(x, "lsa_data")) return(x)

  # Recover sequences from an external object (tna / group_tna,
  # nestimate_data, stslist, ...) before the generic shape
  # heuristics run. Without this, a fitted object -- which is just a
  # list -- would be silently misread as "a list of sequences".
  if (.is_seq_object(x)) {
    seqs <- .sequences_from_object(x)
    return(.lsa_data_from_sequences(seqs, labels = labels))
  }

  if (.is_transition_matrix(x)) {
    return(.lsa_data_from_matrix(x, labels = labels))
  }

  seqs <- .as_sequence_list(x)
  .lsa_data_from_sequences(seqs, labels = labels)
}

# Detect a square numeric matrix (possibly with row/colnames). A
# data.frame is never treated as a transition matrix even if square,
# because users routinely pass wide-format sequence data.frames.
.is_transition_matrix <- function(x) {
  is.matrix(x) && is.numeric(x) && nrow(x) == ncol(x) && nrow(x) >= 2L
}

# Coerce vector / list / wide matrix-or-df into a list of vectors.
# NA and empty-string cells are dropped wherever they occur, not just
# at the ends: they represent missingness, not a state, so nothing
# transitions into or out of them. To model missingness as a state,
# recode (e.g. NA -> "missing") before calling lsa().
.as_sequence_list <- function(x) {
  if (is.atomic(x) && is.null(dim(x))) {
    return(list(.clean_sequence(x)))
  }
  if (is.list(x) && !is.data.frame(x)) {
    out <- lapply(x, .clean_sequence)
    return(out[vapply(out, length, integer(1)) > 0L])
  }
  if (is.matrix(x) || is.data.frame(x)) {
    nr <- nrow(x)
    if (is.null(nr)) {
      stop("Unable to determine row count of input.", call. = FALSE)
    }
    out <- lapply(seq_len(nr), function(i) .clean_sequence(unlist(x[i, ])))
    return(out[vapply(out, length, integer(1)) > 0L])
  }
  stop(sprintf(
    "Unsupported input type: %s. Provide a vector, list of vectors, ",
    paste(class(x), collapse = "/")
  ), "wide matrix/data.frame, or square transition matrix.",
       call. = FALSE)
}

# Drop NA and empty-string cells anywhere in the vector. NA is treated
# as missingness, not as a state, so the observed-state stream is the
# raw vector with sentinels removed and the remaining order preserved.
.clean_sequence <- function(v) {
  v <- unname(v)
  if (is.character(v)) v <- v[!is.na(v) & nzchar(v)]
  else                 v <- v[!is.na(v)]
  v
}

# Build canonical lsa_data from a list of cleaned sequences.
.lsa_data_from_sequences <- function(seqs, labels = NULL) {
  if (length(seqs) == 0L) {
    stop("No usable sequences in input (all empty after NA removal).",
         call. = FALSE)
  }
  flat <- unlist(seqs)
  if (length(flat) == 0L) {
    stop("Input contains zero events.", call. = FALSE)
  }
  resolved <- .resolve_and_match(flat, labels)
  events_int <- resolved$events_int
  seq_id <- rep.int(seq_along(seqs),
                    times = vapply(seqs, length, integer(1)))
  per_seq_len <- vapply(seqs, length, integer(1))
  transitions_per_seq <- pmax(per_seq_len - 1L, 0L)
  structure(
    list(
      events = events_int,
      seq_id = seq_id,
      labels = resolved$labels,
      n_states = length(resolved$labels),
      n_sequences = length(seqs),
      n_events = length(events_int),
      transitions_per_seq = transitions_per_seq,
      source = "events",
      obs_input = NULL
    ),
    class = c("lsa_data", "list")
  )
}

# Build canonical lsa_data from a pre-computed transition matrix.
.lsa_data_from_matrix <- function(m, labels = NULL) {
  # Reject non-finite cells before the sign check: an NA would make
  # any(m < 0) itself NA ("missing value where TRUE/FALSE needed"), and
  # an Inf count would pass straight through into the engine and yield a
  # plausible-looking but invalid fit (NaN expecteds, all-NA residuals,
  # a tablewise test that spuriously reports perfect independence).
  if (anyNA(m) || !all(is.finite(m))) {
    stop("Transition matrix must contain only finite counts ",
         "(no NA, NaN, or Inf).", call. = FALSE)
  }
  if (any(m < 0)) {
    stop("Transition matrix must be non-negative.", call. = FALSE)
  }
  K <- nrow(m)
  dn <- dimnames(m)
  default_labels <- if (!is.null(dn) && !is.null(dn[[1]])) {
    dn[[1]]
  } else if (!is.null(dn) && !is.null(dn[[2]])) {
    dn[[2]]
  } else {
    paste("Code", seq_len(K))
  }
  if (is.null(labels)) {
    labels <- default_labels
  } else if (length(labels) != K) {
    stop(sprintf("labels has length %d but transition matrix is %d x %d.",
                 length(labels), K, K), call. = FALSE)
  }
  rownames(m) <- colnames(m) <- as.character(labels)
  structure(
    list(
      events = NULL,
      seq_id = NULL,
      labels = as.character(labels),
      n_states = K,
      n_sequences = NA_integer_,
      n_events = NA_integer_,
      transitions_per_seq = NA_integer_,
      source = "transitions",
      obs_input = m
    ),
    class = c("lsa_data", "list")
  )
}

# Coerce a numeric/logical event vector to 1-based integer codes, but
# only when every value is a whole number >= 1. Without the whole-number
# guard, as.integer() would silently truncate fractional codes
# (1.9 -> 1, 2.8 -> 2), relabelling the sequence into a different one
# than the user supplied.
.whole_event_codes <- function(flat) {
  if (anyNA(flat)) stop("Integer events contain NA.", call. = FALSE)
  if (!all(is.finite(flat))) {
    stop("Event codes must be finite (no NaN or Inf).", call. = FALSE)
  }
  if (!isTRUE(all(flat == floor(flat)))) {
    stop("Numeric event codes must be whole numbers; fractional values ",
         "such as 1.9 are not valid category codes. Pass character ",
         "labels instead if these are not integer codes.", call. = FALSE)
  }
  ints <- as.integer(flat)
  if (min(ints) < 1L) {
    stop("Integer event codes must be >= 1.", call. = FALSE)
  }
  ints
}

# Unified label resolution and integer matching. When `labels` is
# supplied, integers are treated as 1-based indices into `labels`;
# character/factor input is matched by name. When `labels` is NULL,
# derive a canonical set from the event vector itself.
.resolve_and_match <- function(flat, labels) {
  if (!is.null(labels)) {
    if (!is.character(labels)) labels <- as.character(labels)
    if (anyDuplicated(labels)) {
      stop("`labels` must be unique.", call. = FALSE)
    }
    if (is.numeric(flat) || is.integer(flat) || is.logical(flat)) {
      ints <- .whole_event_codes(flat)
      if (max(ints) > length(labels)) {
        stop(sprintf(
          "Maximum event code %d exceeds `labels` length %d.",
          max(ints), length(labels)
        ), call. = FALSE)
      }
      return(list(labels = labels, events_int = ints))
    }
    ints <- match(as.character(flat), labels)
    if (anyNA(ints)) {
      bad <- unique(as.character(flat)[is.na(ints)])
      stop(sprintf("Some events do not match supplied labels: %s",
                   paste(shQuote(bad), collapse = ", ")),
           call. = FALSE)
    }
    return(list(labels = labels, events_int = ints))
  }
  if (is.character(flat)) {
    lab <- sort(unique(flat))
    return(list(labels = lab, events_int = match(flat, lab)))
  }
  if (is.factor(flat)) {
    return(list(labels = levels(flat), events_int = as.integer(flat)))
  }
  if (is.numeric(flat) || is.integer(flat) || is.logical(flat)) {
    ints <- .whole_event_codes(flat)
    K <- max(ints)
    return(list(labels = paste("Code", seq_len(K)), events_int = ints))
  }
  stop(sprintf("Unsupported event vector type: %s",
               paste(class(flat), collapse = "/")), call. = FALSE)
}

#' @export
print.lsa_data <- function(x, ...) {
  cat("<lsa_data>\n")
  cat(sprintf("  source:        %s\n", x$source))
  cat(sprintf("  states (K):    %d\n", x$n_states))
  cat(sprintf("  labels:        %s\n",
              paste(utils::head(x$labels, 10), collapse = ", ")))
  if (identical(x$source, "events")) {
    cat(sprintf("  sequences:     %d\n", x$n_sequences))
    cat(sprintf("  events:        %d\n", x$n_events))
    cat(sprintf("  transitions:   %d (lag 1, summed across sequences)\n",
                sum(x$transitions_per_seq)))
  } else {
    cat(sprintf("  transitions:   %d (from input matrix)\n",
                as.integer(sum(x$obs_input))))
  }
  invisible(x)
}
