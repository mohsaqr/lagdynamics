# Canonical input form. Every entry point in lagdynamics accepts diverse
# user input shapes, but everything downstream reads only an `lsa_data`
# object. Centralizing input handling here means engines, bootstrap,
# permutation, and plotting all see the same simple structure.

#' Convert sequence input to the canonical form used by lag sequential analysis
#'
#' @description
#' `lsa_data()` converts the forms in which sequence data usually arrive into
#' the canonical representation used by lag sequential analysis, and it is what
#' [lsa()] calls on its input. A lag sequential analysis needs, for every
#' sequence, the ordered list of coded events, and a common set of state labels
#' across sequences. Transitions are counted within sequences only, so the
#' boundaries between sequences must be kept.
#'
#' Accepted forms are:
#' * an atomic vector of codes (one sequence);
#' * a list of vectors (several sequences);
#' * a wide matrix or data frame, one row per sequence and one column per
#'   ordered time point;
#' * a square numeric matrix of transition counts, row `i` and column `j`
#'   holding the count of `i -> j` (any square numeric matrix is read this
#'   way, so pass square integer-coded sequences as a data frame);
#' * a sequence object: a `tna` or `group_tna` model with its data, a
#'   `tna_seq_data`, `tna_data` or `nestimate_data` object, a
#'   \pkg{Nestimate} `netobject`, or a \pkg{TraMineR} `stslist`.
#'
#' `NA` and empty strings are treated as missing observations, not as a
#' state: they are dropped wherever they occur and the surrounding events
#' close up. To model missingness as a state, recode it (for example to
#' `"missing"`) first. A transition-count matrix carries no sequences, so
#' the functions that resample sequences ([bootstrap_lsa()],
#' [permute_lsa()], [stability_lsa()], [reliability_lsa()],
#' [compare_lsa()]) cannot be used on a fit built from it. Long event logs
#' are read by [lsa()] through its `actor` and `action` arguments.
#'
#' @param x Sequence input in one of the forms listed above. An `lsa_data`
#'   object is returned unchanged.
#' @param labels Optional character vector of state names. Default `NULL`
#'   takes the sorted distinct values of character input, the levels of a
#'   factor, or `"Code 1"`, `"Code 2"`, ... for integer codes. When
#'   supplied, character codes are matched to it by name and integer codes
#'   index into it.
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and `lagdynamics_bad_input`:
#' an unsupported input type, no usable events, invalid counts in a
#' transition matrix, invalid integer codes, or codes that do not match
#' `labels`.
#'
#' @return An object of class `c("lsa_data", "list")`. Print it for the
#'   input source, the number of states, sequences, events and lag-1
#'   transitions; read the events with `as.data.frame()` (one row per
#'   event, or per cell for a transition matrix). The object stores the
#'   integer-coded events and their sequence membership (`NULL` for a
#'   transition matrix), the state labels, the counts of states,
#'   sequences and events, and the source (`"events"` or `"transitions"`).
#'
#' @examples
#' lsa_data(engagement)
#'
#' # Several sequences: no transition is counted across the boundary
#' lsa_data(list(c("a", "b", "a"), c("b", "c", "a", "b")))
#'
#' # A transition-count matrix
#' counts <- matrix(c(0, 3, 1, 2, 0, 4, 5, 1, 0), 3, 3,
#'                  dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
#' lsa_data(counts)
#'
#' @seealso [as.data.frame.lsa_data()] for the event table;
#'   [lsa_transitions()] for transition counts; [lsa()].
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
      .lsa_stop("Unable to determine row count of input.", call. = FALSE)
    }
    out <- lapply(seq_len(nr), function(i) .clean_sequence(unlist(x[i, ])))
    return(out[vapply(out, length, integer(1)) > 0L])
  }
  .lsa_stop(sprintf(
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
    .lsa_stop("No usable sequences in input (all empty after NA removal).",
         call. = FALSE)
  }
  flat <- unlist(seqs)
  if (length(flat) == 0L) {
    .lsa_stop("Input contains zero events.", call. = FALSE)
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
    .lsa_stop("Transition matrix must contain only finite counts ",
         "(no NA, NaN, or Inf).", call. = FALSE)
  }
  if (any(m < 0)) {
    .lsa_stop("Transition matrix must be non-negative.", call. = FALSE)
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
    .lsa_stop(sprintf("labels has length %d but transition matrix is %d x %d.",
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
  if (anyNA(flat)) .lsa_stop("Integer events contain NA.", call. = FALSE)
  if (!all(is.finite(flat))) {
    .lsa_stop("Event codes must be finite (no NaN or Inf).", call. = FALSE)
  }
  if (!isTRUE(all(flat == floor(flat)))) {
    .lsa_stop("Numeric event codes must be whole numbers; fractional values ",
         "such as 1.9 are not valid category codes. Pass character ",
         "labels instead if these are not integer codes.", call. = FALSE)
  }
  ints <- as.integer(flat)
  if (min(ints) < 1L) {
    .lsa_stop("Integer event codes must be >= 1.", call. = FALSE)
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
      .lsa_stop("`labels` must be unique.", call. = FALSE)
    }
    if (is.numeric(flat) || is.integer(flat) || is.logical(flat)) {
      ints <- .whole_event_codes(flat)
      if (max(ints) > length(labels)) {
        .lsa_stop(sprintf(
          "Maximum event code %d exceeds `labels` length %d.",
          max(ints), length(labels)
        ), call. = FALSE)
      }
      return(list(labels = labels, events_int = ints))
    }
    ints <- match(as.character(flat), labels)
    if (anyNA(ints)) {
      bad <- unique(as.character(flat)[is.na(ints)])
      .lsa_stop(sprintf("Some events do not match supplied labels: %s",
                   paste(sprintf("'%s'", bad), collapse = ", ")),
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
  .lsa_stop(sprintf("Unsupported event vector type: %s",
               paste(class(flat), collapse = "/")), call. = FALSE)
}

#' Read the events of canonical sequence data as a table
#'
#' @description
#' `as.data.frame()` returns the canonical sequence data from [lsa_data()],
#' every event with its sequence and position, as a tidy table, which shows how
#' the input was read: which rows became sequences, where missing values were
#' dropped, and which labels were assigned. For input given as a
#' transition-count matrix, it returns one row per cell instead. `print()` shows
#' the input source, the number of states, sequences, events and lag-1
#' transitions.
#'
#' @param x An `lsa_data` object from [lsa_data()].
#' @param row.names,optional,... Ignored; present for consistency with the
#'   generic.
#'
#' @return `as.data.frame()`: a `data.frame` with one row per event and
#'   the columns `seq_id`, `index` (position within the sequence) and
#'   `state`; or, for transition-matrix input, one row per cell and the
#'   columns `from`, `to` and `count`. `print()` returns `x` invisibly.
#'
#' @examples
#' d <- lsa_data(list(c("a", "b", NA, "a"), c("b", "c")))
#' d
#' as.data.frame(d)
#'
#' @seealso [lsa_data()], [lsa_transitions()]
#'
#' @export
as.data.frame.lsa_data <- function(x, row.names = NULL, optional = FALSE,
                                   ...) {
  if (identical(x$source, "events")) {
    index <- stats::ave(seq_along(x$events), x$seq_id,
                        FUN = seq_along)
    data.frame(seq_id = x$seq_id, index = index,
               state = x$labels[x$events],
               stringsAsFactors = FALSE, row.names = NULL)
  } else {
    m <- x$obs_input
    grid <- expand.grid(from = rownames(m), to = colnames(m),
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    data.frame(from = grid$from, to = grid$to, count = as.vector(m),
               stringsAsFactors = FALSE, row.names = NULL)
  }
}

#' @rdname as.data.frame.lsa_data
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
