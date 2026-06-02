# Long-format event-log sequencing. Raw interaction logs arrive with
# one row per event and columns identifying who acted (actor), what
# they did (action), and when (time / order). lagseq's analysis works
# on already-grouped sequences, so this base-R helper turns a long log
# into a list of event sequences: group by actor (and optionally an
# explicit session column), order within each group by time/order, and
# -- when a `time` column is given -- start a new session whenever the
# gap between consecutive events exceeds `time_threshold` seconds.
#
# The argument names and the 900-second default match the sibling
# ecosystem (tna::prepare_data, Nestimate::prepare) for familiarity,
# but the implementation is independent base R with no for-loops.

# Build a list of action sequences from a long-format data.frame.
.prepare_long <- function(data, actor, action, time = NULL, order = NULL,
                          session = NULL, time_threshold = 900,
                          custom_format = NULL, is_unix_time = FALSE,
                          unix_time_unit = "seconds") {
  if (!is.data.frame(data)) {
    stop("Long-format sequencing needs a data.frame; got ",
         paste(class(data), collapse = "/"), ".", call. = FALSE)
  }
  .check_col(data, actor, "actor")
  .check_col(data, action, "action")
  if (!is.null(time))    .check_col(data, time, "time")
  if (!is.null(order))   .check_col(data, order, "order")
  if (!is.null(session)) .check_col(data, session, "session")

  ev <- as.character(data[[action]])

  tvec <- if (!is.null(time)) {
    .parse_time(data[[time]], custom_format, is_unix_time, unix_time_unit)
  } else {
    NULL
  }
  # Ordering key: explicit `order` wins, else parsed `time`, else the
  # row order in which events were supplied.
  ordkey <- if (!is.null(order)) {
    data[[order]]
  } else if (!is.null(tvec)) {
    tvec
  } else {
    seq_len(nrow(data))
  }

  # Grouping key: actor, optionally crossed with an explicit session id.
  grp <- if (!is.null(session)) {
    interaction(as.character(data[[actor]]),
                as.character(data[[session]]),
                drop = TRUE, lex.order = TRUE)
  } else {
    factor(as.character(data[[actor]]))
  }

  split_on_gap <- !is.null(tvec) && is.null(session) &&
    is.finite(time_threshold)

  idx_by_grp <- split(seq_len(nrow(data)), grp)
  seqs <- lapply(idx_by_grp, function(ix) {
    o <- ix[order(ordkey[ix])]
    evs <- ev[o]
    if (split_on_gap) {
      tt <- as.numeric(tvec[o])
      cut <- cumsum(c(0, diff(tt)) > time_threshold)
      unname(split(evs, cut))
    } else {
      list(evs)
    }
  })
  # Flatten the one-actor-many-sessions nesting into a flat sequence
  # list and drop any empty sequences.
  out <- unlist(seqs, recursive = FALSE, use.names = FALSE)
  out[vapply(out, length, integer(1L)) > 0L]
}

# Coerce a time column to numeric seconds so gaps are comparable to
# `time_threshold`. Handles Unix epochs, custom strptime formats,
# native date/time classes, and already-numeric columns.
.parse_time <- function(v, custom_format, is_unix_time, unix_time_unit) {
  if (isTRUE(is_unix_time)) {
    div <- switch(match.arg(unix_time_unit,
                            c("seconds", "milliseconds", "microseconds")),
                  seconds = 1, milliseconds = 1e3, microseconds = 1e6)
    return(as.numeric(v) / div)
  }
  if (!is.null(custom_format)) {
    return(as.numeric(as.POSIXct(as.character(v),
                                 format = custom_format, tz = "UTC")))
  }
  if (inherits(v, c("POSIXct", "POSIXt", "Date"))) return(as.numeric(v))
  if (is.numeric(v)) return(as.numeric(v))
  parsed <- suppressWarnings(tryCatch(
    as.numeric(as.POSIXct(as.character(v), tz = "UTC")),
    error = function(e) rep(NA_real_, length(v))
  ))
  if (anyNA(parsed)) {
    stop("Could not parse the `time` column as date/time. Supply ",
         "`custom_format` (a strptime format) or `is_unix_time = TRUE`.",
         call. = FALSE)
  }
  parsed
}

# Validate that a named column exists, with an actionable message.
.check_col <- function(data, col, role) {
  if (!is.character(col) || length(col) != 1L) {
    stop(sprintf("`%s` must be a single column name (character).", role),
         call. = FALSE)
  }
  if (!col %in% names(data)) {
    stop(sprintf("`%s` column %s not found in the data.",
                 role, shQuote(col)), call. = FALSE)
  }
  invisible(TRUE)
}
