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
  parsed <- suppressWarnings(
    if (isTRUE(is_unix_time)) {
      div <- switch(match.arg(unix_time_unit,
                              c("seconds", "milliseconds", "microseconds")),
                    seconds = 1, milliseconds = 1e3, microseconds = 1e6)
      as.numeric(v) / div
    } else if (!is.null(custom_format)) {
      as.numeric(as.POSIXct(as.character(v),
                            format = custom_format, tz = "UTC"))
    } else if (inherits(v, c("POSIXct", "POSIXt", "Date"))) {
      as.numeric(v)
    } else if (is.numeric(v)) {
      as.numeric(v)
    } else {
      tryCatch(as.numeric(as.POSIXct(as.character(v), tz = "UTC")),
               error = function(e) rep(NA_real_, length(v)))
    }
  )
  # Every branch is NA-checked, not just the character fallback: an
  # unparseable custom format, a missing POSIXct/numeric value, or a bad
  # Unix epoch would otherwise yield NA times that silently drop their
  # events during ordering/session splitting and corrupt the counts.
  bad <- which(!is.finite(parsed))
  if (length(bad)) {
    stop(sprintf(
      paste0("`time` could not be parsed for %d row(s) (e.g. row %s). ",
             "Fix the values, supply `custom_format`, or set ",
             "`is_unix_time = TRUE`; lagseq will not silently drop ",
             "events with unparseable times."),
      length(bad), paste(utils::head(bad, 5L), collapse = ", ")),
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
