# Classed conditions and console glyphs.
#
# Every user-facing error is raised through .lsa_stop(), which builds an
# errorCondition carrying two classes: a specific one
# ("lagdynamics_<class>") and the umbrella "lagdynamics_error". Callers and
# tests catch the specific class; the message text is free to change.
#
# Classes in use (documented in ?lsa, section "Conditions"):
#   lagdynamics_bad_input           an argument or input value is invalid
#   lagdynamics_needs_events        the fit came from a transition matrix,
#                                   but the verb resamples sequences
#   lagdynamics_no_transitions      the data contain zero transitions
#   lagdynamics_unknown_engine      the engine name is not registered
#   lagdynamics_unsupported_option  the engine cannot honour an option
#                                   (e.g. structural zeros)
#   lagdynamics_incompatible_fits   fits to be compared differ in states,
#                                   engine, lag or structural zeros
#   lagdynamics_missing_package     a suggested package is not installed
#   lagdynamics_ipf_infeasible      iterative proportional fitting has no
#                                   solution for the structural zeros
#
# Classed warnings (warningCondition):
#   lagdynamics_structural_zero_counts  observed transitions fall in cells
#                                   declared structural zeros (lsa())
#   lagdynamics_tna_no_data         lsa_to_tna() with weights that tna
#                                   cannot resample omits sequence data

# Signal a classed lagdynamics error. `...` is pasted as stop() would.
.lsa_stop <- function(..., class = "bad_input", call. = FALSE) {
  stop(errorCondition(
    .makeMessage(...),
    class = c(paste0("lagdynamics_", class), "lagdynamics_error"),
    call = NULL
  ))
}

# TRUE when the session can print UTF-8 (the console bar glyph, G-squared,
# the middle dot). Under LC_ALL=C or another non-UTF-8 locale the ASCII
# fallback is used so output never degrades to "<U+2588>".
.lsa_utf8 <- function() isTRUE(l10n_info()[["UTF-8"]])

# Pick the UTF-8 glyph or its ASCII fallback.
.lsa_glyph <- function(utf8, ascii) if (.lsa_utf8()) utf8 else ascii

# Scoped seeding. A verb that takes `seed` sets it for its own draws and
# restores the caller's RNG state on exit, so calling it never changes the
# random stream of the surrounding script.
.lsa_get_seed <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

.lsa_restore_seed <- function(old_seed) {
  if (is.null(old_seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", old_seed, envir = globalenv())
  }
  invisible(NULL)
}

# Replace the middle dot used in plot subtitles by "|" outside UTF-8, so a
# non-UTF-8 graphics device never meets a character it cannot encode.
.lsa_ascii <- function(x) {
  if (.lsa_utf8()) x else gsub("\u00b7", "|", x, fixed = TRUE)
}
