# Engine registry. The set of LSA engines is open: users can add their
# own via `register_lsa_engine()`. Built-in engines register themselves
# in `R/zzz.R` at load time. The registry lives in a package-private
# environment so it survives the lifetime of the session.

.lsa_engine_registry <- new.env(parent = emptyenv())

#' Register a lag sequential analysis engine
#'
#' @description
#' `register_lsa_engine()` adds a lag sequential analysis engine under a name,
#' after which `lsa(data, engine = name)` fits it; the fit is read, resampled
#' and plotted like any other. An engine is the part of a lag sequential
#' analysis that turns a table of observed transitions into expected counts and
#' cell-wise tests. The built-in engines implement the classical analysis of
#' Bakeman and Quera (1995), a two-cell test, Sackett's (1979) bidirectional
#' test, and the parallel and non-parallel dominance tests (Wampold, 1984).
#' Other tests of the same table, such as a different residual or a model-based
#' expectation, can be added as engines and then used with every function of the
#' package.
#'
#' Re-registering a name replaces the engine for the session. The built-in
#' engines are registered when the package loads.
#'
#' @param name Name of the engine, used as `lsa(engine = name)`.
#' @param fn The engine function. It is called with `transitions` (the
#'   [lsa_transitions()] object), `structural_zeros` (`NULL` or a `K x K`
#'   0/1 matrix), `alternative`, `n_events` (`NULL` for transition-matrix
#'   input) and any engine parameters given to [lsa()], so it should accept
#'   `transitions` and `...`. It must return a named list holding at least
#'   the `K x K` matrices `obs`, `exp`, `prob`, `adj_res` and `p`; further
#'   `K x K` matrices become extra columns of [transitions()].
#' @param description One-line description shown by [list_lsa_engines()].
#' @param requires Names of packages the engine needs (default none);
#'   [lsa()] checks that they are installed before fitting.
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and `lagdynamics_bad_input`
#' when `fn` has no `transitions` argument.
#'
#' @return `name`, invisibly.
#'
#' @examples
#' # An engine that tests Pearson residuals under independence
#' pearson_engine <- function(transitions, ...) {
#'   obs <- transitions$obs
#'   expected <- outer(rowSums(obs), colSums(obs)) / sum(obs)
#'   residual <- (obs - expected) / sqrt(expected)
#'   list(obs = obs, exp = expected, prob = obs / rowSums(obs),
#'        adj_res = residual, p = 2 * stats::pnorm(-abs(residual)))
#' }
#' register_lsa_engine("pearson", pearson_engine,
#'                     description = "Pearson residuals under independence")
#' list_lsa_engines()
#' fit <- lsa(engagement, engine = "pearson")
#' transitions(fit, significant = TRUE)
#' unregister_lsa_engine("pearson")
#'
#' @references
#' Bakeman, R., & Quera, V. (1995). \emph{Analyzing interaction:
#' Sequential analysis with SDIS and GSEQ}. Cambridge University Press.
#'
#' Sackett, G. P. (1979). The lag sequential analysis of contingency and
#' cyclicity in behavioral interaction research. In J. D. Osofsky (Ed.),
#' \emph{Handbook of infant development} (pp. 623-649). Wiley.
#'
#' Wampold, B. E. (1984). Tests of dominance in sequential categorical
#' data. \emph{Psychological Bulletin, 96}(2), 424-429.
#'
#' @seealso [list_lsa_engines()], [get_lsa_engine()],
#'   [unregister_lsa_engine()], [lsa()]
#'
#' @export
register_lsa_engine <- function(name, fn, description, requires = character()) {
  stopifnot(
    "`name` must be a single non-empty string" =
      is.character(name) && length(name) == 1L && nzchar(name),
    "`fn` must be a function" = is.function(fn),
    "`description` must be a single string" =
      is.character(description) && length(description) == 1L,
    "`requires` must be a character vector" = is.character(requires)
  )
  if (!"transitions" %in% names(formals(fn))) {
    .lsa_stop("Engine function must accept a `transitions` argument.")
  }
  entry <- list(name = name, fn = fn, description = description,
                requires = requires)
  assign(name, entry, envir = .lsa_engine_registry)
  invisible(name)
}

#' Retrieve a registered engine
#'
#' @description
#' `get_lsa_engine()` returns the registry entry of one engine: its name, its
#' function, its description and the packages it requires. It is the lookup
#' [lsa()] performs before fitting, and it is useful for inspecting a built-in
#' engine's function or for checking that a registered engine can run in the
#' current session. To see which engines exist, use [list_lsa_engines()].
#'
#' @param name Name of the engine.
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and one of:
#' * `lagdynamics_unknown_engine`: no engine is registered under `name`.
#' * `lagdynamics_missing_package`: a package the engine requires is not
#'   installed.
#'
#' @return A list with the elements `name`, `fn`, `description` and
#'   `requires`, as given to [register_lsa_engine()].
#'
#' @examples
#' str(get_lsa_engine("classical"), give.attr = FALSE)
#'
#' @seealso [list_lsa_engines()], [register_lsa_engine()]
#'
#' @export
get_lsa_engine <- function(name) {
  stopifnot(
    "`name` must be a single string" = is.character(name) && length(name) == 1L
  )
  if (!exists(name, envir = .lsa_engine_registry, inherits = FALSE)) {
    available <- ls(envir = .lsa_engine_registry)
    .lsa_stop(class = "unknown_engine", sprintf(
      "Engine '%s' is not registered. Available engines: %s",
      name, paste(sprintf("'%s'", available), collapse = ", ")
    ))
  }
  entry <- get(name, envir = .lsa_engine_registry, inherits = FALSE)
  installed <- vapply(entry$requires, requireNamespace, logical(1L),
                      quietly = TRUE)
  if (!all(installed)) {
    .lsa_stop(class = "missing_package", sprintf(
      "Engine '%s' requires the package(s) %s, which are not installed.",
      name, paste(sprintf("'%s'", entry$requires[!installed]), collapse = ", ")
    ))
  }
  entry
}

#' List the registered engines
#'
#' @description
#' `list_lsa_engines()` lists the engines available in the session, the five
#' built-in ones and any added with [register_lsa_engine()], with a one-line
#' description of each. Each engine of [lsa()] implements one test of the
#' transition table. Use the list to choose the `engine` argument of [lsa()].
#'
#' @return A `data.frame` with one row per registered engine, in
#'   alphabetical order, and the columns `name`, `description` and
#'   `requires` (comma-separated package names, `""` when none).
#'
#' @examples
#' list_lsa_engines()
#'
#' @seealso [lsa()], [register_lsa_engine()], [get_lsa_engine()]
#'
#' @export
list_lsa_engines <- function() {
  nms <- ls(envir = .lsa_engine_registry)
  if (length(nms) == 0L) {
    return(data.frame(name = character(0),
                      description = character(0),
                      requires = character(0),
                      stringsAsFactors = FALSE))
  }
  entries <- lapply(nms, function(nm) {
    get(nm, envir = .lsa_engine_registry, inherits = FALSE)
  })
  data.frame(
    name = vapply(entries, `[[`, character(1), "name"),
    description = vapply(entries, `[[`, character(1), "description"),
    requires = vapply(entries,
                      function(e) paste(e$requires, collapse = ", "),
                      character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Remove a registered engine
#'
#' @description
#' `unregister_lsa_engine()` removes an engine from the registry for the rest of
#' the session, so that `lsa(engine = name)` no longer finds it. It is mainly
#' used to clean up an engine added with [register_lsa_engine()]. A removed
#' built-in engine returns when the package is loaded again.
#'
#' @param name Name of the engine.
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and
#' `lagdynamics_unknown_engine` when no engine is registered under `name`.
#'
#' @return `NULL`, invisibly.
#'
#' @examples
#' placeholder <- function(transitions, ...) {
#'   obs <- transitions$obs
#'   expected <- outer(rowSums(obs), colSums(obs)) / sum(obs)
#'   list(obs = obs, exp = expected, prob = obs / rowSums(obs),
#'        adj_res = obs * 0, p = obs * 0 + 1)
#' }
#' register_lsa_engine("placeholder", placeholder,
#'                     description = "Placeholder engine")
#' unregister_lsa_engine("placeholder")
#' list_lsa_engines()
#'
#' @seealso [register_lsa_engine()], [list_lsa_engines()]
#'
#' @export
unregister_lsa_engine <- function(name) {
  stopifnot(
    "`name` must be a single string" = is.character(name) && length(name) == 1L
  )
  if (!exists(name, envir = .lsa_engine_registry, inherits = FALSE)) {
    .lsa_stop(class = "unknown_engine",
              sprintf("Engine '%s' is not registered.", name))
  }
  rm(list = name, envir = .lsa_engine_registry)
  invisible(NULL)
}
