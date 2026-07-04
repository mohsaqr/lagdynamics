# Engine registry. The set of LSA engines is open: users can add their
# own via `register_lsa_engine()`. Built-in engines register themselves
# in `R/zzz.R` at load time. The registry lives in a package-private
# environment so it survives the lifetime of the session.

.lsa_engine_registry <- new.env(parent = emptyenv())

#' Register a Lag Sequential Analysis Engine
#'
#' Adds a new engine to the lagdynamics registry so it can be referenced by
#' name via `lsa(..., engine = "<name>")`. Built-in engines
#' (`"classical"`, `"two_cell"`, `"bidirectional"`,
#' `"parallel_dominance"`, `"nonparallel_dominance"`) are registered
#' automatically when the package loads.
#'
#' @param name Character scalar. The engine's identifier as used in
#'   `lsa(engine = name)`.
#' @param fn A function. Must accept a `transitions` argument (a tidy
#'   transition table produced by [lsa_transitions()]) and arbitrary
#'   named `...` arguments forwarded from `lsa(params = list(...))`.
#'   Must return a named list with at least the matrix elements `obs`,
#'   `exp`, `prob`, `adj_res`, and `p` (each `K x K`). Additional `K x K`
#'   matrices are preserved as engine-specific edge statistics.
#' @param description Character scalar. One-line human-readable
#'   description shown by [list_lsa_engines()].
#' @param requires Character vector. Names of packages the engine
#'   depends on. Empty by default.
#'
#' @return Invisibly returns `name`.
#'
#' @examples
#' my_engine <- function(transitions, ...) {
#'   get_lsa_engine("classical")$fn(transitions, ...)
#' }
#' register_lsa_engine("my_classical", my_engine, "Classical test alias")
#' fit <- lsa(engagement, engine = "my_classical")
#' unregister_lsa_engine("my_classical")
#' fit
#'
#' @seealso [get_lsa_engine()], [list_lsa_engines()],
#'   [unregister_lsa_engine()], [lsa()]
#'
#' @export
register_lsa_engine <- function(name, fn, description, requires = character()) {
  stopifnot(
    is.character(name), length(name) == 1L, nzchar(name),
    is.function(fn),
    is.character(description), length(description) == 1L,
    is.character(requires)
  )
  if (!"transitions" %in% names(formals(fn))) {
    stop("Engine function must accept a `transitions` argument.",
         call. = FALSE)
  }
  entry <- list(name = name, fn = fn, description = description,
                requires = requires)
  assign(name, entry, envir = .lsa_engine_registry)
  invisible(name)
}

#' Retrieve a Registered LSA Engine
#'
#' @param name Character scalar. The engine's identifier.
#'
#' @return The registry entry: a list with elements `name`, `fn`,
#'   `description`, `requires`.
#'
#' @examples
#' classical <- get_lsa_engine("classical")
#' classical$description
#'
#' @seealso [register_lsa_engine()], [list_lsa_engines()]
#'
#' @export
get_lsa_engine <- function(name) {
  stopifnot(is.character(name), length(name) == 1L)
  if (!exists(name, envir = .lsa_engine_registry, inherits = FALSE)) {
    available <- ls(envir = .lsa_engine_registry)
    stop(sprintf(
      "Engine '%s' is not registered. Available engines: %s",
      name, paste(sprintf("'%s'", available), collapse = ", ")
    ), call. = FALSE)
  }
  entry <- get(name, envir = .lsa_engine_registry, inherits = FALSE)
  for (pkg in entry$requires) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf(
        "Engine '%s' requires the '%s' package, which is not installed.",
        name, pkg
      ), call. = FALSE)
    }
  }
  entry
}

#' List All Registered LSA Engines
#'
#' @return A data.frame with columns `name`, `description`, `requires`.
#'
#' @examples
#' list_lsa_engines()
#'
#' @seealso [register_lsa_engine()], [get_lsa_engine()]
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

#' Remove a Registered LSA Engine
#'
#' @param name Character scalar. The engine's identifier.
#'
#' @return Invisibly `NULL`.
#'
#' @examples
#' my_engine <- function(transitions, ...) {
#'   get_lsa_engine("classical")$fn(transitions, ...)
#' }
#' register_lsa_engine("temporary_engine", my_engine, "Temporary alias")
#' unregister_lsa_engine("temporary_engine")
#'
#' @seealso [register_lsa_engine()]
#'
#' @export
unregister_lsa_engine <- function(name) {
  stopifnot(is.character(name), length(name) == 1L)
  if (!exists(name, envir = .lsa_engine_registry, inherits = FALSE)) {
    stop(sprintf("Engine '%s' is not registered.", name), call. = FALSE)
  }
  rm(list = name, envir = .lsa_engine_registry)
  invisible(NULL)
}
