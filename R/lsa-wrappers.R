# Convenience wrappers. Each one pins a specific engine.

#' @rdname lsa
#' @export
lsa_classical <- function(data, ...) {
  lsa(data, engine = "classical", ...)
}

#' @rdname lsa
#' @export
lsa_two_cell <- function(data, ...) {
  lsa(data, engine = "two_cell", ...)
}

#' @rdname lsa
#' @export
lsa_bidirectional <- function(data, ...) {
  lsa(data, engine = "bidirectional", ...)
}

#' @rdname lsa
#' @export
lsa_parallel_dominance <- function(data, ...) {
  lsa(data, engine = "parallel_dominance", ...)
}

#' @rdname lsa
#' @export
lsa_nonparallel_dominance <- function(data, ...) {
  lsa(data, engine = "nonparallel_dominance", ...)
}
