# Implementation notes removed from `?register_lsa_engine` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/engine-registry.R); it is kept here for maintainers.

- **Registry**: a package-private environment `.lsa_engine_registry`;
  built-in engines are registered in `.onLoad()` (R/zzz.R).
- **Argument checks** other than the `transitions` formal are plain
  `stopifnot()` (unclassed).
- **`requires`** is checked by `get_lsa_engine()` with
  `requireNamespace()`, i.e. at fit time, not at registration.
