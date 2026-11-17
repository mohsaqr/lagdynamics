# Implementation notes removed from `?plot_chords` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/plot-chords.R); it is kept here for maintainers.

- **Self-loops**: `cograph::plot_chord()` keeps self-loops regardless of
  its `self_loop` flag, so `self_loops = FALSE` zeros the diagonal of the
  width matrix before drawing.
- **Colour order**: cograph draws chords row-major over non-zero widths;
  the colour vector is built in that order (`.lsa_chord_colors()`),
  opaque, with transparency from `chord_alpha`. Diverging scales are
  symmetric about 0 (`#2166AC` / `#F7F7F7` / `#B2182B`); sequential
  scales span the observed range (`#DEEBF7` -> `#08519C`); non-finite
  values are `grey85`.
- **compare**: states matched by name; `significant` is ignored.
- **API gap**: comparing two groups needs two separate fits built with
  `subset()`; a grouped fit (`lsa_group`) cannot be passed directly.
