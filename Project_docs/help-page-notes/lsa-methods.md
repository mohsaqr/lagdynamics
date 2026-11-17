# Implementation notes removed from `?lsa-methods` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/lsa.R); it is kept here for maintainers.

- **ASCII fallback**: printed output uses a block glyph for bars and a
  superscript in G-squared when the session is UTF-8, and plain ASCII
  (`#`, `G^2`) otherwise (`.lsa_utf8()`, `.lsa_glyph()` in
  R/conditions.R).
- **print.lsa()** lists at most five strongest over-represented
  transitions (by `adj_res`, significant only) and the initial-state
  bar chart sorted by probability; it shows at most 12 state names.
- **print.lsa_group()** shows at most 10 state names.
