# Implementation notes removed from `?plot_polar` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the code (R/plot-polar.R); it is kept here for maintainers.

- **Geometry** (`.lsa_polar_geom`): source band radii 0.34-0.50, target
  ring 0.54-1; labels set tangentially and flipped on the lower half
  (`.lsa_tangential()`).
- **Rose**: equal sectors per source and equal slots per target; bar
  height is the count or probability, square-root scaled.
- **Wedge**: sectors sized by outgoing volume; wedges below `min_show`
  of their source's outflow are omitted; `labels = "auto"` names only
  wide wedges.
- **Errors**: "Fit has no transitions to draw" and "Fit has no outgoing
  transitions to draw" (`lagdynamics_bad_input`).
