# Implementation notes removed from `?engagement`, `?group_regulation`, `?group_regulation_long` and `?ai_long` (2026-11-17)

Material cut from the help page when it was rewritten for readers. It is
still true of the data (R/data-*.R); it is kept here for maintainers.

- The previous page described the time points as weekly; the source
  (`Nestimate` `trajectories`, synthetic) does not establish the unit, so
  the page now says "time points". Mean observed length is 13.75 (136
  non-empty rows).
- Data provenance claim "`ai_long` derived from `Nestimate::ai_long`
  (Nestimate 0.7.7)" and "`group_regulation_long` from
  `tna::group_regulation_long`" are kept from the old pages. Note that
  the sequences of `group_regulation_long` (ordered by Actor and Time) are
  not identical row by row to `group_regulation`, so the pages do not
  claim they are the same events.
