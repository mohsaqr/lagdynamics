# Coded AI behaviours in human-AI programming sessions

`ai_long` is a coded dataset, in long format, of the behaviours of an AI
coding assistant in "vibe-coding" sessions, in which a person builds
software by instructing the assistant in natural language. Each row is
one coded AI behaviour, with the project and session it belongs to and
its position in the session. A single assistant message can carry
several codes, which are ordered within the message.

The log contains 8,551 coded events in 428 sessions from 34 projects.
Eight behaviour codes are grouped into three clusters: `Action`
(`Delegate`, `Execute`, `Investigate`), `Communication` (`Ask`,
`Explain`, `Report`) and `Repair` (`Plan`, `Repair`). With
`session = "session_id"`, `action = "code"` and
`order = "order_in_session"`,
[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) reads
each session as one sequence and tests which behaviour tends to follow
which.

## Usage

``` r
ai_long
```

## Format

A data frame with 8,551 rows (coded events) and 9 columns:

- `message_id`:

  Integer. Identifier of the assistant message.

- `project`:

  Character. Project identifier (34 projects).

- `session_id`:

  Character. Session identifier (428 sessions, each within one project).

- `timestamp`:

  Integer. Unix time in seconds.

- `session_date`:

  Character. Session date as `YYYY-MM-DD`.

- `code`:

  Character. AI behaviour code (8 values).

- `cluster`:

  Character. Behaviour cluster (3 values).

- `code_order`:

  Integer. Position of the code within its message.

- `order_in_session`:

  Integer. Position of the event within its session.

## Source

Taken without modification from
[`Nestimate::ai_long`](https://pak.dynasite.org/Nestimate/reference/long-data.html)
(Nestimate 0.7.7, MIT license).

## See also

[`lsa()`](https://pak.dynasite.org/lagdynamics/reference/lsa.md) and its
long-format arguments;
[group_regulation_long](https://pak.dynasite.org/lagdynamics/reference/group_regulation_long.md)
for another event log.

## Examples

``` r
fit <- lsa(ai_long, session = "session_id", action = "code",
           order = "order_in_session")
transitions(fit, direction = "over", sort = "strength")
#>           from          to lag count    expected       prob   prob_col
#> 1  Investigate        Plan   1   977  441.607165 0.43713647 0.60872274
#> 2     Delegate        Plan   1   174   55.324388 0.62142857 0.10841121
#> 3          Ask     Explain   1    44    5.859042 0.48351648 0.08413002
#> 4      Execute     Execute   1  1575 1229.837498 0.50970874 0.48716363
#> 5      Explain      Report   1    56   11.230334 0.11111111 0.30939227
#> 6       Repair     Explain   1    43   16.160655 0.17131474 0.08221797
#> 7      Explain      Repair   1    42   15.945833 0.08333333 0.16342412
#> 8         Plan     Execute   1   709  600.988551 0.46953642 0.21930096
#> 9      Explain     Explain   1    66   32.450080 0.13095238 0.12619503
#> 10      Repair         Ask   1    13    2.997292 0.05179283 0.13402062
#> 11        Plan    Delegate   1    73   40.152653 0.04834437 0.33796296
#> 12        Plan Investigate   1   458  373.828635 0.30331126 0.22774739
#> 13        Plan      Report   1    56   33.646436 0.03708609 0.30939227
#> 14      Report      Repair   1    14    5.125446 0.08641975 0.05447471
#> 15     Explain Investigate   1   158  124.774591 0.31349206 0.07856788
#> 16      Report    Delegate   1    10    4.307768 0.06172840 0.04629630
#> 17      Repair     Execute   1   119   99.899421 0.47410359 0.03680792
#> 18      Report     Explain   1    18   10.430383 0.11111111 0.03441683
#> 19        Plan      Repair   1    60   47.774221 0.03973510 0.23346304
#>      adj_res             p   yules_q      kappa   kappa_z       kappa_p
#> 1  33.406415 1.106324e-244 0.7335071 0.45556496 33.405141 1.154465e-244
#> 2  18.126808  1.957990e-73 0.7606295 0.49396057 17.858375  2.487761e-71
#> 3  16.382346  2.556890e-60 0.8731108 0.40817795 15.987533  1.560920e-57
#> 4  16.115668  1.980139e-58 0.3581931 0.16545829 15.299301  7.728241e-53
#> 5  13.950509  3.123680e-44 0.7645397 0.26430961 14.066419  6.108484e-45
#> 6   7.011478  2.358143e-12 0.5219490 0.11295778  7.196059  6.197777e-13
#> 7   6.846113  7.588334e-12 0.5158153 0.10881272  6.931995  4.149462e-12
#> 8   6.293605  3.101773e-10 0.1782988 0.09151372  5.214398  1.844147e-07
#> 9   6.287012  3.236349e-10 0.4050275 0.06889677  6.370623  1.882619e-10
#> 10  5.904423  3.538838e-09 0.6701970 0.10439578  5.934981  2.938803e-09
#> 11  5.823103  5.776481e-09 0.3936539 0.07156401  2.587285  9.673542e-03
#> 12  5.562290  2.662577e-08 0.1730405 0.01612206  1.182180  2.371343e-01
#> 13  4.319460  1.564112e-05 0.3331230 0.14797479  4.161863  3.156616e-05
#> 14  4.023800  5.726663e-05 0.5005555 0.04875854  3.766394  1.656224e-04
#> 15  3.540650  3.991420e-04 0.1738796 0.04192412  1.624637  1.042399e-01
#> 16  2.807910  4.986414e-03 0.4247396 0.02149115  1.545942  1.221186e-01
#> 17  2.501984  1.234995e-02 0.1587659 0.13251699  2.749252  5.973140e-03
#> 18  2.447652  1.437905e-02 0.2971468 0.04065973  2.163890  3.047281e-02
#> 19  1.992140  4.635566e-02 0.1480891 0.05429843  1.828080  6.753747e-02
#>        lift sign significant
#> 1  2.212374 over        TRUE
#> 2  3.145087 over        TRUE
#> 3  7.509760 over        TRUE
#> 4  1.280657 over        TRUE
#> 5  4.986495 over        TRUE
#> 6  2.660783 over        TRUE
#> 7  2.633917 over        TRUE
#> 8  1.179723 over        TRUE
#> 9  2.033893 over        TRUE
#> 10 4.337249 over        TRUE
#> 11 1.818062 over        TRUE
#> 12 1.225160 over        TRUE
#> 13 1.664366 over        TRUE
#> 14 2.731469 over        TRUE
#> 15 1.266283 over        TRUE
#> 16 2.321388 over        TRUE
#> 17 1.191198 over        TRUE
#> 18 1.725728 over        TRUE
#> 19 1.255907 over        TRUE
```
