# Long-format event-log sequencing inside lsa() (actor / action / time
# / order / session / time_threshold). The log is sequenced into event
# sequences before analysis; these tests pin the grouping, ordering,
# session-splitting, time parsing, and validation behaviour.

make_log <- function() {
  data.frame(
    user = c("u1", "u1", "u1", "u1", "u2", "u2", "u2"),
    act  = c("login", "read", "quiz", "logout",
             "login", "read", "quiz"),
    ts   = as.POSIXct(c("2025-01-01 10:00:00", "2025-01-01 10:01:00",
                        "2025-01-01 10:02:00",
                        "2025-01-01 11:00:00",   # 58-min gap -> split
                        "2025-01-02 09:00:00", "2025-01-02 09:01:00",
                        "2025-01-02 09:02:00"), tz = "UTC"),
    stringsAsFactors = FALSE
  )
}

test_that("time-ordered log splits sessions on a gap > time_threshold", {
  fit <- lsa(make_log(), actor = "user", action = "act", time = "ts",
             time_threshold = 900)
  # u1's 58-min gap splits it into 2 sessions; u2 stays 1 -> 3 total.
  expect_equal(fit$data$n_sequences, 3L)
  # Session lengths: u1 {login,read,quiz}=3, u1 {logout}=1, u2=3.
  expect_setequal(fit$data$transitions_per_seq + 1L, c(3L, 1L, 3L))
})

test_that("a large time_threshold keeps each actor as one session", {
  fit <- lsa(make_log(), actor = "user", action = "act", time = "ts",
             time_threshold = 1e9)
  expect_equal(fit$data$n_sequences, 2L)
})

test_that("no time column yields one session per actor", {
  fit <- lsa(make_log(), actor = "user", action = "act")
  expect_equal(fit$data$n_sequences, 2L)
})

test_that("an explicit session column overrides gap splitting", {
  log <- make_log()
  log$sess <- c(1, 1, 1, 1, 1, 1, 1)   # force u1 into a single session
  fit <- lsa(log, actor = "user", action = "act", time = "ts",
             session = "sess")
  # No gap splitting when session is given: 2 actors x 1 session = 2.
  expect_equal(fit$data$n_sequences, 2L)
})

test_that("unix-time and custom-format columns parse equivalently", {
  log <- make_log()
  log$ux <- as.numeric(log$ts)
  fit_unix <- lsa(log, actor = "user", action = "act", time = "ux",
                  is_unix_time = TRUE, time_threshold = 900)
  expect_equal(fit_unix$data$n_sequences, 3L)

  log$str <- format(log$ts, "%Y-%m-%d %H:%M:%S")
  fit_fmt <- lsa(log, actor = "user", action = "act", time = "str",
                 custom_format = "%Y-%m-%d %H:%M:%S", time_threshold = 900)
  expect_equal(fit_fmt$data$n_sequences, 3L)
})

test_that("order column overrides supplied row order", {
  log <- make_log()[c(3, 1, 2, 4, 7, 5, 6), ]   # scramble rows
  log$ord <- c(3, 1, 2, 4, 3, 1, 2)             # true within-actor order
  fit <- lsa(log, actor = "user", action = "act", order = "ord")
  # u1 reconstructed as login->read->quiz->logout (one session, no time).
  e <- fit$edges
  expect_true(any(e$from_label == "login" & e$to_label == "read" &
                    e$count >= 1))
})

test_that("long-format validation errors are actionable", {
  log <- make_log()
  expect_error(lsa(log, action = "act"),
               "needs both `actor` and `action`")
  expect_error(lsa(log, actor = "user", action = "missing"),
               "column 'missing' not found")
  expect_error(lsa(log, actor = "user", action = "act",
                   group = c("a", "b")),
               "cannot be combined with long-format")
})

test_that("unparseable time without a format errors clearly", {
  log <- make_log()
  log$bad <- c("foo", "bar", "baz", "qux", "a", "b", "c")
  expect_error(lsa(log, actor = "user", action = "act", time = "bad"),
               "Could not parse the `time` column")
})
