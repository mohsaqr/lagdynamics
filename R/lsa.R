# Main constructor. Reads canonical lsa_data + lsa_transitions, dispatches
# to the named engine via the registry, and assembles the S3 fit object
# with stable named slots.

#' Fit a lag sequential analysis model
#'
#' @description
#' `lsa()` fits a lag sequential analysis model to event sequences, wide
#' sequence tables, transition-count matrices, long event logs, or sequence
#' objects from other packages. Lag sequential analysis (Sackett, 1979; Bakeman
#' & Gottman, 1997) tests whether, in a stream of coded events, one state
#' follows another more or less often than chance predicts.
#'
#' The transitions observed at a given lag are counted into a table of
#' antecedent (`from`) by consequent (`to`) states; an expected table is derived
#' under the hypothesis that the two are independent; and each cell is tested
#' with an adjusted residual (Haberman, 1979), a standard normal deviate that is
#' positive when the transition is over-represented and negative when it is
#' avoided. The significant cells define a lag transition network, and the
#' tablewise likelihood-ratio and Pearson tests show whether there is any
#' sequential dependence at all. For every cell, `lsa()` estimates the observed
#' and expected counts, the transition probability, the adjusted residual and
#' its p-value, Yule's Q, kappa and lift. The `engine` selects the test: the
#' classical analysis of Bakeman and Quera (1995) by default; a two-cell (2 x 2)
#' test with odds ratios; Sackett's bidirectional test on the symmetrised table;
#' or the parallel and non-parallel dominance tests of Sackett (1979) and
#' Wampold (1984). Self-transitions or any other cells that cannot occur by
#' construction are declared as structural zeros, and the expected table is then
#' fitted by iterative proportional fitting (Christensen, 1997). With `group`,
#' one model is fitted per group on a shared state set so the groups can be
#' compared.
#'
#' It returns an `lsa` object, or with `group` an `lsa_group` object holding one
#' fit per group.
#'
#' The fit is read with [transitions()], [nodes()], [tests()], [initial()]
#' and [transition_probabilities()], and drawn with [plot.lsa()]. To fit
#' several lags at once use [lsa_lags()]; to assess the transitions use
#' [bootstrap_lsa()], [certainty_lsa()] or [permute_lsa()]; to compare
#' groups use [compare_lsa()]. The engine wrappers `lsa_classical()`,
#' `lsa_two_cell()`, `lsa_bidirectional()`, `lsa_parallel_dominance()` and
#' `lsa_nonparallel_dominance()` are `lsa()` with `engine` fixed and return
#' the same object.
#'
#' @param data Sequences in any form accepted by [lsa_data()]: a vector, a
#'   list of sequences, a wide table (one row per sequence), a square
#'   transition-count matrix, or a `tna`, `group_tna`, `tna_seq_data`,
#'   `tna_data`, `nestimate_data`, `netobject` or `stslist` object. With
#'   `action`, a long event log (one row per event). `NA` and empty
#'   strings are missing values, not a state; to model missingness as a
#'   state, recode it first.
#' @param lag Integer lag (default `1`). Positive lags pair each event with
#'   the one `lag` steps later, negative lags with the one `|lag|` steps
#'   earlier. A transition-count matrix supports `lag = 1` only.
#' @param engine Name of a registered engine (default `"classical"`); the
#'   built-in engines are listed by [list_lsa_engines()].
#' @param alternative Alternative hypothesis for the residual and kappa
#'   p-values: `"two.sided"` (default), `"greater"` or `"less"`.
#' @param alpha Significance level that flags a transition as significant
#'   (default `0.05`).
#' @param loops Keep self-transitions in the model (default `TRUE`).
#'   `FALSE` declares every self-transition a structural zero.
#' @param structural_zeros Optional `K x K` 0/1 matrix in which `0` marks a
#'   transition that cannot occur; default `NULL` (none). Combined with
#'   `loops = FALSE`, the diagonal is also forbidden. Supported by the
#'   classical engine.
#' @param labels Optional character vector of state names. Default `NULL`
#'   takes them from the data (sorted values of character input).
#' @param group Optional grouping: one value per sequence, or, with a long
#'   event log, the name of a column that is constant within each
#'   sequence. Default `NULL` (one model).
#' @param actor,action,time,order,session Column names for a long event
#'   log (all default `NULL`). Supplying `action` (the state column)
#'   switches to long-format input. Rows are grouped into sequences by
#'   `actor`, `session`, or both, ordered by `order` or else by `time`.
#' @param time_threshold Gap in seconds that starts a new sequence when
#'   `time` is given and `session` is not (default `900`).
#' @param custom_format Optional `strptime()` format for the `time` column
#'   (default `NULL`).
#' @param is_unix_time Read the `time` column as a Unix epoch (default
#'   `FALSE`).
#' @param unix_time_unit Unit of the epoch: `"seconds"` (default),
#'   `"milliseconds"` or `"microseconds"`.
#' @param params Named list of engine-specific parameters (default
#'   `list()`).
#' @param ... Engine-specific parameters, merged into `params`. For the
#'   engine wrappers, any argument of `lsa()` other than `engine`.
#'
#' @section Long-format input:
#' A long event log is turned into sequences before the model is fitted:
#' one sequence per `actor`, per `session`, or per actor-session pair, or a
#' single sequence when neither is given. Events are ordered within a
#' sequence by `order`, or else by `time`; when `time` is given without
#' `session`, a gap longer than `time_threshold` seconds starts a new
#' sequence. If the session identifier is already unique across actors,
#' pass it as `actor` and leave `session` unset. Check the number of
#' sequences reported by `summary(fit)` against the number expected.
#'
#' @section Conditions:
#' Errors carry the class `lagdynamics_error` and one of:
#' * `lagdynamics_bad_input`: invalid data, labels, group, structural-zero
#'   matrix or long-format column.
#' * `lagdynamics_unknown_engine`: `engine` is not registered.
#' * `lagdynamics_missing_package`: the engine needs a package that is not
#'   installed.
#' * `lagdynamics_no_transitions`: the data contain no transition at the
#'   requested lag.
#' * `lagdynamics_unsupported_option`: the engine does not support
#'   structural zeros.
#' * `lagdynamics_ipf_infeasible`: the structural zeros leave a state with
#'   observed transitions but no permitted cell.
#' * `lagdynamics_needs_events`: `group` was given with a
#'   transition-count matrix, which cannot be split.
#'
#' With structural zeros (`loops = FALSE` or `structural_zeros`), the
#' expected counts, adjusted residuals and the tablewise G^2 come from the
#' quasi-independence model fitted to the possible cells (Bakeman & Quera,
#' 1995, \emph{Psychological Bulletin}): its margins are the totals of the
#' possible cells, so transitions that the data contain in a forbidden cell
#' are left out of the fit. The observed counts and transition probabilities
#' still report them. A warning of class `lagdynamics_structural_zero_counts`
#' is raised when this happens, naming the cells, because the data
#' contradict the declaration.
#'
#' @return An object of class `c("lsa", "cograph_network")`. Read it with
#'   [transitions()] (one row per `from -> to` pair: `from`, `to`, `lag`,
#'   `count`, `expected`, `prob`, `prob_col`, `adj_res`, `p`, `yules_q`,
#'   `kappa`, `kappa_z`, `kappa_p`, `lift`, `sign`, `significant`, plus any
#'   engine-specific columns), [nodes()] (one row per state), [tests()]
#'   (one row per tablewise test), [initial()] (one row per state),
#'   [transition_probabilities()] and `summary()` (see [lsa-methods]).
#'   The object also stores the per-cell quantities as `K x K` matrices,
#'   the canonical [lsa_data()] input and the parameters of the call, for
#'   the inference and plotting functions.
#'
#'   With `group`, an object of class `c("lsa_group", "list")`: one `lsa`
#'   fit per group level. The reading verbs return the same tables with a
#'   leading `group` column, and `summary()` returns one row per group.
#'
#' @examples
#' fit <- lsa(engagement)
#' fit
#' transitions(fit, significant = TRUE, sort = "strength")
#' tests(fit)
#'
#' # The same model through an engine wrapper, at lag 2
#' lsa_classical(engagement, lag = 2)
#'
#' # A long event log, one model per achievement group
#' gfit <- lsa(group_regulation_long, actor = "Actor", action = "Action",
#'             time = "Time", group = "Achiever")
#' summary(gfit)
#'
#' @references
#' Bakeman, R., & Gottman, J. M. (1997). \emph{Observing interaction: An
#' introduction to sequential analysis} (2nd ed.). Cambridge University
#' Press.
#'
#' Bakeman, R., & Quera, V. (1995). \emph{Analyzing interaction:
#' Sequential analysis with SDIS and GSEQ}. Cambridge University Press.
#'
#' Bakeman, R., & Quera, V. (1995). Log-linear approaches to lag-sequential
#' analysis when consecutive codes may and cannot repeat. \emph{Psychological
#' Bulletin, 118}(2), 272-284. \doi{10.1037/0033-2909.118.2.272}
#'
#' Christensen, R. (1997). \emph{Log-linear models and logistic
#' regression} (2nd ed.). Springer.
#'
#' Haberman, S. J. (1979). \emph{Analysis of qualitative data: Volume 2,
#' New developments}. Academic Press.
#'
#' Sackett, G. P. (1979). The lag sequential analysis of contingency and
#' cyclicity in behavioral interaction research. In J. D. Osofsky (Ed.),
#' \emph{Handbook of infant development} (pp. 623-649). Wiley.
#'
#' Wampold, B. E. (1984). Tests of dominance in sequential categorical
#' data. \emph{Psychological Bulletin, 96}(2), 424-429.
#'
#' The formula-by-formula source map is installed with the package at
#' `system.file("REFERENCES.md", package = "lagdynamics")`.
#'
#' @seealso [transitions()], [nodes()], [tests()], [initial()] to read the
#'   fit; [lsa_lags()] for several lags; [bootstrap_lsa()],
#'   [permute_lsa()] and [compare_lsa()] for inference; [plot.lsa()];
#'   [list_lsa_engines()] and [register_lsa_engine()] for engines;
#'   [lsa_data()] for the accepted input forms.
#'
#' @export
lsa <- function(data,
                lag = 1,
                engine = "classical",
                alternative = c("two.sided", "greater", "less"),
                alpha = 0.05,
                loops = TRUE,
                structural_zeros = NULL,
                labels = NULL,
                group = NULL,
                actor = NULL,
                action = NULL,
                time = NULL,
                order = NULL,
                session = NULL,
                time_threshold = 900,
                custom_format = NULL,
                is_unix_time = FALSE,
                unix_time_unit = "seconds",
                params = list(),
                ...) {
  call <- match.call()
  alternative <- match.arg(alternative)
  stopifnot(
    "`alpha` must be a single number in (0, 1)" =
      is.numeric(alpha) && length(alpha) == 1L && alpha > 0 && alpha < 1
  )

  # Long-format mode: a raw event log is sequenced into a list of event
  # sequences before any analysis. Triggered by supplying `action`.
  if (!is.null(action)) {
    # With long-format input, `group` (if given) is the NAME of a grouping
    # column in the log; .prepare_long() derives one label per recovered
    # sequence and returns them as an attribute, which then drives the
    # grouped fit below -- no manual split-then-relabel ritual.
    long_group <- NULL
    if (!is.null(group)) {
      if (!is.character(group) || length(group) != 1L) {
        .lsa_stop("With long-format input, `group` must be the name of a ",
             "grouping column in the log (a single string).",
             call. = FALSE)
      }
      long_group <- group
    }
    data <- .prepare_long(
      data, actor = actor, action = action, time = time, order = order,
      session = session, group = long_group, time_threshold = time_threshold,
      custom_format = custom_format, is_unix_time = is_unix_time,
      unix_time_unit = unix_time_unit
    )
    if (!is.null(long_group)) {
      group <- attr(data, "group")
      attr(data, "group") <- NULL
    }
  }

  if (!is.null(group)) {
    return(.lsa_grouped(
      data = data, group = group, lag = lag, engine = engine,
      alternative = alternative, alpha = alpha, loops = loops,
      structural_zeros = structural_zeros, labels = labels,
      params = params, call = call, ...
    ))
  }

  d  <- lsa_data(data, labels = labels)
  structural_zeros <- .resolve_structural_zeros(
    structural_zeros, loops = loops, K = d$n_states, labels = d$labels)
  tx <- lsa_transitions(d, lag = lag)

  # Attach event-level totals for engines that need them (kappa).
  if (identical(d$source, "events")) {
    n_events <- d$n_events
    event_totals_col <- tabulate(d$events, nbins = d$n_states)
    attr(tx, "event_totals_col") <- event_totals_col
  } else {
    n_events <- NULL
  }

  entry <- get_lsa_engine(engine)
  fn <- entry$fn

  # Merge `...` into params (caller can use either)
  extra <- list(...)
  if (length(extra)) params <- utils::modifyList(params, extra)

  result <- do.call(fn, c(
    list(
      transitions = tx,
      structural_zeros = structural_zeros,
      alternative = alternative,
      n_events = n_events
    ),
    params
  ))

  # After the engine: an engine that rejects structural zeros errors first.
  .warn_structural_zero_counts(tx$obs, structural_zeros)

  .build_lsa_object(
    result = result, data = d,
    engine = engine, params = params, alpha = alpha,
    structural_zeros = structural_zeros, alternative = alternative,
    lag = lag, call = call
  )
}

# Convenience wrappers. Each one pins a specific engine. They live here,
# after lsa(), so the shared help topic is named "lsa".

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

# Structural zeros declare cells that cannot occur. When the data contain
# transitions in such cells, the quasi-independence fit leaves them out
# (its margins are those of the possible cells; Bakeman & Quera, 1995), and
# the fit warns with a classed condition because the data contradict the
# declaration.
.warn_structural_zero_counts <- function(obs, structural_zeros) {
  if (is.null(structural_zeros)) return(invisible(NULL))
  forbidden <- structural_zeros == 0 & obs > 0
  if (!any(forbidden)) return(invisible(NULL))
  cells <- which(forbidden, arr.ind = TRUE)
  labels <- rownames(obs)
  shown <- utils::head(sprintf("%s -> %s (%d)", labels[cells[, 1L]],
                               labels[cells[, 2L]],
                               as.integer(obs[forbidden])), 5L)
  warning(warningCondition(
    sprintf(paste0(
      "%d transitions fall in cells declared structural zeros: %s%s. ",
      "They are left out of the quasi-independence fit (expected counts, ",
      "residuals and the tablewise test use the possible cells only). ",
      "If the cells are possible, do not declare them structural zeros; ",
      "if repeated codes should be merged, collapse them first."),
      as.integer(sum(obs[forbidden])), paste(shown, collapse = ", "),
      if (sum(forbidden) > 5L) ", ..." else ""),
    class = "lagdynamics_structural_zero_counts", call = NULL))
  invisible(NULL)
}

# Resolve the structural-zeros spec to a K x K 0/1 matrix (or NULL when
# nothing is forbidden). `loops = FALSE` forbids every self-transition
# (the diagonal), so callers never hand-build `1 - diag(K)`; an explicit
# 0/1 matrix forbids an arbitrary set of cells. The two combine: with
# loops = FALSE the diagonal of a supplied matrix is also zeroed.
.resolve_structural_zeros <- function(structural_zeros, loops, K, labels) {
  if (is.null(structural_zeros) && isTRUE(loops)) return(NULL)
  m <- if (is.null(structural_zeros)) matrix(1, K, K) else structural_zeros
  if (!is.matrix(m) || nrow(m) != K || ncol(m) != K) {
    .lsa_stop(sprintf("`structural_zeros` must be a %d x %d 0/1 matrix.", K, K),
         call. = FALSE)
  }
  if (!isTRUE(loops)) diag(m) <- 0
  dimnames(m) <- list(labels, labels)
  m
}

# Multi-group fit. Partition the input sequences by `group`, then fit
# the engine once per group on a SHARED global label set so every
# group's K x K matrices index the same states (a group that never
# visits a state still gets a full-size matrix with zeros). Returns a
# named list of `lsa` fits with class "lsa_group" (a named list of
# single-group fits) so the same downstream verbs can dispatch on it.
.lsa_grouped <- function(data, group, lag, engine, alternative, alpha,
                         loops = TRUE, structural_zeros, labels, params,
                         call, ...) {
  d <- lsa_data(data, labels = labels)
  if (!identical(d$source, "events")) {
    .lsa_stop(class = "needs_events", "Grouped lsa() requires event-level sequence data; a ",
         "pre-computed transition matrix cannot be split by group.",
         call. = FALSE)
  }
  # One element per sequence, in original sequence order. Forcing the
  # factor levels to seq_len(n_sequences) avoids the lexicographic
  # reordering that split() would otherwise apply to integer ids.
  per_seq <- split(d$events,
                   factor(d$seq_id, levels = seq_len(d$n_sequences)))
  S <- length(per_seq)
  g <- .resolve_group(group, S)
  levs <- levels(g)
  global_labels <- d$labels

  fits <- lapply(levs, function(lv) {
    idx <- which(g == lv)
    # per_seq[idx] is a list of integer-coded sequences; passing
    # labels = global_labels makes lsa() read those integers as indices
    # into the shared label set rather than re-deriving labels per group.
    lsa(data = per_seq[idx], lag = lag, engine = engine,
        alternative = alternative, alpha = alpha, loops = loops,
        structural_zeros = structural_zeros, labels = global_labels,
        params = params, ...)
  })
  names(fits) <- levs

  structure(
    fits,
    levels = levs,
    group_sizes = as.integer(table(g)[levs]),
    labels = global_labels,
    engine = engine,
    call = call,
    class = c("lsa_group", "list")
  )
}

# Validate and normalise the grouping vector to a factor of length S
# (one entry per sequence). Empty levels are dropped so downstream
# fits and prints never carry a zero-sequence group.
.resolve_group <- function(group, S) {
  if (length(group) != S) {
    .lsa_stop(sprintf(
      paste0("`group` must have one entry per sequence: got length %d ",
             "but the data has %d sequences."),
      length(group), S), call. = FALSE)
  }
  if (anyNA(group)) {
    .lsa_stop("`group` must not contain NA.", call. = FALSE)
  }
  g <- if (is.factor(group)) droplevels(group) else as.factor(group)
  g
}

#' Print, summarise and tidy a lag sequential analysis fit
#'
#' @description
#' `print()`, `summary()`, `as.data.frame()` and `nobs()` describe a fitted lag
#' sequential analysis, an `lsa` object from [lsa()]. A fit has three levels of
#' description: the model as a whole (its size and the tablewise test of
#' independence), the states (how often each is left and entered, and how often
#' a sequence starts in it), and the individual transitions. These methods give
#' the first two levels; [transitions()] gives the third.
#'
#' `print()` shows a compact overview: the engine and lag, the numbers of
#' states, transitions, events and sequences, the likelihood-ratio test of
#' independence, the number of significant transitions out of the cells
#' tested (a structural zero, such as a self-transition with `loops = FALSE`,
#' has no p-value and is not counted), the strongest over-represented
#' transitions and the initial-state distribution. `summary()` returns a `summary.lsa` object
#' holding two tables, a one-row overview of the fit and one row per state,
#' which are read with `as.data.frame(summary(fit), what = )`. For a
#' grouped fit, `summary()` returns one overview row per group.
#' `as.data.frame()` on a fit returns the full edge table, including the
#' integer state ids and the network columns that [transitions()] leaves
#' out. `nobs()` is the number of transitions the model was fitted to.
#'
#' For the per-transition, per-state and per-test tables, prefer the
#' reading verbs [transitions()], [nodes()], [tests()] and [initial()].
#'
#' @param x,object An `lsa` fit from [lsa()], or a grouped fit from
#'   `lsa(..., group = )`. For `print()` and `as.data.frame()` on a
#'   summary, a `summary.lsa` object.
#' @param row.names Optional row names for the returned data frame;
#'   ignored by the grouped method.
#' @param optional Ignored; present for consistency with the generic.
#' @param ... Ignored.
#'
#' @return
#' * `print()` returns `x` invisibly.
#' * `summary()` on an `lsa` fit returns an object of class `summary.lsa`
#'   holding two tables. The overview has one row and the columns
#'   `engine`, `lag`, `n_states`, `n_sequences`, `n_events`,
#'   `n_transitions`, `n_significant`, `alpha`, `lrx2`, `lrx2_df`,
#'   `lrx2_p`, `x2`, `x2_df`, `x2_p` (`NA` where the engine computed no
#'   tablewise test or the input was a transition matrix). The node table
#'   has one row per state and the columns `state`, `outgoing`,
#'   `incoming`, `out_share`, `in_share`, `init_prob`. Read them with
#'   `as.data.frame(s)` and `as.data.frame(s, what = "nodes")`.
#' * `summary()` on a grouped fit returns a `data.frame` with one row per
#'   group: a leading `group` column followed by the overview columns.
#' * `as.data.frame()` on an `lsa` fit returns a `data.frame` with one row
#'   per `from -> to` cell: integer `from` and `to` ids, `from_label`,
#'   `to_label`, the statistic columns of [transitions()], and the network
#'   columns `edge` and `weight`. On a grouped fit the tables are stacked
#'   under a leading `group` column.
#' * `as.data.frame()` on a summary returns the overview or the node table.
#' * `nobs()` returns a single number, the total transition count.
#'
#' @examples
#' fit <- lsa(engagement)
#' fit
#' s <- summary(fit)
#' s
#' as.data.frame(s)
#' as.data.frame(s, what = "nodes")
#' nobs(fit)
#' head(as.data.frame(fit))
#'
#' gfit <- lsa(group_regulation_long, actor = "Actor", action = "Action",
#'             time = "Time", group = "Achiever")
#' gfit
#' summary(gfit)
#'
#' @seealso [lsa()]; [transitions()], [nodes()], [tests()] and [initial()]
#'   for the reading verbs.
#' @name lsa-methods
NULL

#' @rdname lsa-methods
#' @export
print.lsa_group <- function(x, ...) {
  labels <- attr(x, "labels")
  sizes <- attr(x, "group_sizes")
  cat("<lsa_group>\n")
  cat(sprintf("  engine:    %s\n", attr(x, "engine")))
  cat(sprintf("  states:    %d (%s)\n", length(labels),
              paste(utils::head(labels, 10), collapse = ", ")))
  cat(sprintf("  groups:    %d\n", length(x)))
  for (i in seq_along(x)) {
    cat(sprintf("    - %-12s %d sequences\n",
                paste0(names(x)[i], ":"), sizes[i]))
  }
  invisible(x)
}

# Assemble the public S3 object with stable, fully-documented slots.
.build_lsa_object <- function(result, data, engine,
                              params, alpha, structural_zeros,
                              alternative, lag, call) {
  K <- length(data$labels)
  labels <- data$labels

  obs     <- result$obs
  exp_mat <- result$exp
  z       <- result$adj_res
  p       <- result$p
  prob    <- result$prob
  yulesq  <- result$yules_q
  kappa   <- result$kappa
  k_z     <- result$kappa_z
  k_p     <- result$kappa_p

  dimnames(obs)     <- list(labels, labels)
  dimnames(exp_mat) <- list(labels, labels)
  dimnames(z)       <- list(labels, labels)
  dimnames(p)       <- list(labels, labels)
  dimnames(prob)    <- list(labels, labels)
  if (!is.null(yulesq)) dimnames(yulesq) <- list(labels, labels)
  if (!is.null(kappa))  dimnames(kappa)  <- list(labels, labels)
  if (!is.null(k_z))    dimnames(k_z)    <- list(labels, labels)
  if (!is.null(k_p))    dimnames(k_p)    <- list(labels, labels)

  # Column-conditional probabilities P(from | to): the probability that
  # a transition INTO state j came FROM state i. Complements `prob`,
  # which is the row-conditional P(to | from).
  prob_col <- .col_conditional(obs)
  # Pearson tablewise chi-square alongside the LR G^2, on the same model
  # (so it shares lrx2's df). NULL when the engine has no expected table.
  x2 <- .pearson_x2(obs, exp_mat, lrx2 = result$lrx2)

  # Engine-specific extra matrices (anything the engine returned beyond the
  # standard slots, e.g. the two-cell engine's odds_ratio / log_or /
  # log_or_se) -- surfaced both as tidy edge columns and in meta$extra.
  engine_extra <- result[setdiff(
    names(result),
    c("obs", "exp", "prob", "adj_res", "p", "yules_q", "kappa",
      "kappa_z", "kappa_p", "lrx2", "ipf", "structural_zeros",
      "alternative", "n_events_used")
  )]

  edges <- .build_edges(
    obs = obs, exp_mat = exp_mat, prob = prob, prob_col = prob_col,
    z = z, p = p, yulesq = yulesq, kappa = kappa,
    k_z = k_z, k_p = k_p, labels = labels, lag = lag, alpha = alpha,
    extra = engine_extra
  )

  nodes <- .build_nodes(obs = obs, labels = labels)
  inits <- .initial_dist(data)

  is_directed <- !identical(engine, "bidirectional")

  # Slots are grouped by role. The statistical matrices stay flat
  # (fit$obs, fit$adj_res, ...) because that access is idiomatic and
  # used everywhere; the long `edges` table is the same numbers in tidy
  # form. The engine name is stored ONCE as `method` (the name the
  # cograph_network protocol reads) and again only inside `params` as
  # the immutable recipe -- no third standalone `engine` copy. The
  # `lsa_transitions` helper is not stored: it duplicated `obs`/`edges`
  # and nothing downstream read it.
  fit <- list(
    # --- statistical matrices (K x K) ---
    obs = obs,
    exp = exp_mat,
    prob = prob,
    prob_col = prob_col,
    adj_res = z,
    p = p,
    yules_q = yulesq,
    kappa = kappa,
    kappa_z = k_z,
    kappa_p = k_p,
    # --- tablewise tests of independence ---
    lrx2 = result$lrx2,
    x2 = x2,
    # --- cograph_network protocol ---
    weights = obs,
    nodes = nodes,
    edges = edges,
    directed = is_directed,
    method = engine,
    # --- initial-state distribution (NULL for matrix input) ---
    inits = inits,
    # --- provenance ---
    data = data,
    params = list(
      lag = lag,
      engine = engine,
      alternative = alternative,
      alpha = alpha,
      structural_zeros = structural_zeros,
      params = params
    ),
    meta = list(
      source = data$source,
      ipf = result$ipf,
      n_events_used = result$n_events_used,
      package_version = utils::packageVersion("lagdynamics"),
      call = call,
      extra = engine_extra
    )
  )
  class(fit) <- c("lsa", "cograph_network")
  fit
}

.build_edges <- function(obs, exp_mat, prob, prob_col, z, p, yulesq,
                          kappa, k_z, k_p, labels, lag, alpha,
                          extra = NULL) {
  K <- length(labels)
  # `from`/`to` are INTEGER node ids (matching `nodes$id`) to satisfy
  # the cograph_network protocol; human-readable state names live on
  # `from_label`/`to_label` and `edge`.
  grid_id    <- expand.grid(from = seq_len(K), to = seq_len(K),
                            KEEP.OUT.ATTRS = FALSE)
  grid_label <- expand.grid(from = labels, to = labels,
                            KEEP.OUT.ATTRS = FALSE,
                            stringsAsFactors = FALSE)
  vec <- function(m) if (is.null(m)) rep(NA_real_, K * K) else as.vector(m)
  count <- as.vector(obs)
  exp_v <- vec(exp_mat)
  edges <- data.frame(
    from       = grid_id$from,
    to         = grid_id$to,
    from_label = grid_label$from,
    to_label   = grid_label$to,
    lag        = lag,
    count      = count,
    expected   = exp_v,
    prob       = vec(prob),
    prob_col   = vec(prob_col),
    adj_res    = vec(z),
    p          = vec(p),
    yules_q    = vec(yulesq),
    kappa      = vec(kappa),
    kappa_z    = vec(k_z),
    kappa_p    = vec(k_p),
    stringsAsFactors = FALSE,
    row.names  = NULL
  )
  # Engine-specific statistics (e.g. odds_ratio / log_or / log_or_se from the
  # two-cell engine): surface every K x K matrix the engine returned as its
  # own tidy column, so transitions() reports what the engine actually
  # computed instead of hiding it in the fit object.
  for (nm in names(extra)) {
    m <- extra[[nm]]
    if (is.matrix(m) && nrow(m) == K && ncol(m) == K && is.null(edges[[nm]])) {
      edges[[nm]] <- as.vector(m)
    }
  }
  # Derived columns
  edges$lift <- ifelse(is.finite(exp_v) & exp_v > 0, count / exp_v,
                      NA_real_)
  edges$sign <- ifelse(count > exp_v, "over",
                       ifelse(count < exp_v, "under", "expected"))
  edges$significant <- is.finite(edges$p) & edges$p < alpha
  edges$edge <- paste(edges$from_label, edges$to_label, sep = " -> ")
  # `weight` is the cograph_network protocol's edge-weight column,
  # matching the `weights` matrix slot (counts by default). Downstream
  # renderers (cograph::splot edge labels/widths) read edges$weight.
  edges$weight <- count
  edges
}

.build_nodes <- function(obs, labels) {
  K <- length(labels)
  data.frame(
    id = seq_len(K),
    label = labels,
    name = labels,
    outgoing = as.numeric(rowSums(obs)),
    incoming = as.numeric(colSums(obs)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# Initial-state distribution: the proportion of sequences that start in
# each state (named, sums to 1). Returns NULL when the fit came from a
# transition matrix (no sequences, so no notion of a starting state).
.initial_dist <- function(data) {
  if (!identical(data$source, "events") || is.null(data$events)) {
    return(NULL)
  }
  per <- split(data$events,
               factor(data$seq_id, levels = seq_len(data$n_sequences)))
  first <- vapply(per, function(s) s[1L], integer(1L))
  counts <- tabulate(first, nbins = data$n_states)
  inits <- counts / sum(counts)
  names(inits) <- data$labels
  inits
}

# Column-conditional probabilities P(from | to) = obs[i, j] / sum_i
# obs[i, j]: of the transitions that arrive in target j, what fraction
# came from source i. Columns with no incoming transitions are NA.
.col_conditional <- function(obs) {
  cs <- colSums(obs)
  pc <- sweep(obs, 2L, cs, "/")
  if (any(cs == 0)) pc[, cs == 0] <- NA_real_
  pc
}

# Pearson tablewise chi-square on the same independence model as the LR
# G^2, so it shares lrx2's degrees of freedom (correct under structural
# zeros / quasi-independence). NULL when the engine has no expected
# table or no LR test to borrow the df from.
.pearson_x2 <- function(obs, exp_mat, lrx2) {
  if (is.null(exp_mat) || is.null(lrx2)) return(NULL)
  ok <- is.finite(exp_mat) & exp_mat > 0
  stat <- sum((obs[ok] - exp_mat[ok])^2 / exp_mat[ok])
  df <- lrx2$df
  p <- if (is.finite(df) && df > 0) {
    stats::pchisq(stat, df = df, lower.tail = FALSE)
  } else {
    NA_real_
  }
  list(statistic = stat, df = df, p = p)
}
