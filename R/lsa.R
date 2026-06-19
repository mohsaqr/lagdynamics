# Main constructor. Reads canonical lsa_data + lsa_transitions, dispatches
# to the named engine via the registry, and assembles the S3 fit object
# with stable named slots.

#' Lag Sequential Analysis
#'
#' Fits a lag sequential analysis (LSA) on categorical event sequence
#' data using a registered engine. Returns a tidy S3 object with named
#' slots for observed/expected/probability/residual matrices and a
#' long-format edge table suitable for transition-network visualization.
#'
#' @param data Sequence input (any form accepted by [lsa_data()]),
#'   *or* a raw long-format event-log `data.frame` when the `actor` /
#'   `action` arguments are supplied (see below). Accepted already-
#'   sequenced forms include vectors, lists of sequences, wide
#'   matrices/data.frames, transition-count matrices, and sequence-
#'   bearing objects (`tna`, `group_tna`, `nestimate_data`, `stslist`).
#'   `NA` and empty-string cells are treated as missingness, not as a
#'   state: they are dropped wherever they occur and no transition is
#'   counted into or out of them. To model missingness as its own
#'   state, recode it (e.g. `NA -> "missing"`) before calling `lsa()`.
#' @param lag Integer. The transition lag. Default `1`. Positive lags
#'   count successors (state at `t -> t + lag`); negative lags count
#'   predecessors (what occurred `|lag|` steps before); `0` pairs each
#'   event with itself (degenerate for single-stream event data --
#'   genuine co-occurrence needs concurrent codes, not yet supported).
#'   Pre-computed transition-matrix input supports `lag = 1` only. To
#'   analyse several lags at once, see [lsa_lags()].
#' @param engine Character scalar. The engine name, registered via
#'   [register_lsa_engine()]. Built-in engines: `"classical"`,
#'   `"two_cell"`, `"bidirectional"`, `"parallel_dominance"`,
#'   `"nonparallel_dominance"`. Default `"classical"`.
#' @param alternative Character scalar. The alternative hypothesis for
#'   adjusted-residual and kappa p-values: one of `"two.sided"`
#'   (default), `"greater"`, or `"less"`.
#' @param alpha Numeric. Significance threshold used to mark edges as
#'   significant in `fit$edges$significant`. Default `0.05`.
#' @param loops Logical. Keep self-transitions (the diagonal)? **Default
#'   `TRUE`.** Set `loops = FALSE` to forbid every self-transition -- the
#'   common reason to exclude cells -- without building a matrix by hand.
#' @param structural_zeros Optional `K x K` 0/1 matrix for an *arbitrary*
#'   forbidden-cell pattern, where `0` marks a forbidden (structural-zero)
#'   cell and `1` an estimable one. **Default `NULL`: every cell is part
#'   of the model.** Combines with `loops`: `loops = FALSE` also zeros the
#'   diagonal of a supplied matrix. When any cell is forbidden the engine
#'   switches to iterative proportional fitting and Christensen's
#'   design-matrix residuals (see `inst/REFERENCES.md` §2.2, §4.2).
#' @param labels Optional character vector of state labels.
#' @param group Optional grouping for a multi-group fit. Either a vector
#'   with one entry per input sequence (length `n_sequences`), or --- for
#'   **long-format** input (see `actor`/`action`) --- the **name of a
#'   grouping column** in the log, which must be constant within each
#'   actor/session so each recovered sequence maps to one group. Sequences
#'   are partitioned by group and a separate `lsa` fit is built for each.
#'   All group fits share one global label set (derived from the full
#'   data) so their `K x K` matrices are directly comparable, even when
#'   a group never visits some state. Returns an `lsa_group` object (a
#'   named list of `lsa` fits). Requires event-level input; a
#'   pre-computed transition matrix cannot be split by group. Default
#'   `NULL` (single-group fit).
#' @param actor,action,time,order,session Column names (each a single
#'   string) for **long-format** event-log input. Supplying `action`
#'   (and `actor`) switches `lsa()` into long-format mode: the raw log
#'   in `data` is sequenced into event sequences by grouping rows per
#'   `actor` (optionally crossed with an explicit `session` id),
#'   ordering within each group by `order` if given else by `time`, and
#'   -- when `time` is given and no `session` column is -- starting a
#'   new session whenever the gap between consecutive events exceeds
#'   `time_threshold` seconds. All `NULL` by default (input is taken
#'   as already-sequenced). Cannot be combined with `group`.
#' @param time_threshold Numeric. Maximum gap in seconds between
#'   consecutive events before a new session is started in long-format
#'   mode. Default `900` (15 minutes). Ignored unless `time` is given
#'   and `session` is not.
#' @param custom_format Optional `strptime` format string for parsing
#'   the `time` column (e.g. `"%Y-%m-%d %H:%M:%S"`). Default `NULL`
#'   (native date/time classes and ISO strings are parsed directly).
#' @param is_unix_time Logical. Treat the `time` column as a Unix
#'   epoch. Default `FALSE`.
#' @param unix_time_unit Character. Unit of the Unix epoch when
#'   `is_unix_time = TRUE`: `"seconds"` (default), `"milliseconds"`, or
#'   `"microseconds"`.
#' @param params Optional named list of engine-specific parameters
#'   forwarded to the engine function.
#' @param ... Additional engine-specific parameters (merged into
#'   `params`).
#'
#' @return An object of class `c("lsa", "cograph_network")`. Read it with
#' the verbs rather than by reaching into slots: [transitions()] for the
#' tidy edge table, [nodes()], [tests()], [initial()], and [summary()]
#' for the other results, and [plot()]/[plot_transitions()] to draw it.
#' Every number a verb returns is backed by these slots:
#' \describe{
#'   \item{edges}{The tidy one-row-per-transition frame that backs
#'     [transitions()] (with extra `cograph_network` protocol columns).}
#'   \item{nodes}{Data frame backing [nodes()]: `id, label, name,
#'     outgoing, incoming`.}
#'   \item{obs, exp, prob, prob_col, adj_res, p, yules_q, kappa,
#'     kappa_z, kappa_p}{The same per-cell quantities as `edges`, in
#'     `K x K` matrix form (`prob` is row-conditional P(to | from),
#'     `prob_col` column-conditional P(from | to)). Convenient for
#'     matrix algebra; not the primary interface.}
#'   \item{lrx2, x2}{Lists `(statistic, df, p)` backing [tests()]: the
#'     tablewise likelihood-ratio (G^2) and Pearson chi-square tests of
#'     independence; `NULL` for engines without an expected table.}
#'   \item{inits}{Named numeric vector backing [initial()] (proportion
#'     of sequences starting in each state, sums to 1); `NULL` for
#'     transition-matrix input.}
#'   \item{weights}{`K x K` matrix used as the default edge weight for
#'     plotting. Equal to `obs` (counts) by default.}
#'   \item{directed}{Logical scalar; `TRUE` for directed engines,
#'     `FALSE` for `bidirectional`.}
#'   \item{method}{Engine name (the slot the `cograph_network` protocol
#'     reads). Also recorded in `params$engine`.}
#'   \item{data}{The canonical `lsa_data` object (events + seq_id).}
#'   \item{params}{Immutable snapshot of all parameters used (recipe),
#'     including `params$engine`.}
#'   \item{meta}{List with source, IPF info, version, and call.}
#' }
#'
#' When `group` is supplied, returns an object of class
#' `c("lsa_group", "list")`: a named list of `lsa` fits (one per group
#' level) carrying `levels`, `group_sizes`, `labels`, and `engine`
#' attributes. Downstream verbs ([lsa_to_tna()],
#' [transitions()], [reliability_lsa()], etc.) dispatch on
#' it and return grouped results.
#'
#' @examples
#' seq <- c("Question", "Explain", "Agree",
#'          "Question", "Explain", "Elaborate",
#'          "Agree", "Question", "Explain")
#' fit <- lsa(seq, engine = "classical")
#' fit
#' head(fit$edges)
#'
#' @seealso [lsa_data()], [lsa_transitions()], [register_lsa_engine()],
#'   [list_lsa_engines()]
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
    is.numeric(alpha), length(alpha) == 1L, alpha > 0, alpha < 1
  )

  # Long-format mode: a raw event log is sequenced into a list of event
  # sequences before any analysis. Triggered by supplying `action`.
  if (!is.null(action)) {
    if (is.null(actor)) {
      stop("Long-format sequencing needs both `actor` and `action` ",
           "column names.", call. = FALSE)
    }
    # With long-format input, `group` (if given) is the NAME of a grouping
    # column in the log; .prepare_long() derives one label per recovered
    # sequence and returns them as an attribute, which then drives the
    # grouped fit below -- no manual split-then-relabel ritual.
    long_group <- NULL
    if (!is.null(group)) {
      if (!is.character(group) || length(group) != 1L) {
        stop("With long-format input, `group` must be the name of a ",
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

  .build_lsa_object(
    result = result, data = d,
    engine = engine, params = params, alpha = alpha,
    structural_zeros = structural_zeros, alternative = alternative,
    lag = lag, call = call
  )
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
    stop(sprintf("`structural_zeros` must be a %d x %d 0/1 matrix.", K, K),
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
    stop("Grouped lsa() requires event-level sequence data; a ",
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
    stop(sprintf(
      paste0("`group` must have one entry per sequence: got length %d ",
             "but the data has %d sequences."),
      length(group), S), call. = FALSE)
  }
  if (anyNA(group)) {
    stop("`group` must not contain NA.", call. = FALSE)
  }
  g <- if (is.factor(group)) droplevels(group) else as.factor(group)
  g
}

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
      package_version = utils::packageVersion("lagseq"),
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
