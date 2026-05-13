# Main constructor. Reads canonical lsa_data + lsa_transitions, dispatches
# to the named engine via the registry, and assembles the S3 fit object
# with stable named slots in the style of Nestimate's netobject.

#' Lag Sequential Analysis
#'
#' Fits a lag sequential analysis (LSA) on categorical event sequence
#' data using a registered engine. Returns a tidy S3 object with named
#' slots for observed/expected/probability/residual matrices and a
#' long-format edge table suitable for transition-network visualization.
#'
#' @param data Sequence input (any form accepted by [lsa_data()]).
#' @param lag Positive integer. The transition lag. Default `1`.
#' @param engine Character scalar. The engine name, registered via
#'   [register_lsa_engine()]. Built-in engines: `"classical"`,
#'   `"two_cell"`, `"bidirectional"`, `"parallel_dominance"`,
#'   `"nonparallel_dominance"`. Default `"classical"`.
#' @param alternative Character scalar. The alternative hypothesis for
#'   adjusted-residual and kappa p-values: one of `"two.sided"`
#'   (default), `"greater"`, or `"less"`.
#' @param alpha Numeric. Significance threshold used to mark edges as
#'   significant in `fit$edges$significant`. Default `0.05`.
#' @param structural_zeros Optional `K x K` 0/1 matrix. A `0` marks a
#'   cell as a structural zero (forbidden transition). When supplied,
#'   the engine uses iterative proportional fitting and Christensen's
#'   design-matrix residuals (see `inst/REFERENCES.md` §2.2, §4.2).
#'   A common pattern is `1 - diag(K)` to forbid self-transitions.
#' @param labels Optional character vector of state labels.
#' @param params Optional named list of engine-specific parameters
#'   forwarded to the engine function.
#' @param ... Additional engine-specific parameters (merged into
#'   `params`).
#'
#' @return An object of class `c("lsa", "cograph_network")` with
#' elements:
#' \describe{
#'   \item{obs, exp, prob, adj_res, p, yules_q, kappa, kappa_z, kappa_p}{
#'     `K x K` matrices.}
#'   \item{lrx2}{List `(statistic, df, p)` from the tablewise LR test.}
#'   \item{weights}{`K x K` matrix used as the default edge weight for
#'     plotting. Equal to `obs` (counts) by default.}
#'   \item{nodes}{Data frame: `id, label, name, outgoing, incoming`.}
#'   \item{edges}{Tidy edge frame.}
#'   \item{data}{The canonical `lsa_data` object.}
#'   \item{transitions}{The `lsa_transitions` object.}
#'   \item{engine}{Engine name used.}
#'   \item{params}{Immutable snapshot of all parameters used (recipe).}
#'   \item{directed}{Logical scalar; `TRUE` for directed engines,
#'     `FALSE` for `bidirectional`.}
#'   \item{method}{Equal to `engine` for cograph compatibility.}
#'   \item{meta}{List with engine info, version, and call.}
#' }
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
                structural_zeros = NULL,
                labels = NULL,
                params = list(),
                ...) {
  call <- match.call()
  alternative <- match.arg(alternative)
  stopifnot(
    is.numeric(alpha), length(alpha) == 1L, alpha > 0, alpha < 1
  )

  d  <- lsa_data(data, labels = labels)
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
    result = result, data = d, transitions = tx,
    engine = engine, params = params, alpha = alpha,
    structural_zeros = structural_zeros, alternative = alternative,
    lag = lag, call = call
  )
}

# Assemble the public S3 object with stable, fully-documented slots.
.build_lsa_object <- function(result, data, transitions, engine,
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

  edges <- .build_edges(
    obs = obs, exp_mat = exp_mat, prob = prob,
    z = z, p = p, yulesq = yulesq, kappa = kappa,
    k_z = k_z, k_p = k_p, labels = labels, lag = lag, alpha = alpha
  )

  nodes <- .build_nodes(obs = obs, labels = labels)

  is_directed <- !identical(engine, "bidirectional")

  fit <- list(
    obs = obs,
    exp = exp_mat,
    prob = prob,
    adj_res = z,
    p = p,
    yules_q = yulesq,
    kappa = kappa,
    kappa_z = k_z,
    kappa_p = k_p,
    lrx2 = result$lrx2,
    weights = obs,
    nodes = nodes,
    edges = edges,
    data = data,
    transitions = transitions,
    directed = is_directed,
    method = engine,
    engine = engine,
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
      extra = result[setdiff(
        names(result),
        c("obs", "exp", "prob", "adj_res", "p", "yules_q", "kappa",
          "kappa_z", "kappa_p", "lrx2", "ipf", "structural_zeros",
          "alternative", "n_events_used")
      )]
    )
  )
  class(fit) <- c("lsa", "cograph_network")
  fit
}

.build_edges <- function(obs, exp_mat, prob, z, p, yulesq, kappa,
                          k_z, k_p, labels, lag, alpha) {
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
    adj_res    = vec(z),
    p          = vec(p),
    yules_q    = vec(yulesq),
    kappa      = vec(kappa),
    kappa_z    = vec(k_z),
    kappa_p    = vec(k_p),
    stringsAsFactors = FALSE,
    row.names  = NULL
  )
  # Derived columns
  edges$lift <- ifelse(is.finite(exp_v) & exp_v > 0, count / exp_v,
                      NA_real_)
  edges$sign <- ifelse(count > exp_v, "over",
                       ifelse(count < exp_v, "under", "expected"))
  edges$significant <- is.finite(edges$p) & edges$p < alpha
  edges$edge <- paste(edges$from_label, edges$to_label, sep = " -> ")
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
