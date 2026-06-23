#' @title Basis Profile Curve (\verb{BPC}) identification
#' @description
#' Identifies a small set of canonical temporal response shapes - \emph{basis
#' profile curves} (\verb{BPCs}) - from the single-trial stimulation-evoked
#' responses recorded at one measurement electrode, where the trials are grouped
#' by stimulation site (or any other condition). Each stimulation group is
#' assigned to the \verb{BPC} that best explains its trials, and the projection
#' strength of every group is quantified. This is the across-stimulation-site
#' \verb{BPC} method (see \sQuote{References}); the related
#' \code{\link{crp_cluster}} applies the same idea across electrodes.
#'
#' @details
#' Let \eqn{V} be the (windowed) \eqn{T \times K} matrix of all single trials and
#' let \code{groups} map each of the \eqn{K} columns to one of \eqn{n} stimulation
#' subgroups.
#' \enumerate{
#' \item \strong{Window.} When \code{time} and \code{time_window} are supplied the
#'   rows of \code{x} are cropped to that window; otherwise all rows are used.
#'   Rows with non-finite values are dropped.
#' \item \strong{Internal projections.} Each trial is \eqn{L_2}-normalized into
#'   \eqn{V_0}, and \eqn{P = V_0^\top V} collects the projection of every native
#'   trial onto every normalized trial.
#' \item \strong{Significance matrix.} For each ordered pair of subgroups
#'   \eqn{(k, l)} the set of relevant entries of \eqn{P} is gathered (the
#'   off-diagonal of the within-group block when \eqn{k = l}, the whole cross
#'   block otherwise) and reduced to a one-sample t-statistic versus zero. These
#'   form the \eqn{n \times n} significance matrix \eqn{\Xi}.
#' \item \strong{Rank selection.} \eqn{\Xi} is made non-negative and rescaled, then
#'   factorized with \code{\link{naive_nmf}} at decreasing inner rank while the
#'   degeneracy score \eqn{\zeta} - the sum of the upper off-diagonal of the
#'   row-normalized \eqn{HH^\top} - exceeds \code{zeta_threshold}.
#' \item \strong{Assignment.} Each subgroup takes its winner-take-all \verb{BPC}
#'   over the normalized \verb{NMF} loadings; with \code{null_class}, loadings
#'   below \eqn{1/(2\sqrt{n})} are left unassigned (the \dQuote{null} class).
#' \item \strong{Basis curves.} Per \verb{BPC}, the first linear-kernel \verb{PCA}
#'   component of all trials in its member subgroups, sign-oriented to a positive
#'   mean projection.
#' \item \strong{Weights.} For each member subgroup the per-trial coefficient
#'   \eqn{\alpha} (projection onto the basis curve) is normalized by the residual
#'   magnitude; the mean is the projection weight and a one-sample t-test gives a
#'   significance value.
#' }
#'
#' @param x numeric matrix of single-trial evoked voltages with shape
#' \code{time x trials} (rows are timepoints, columns are trials). All trials
#' across all stimulation groups are stacked column-wise.
#' @param groups vector of length \code{ncol(x)} (integer, character or factor)
#' giving the stimulation subgroup of each column of \code{x}. At least two
#' subgroups are required and each subgroup needs at least two trials.
#' @param time optional numeric vector of length \code{nrow(x)} giving the
#' stimulus-aligned time (in seconds) of each row of \code{x}. Only used to crop
#' to \code{time_window} and as the plotting axis; the method itself does not
#' require it.
#' @param time_window optional numeric \code{c(lo, hi)} analysis window (seconds);
#' requires \code{time}. An \code{NA} extends to the corresponding data edge.
#' \code{NULL} (the default) uses all rows.
#' @param n_bpc integer or \code{NULL}; if given, fixes the number of basis
#' curves and skips automatic rank selection.
#' @param initial_rank integer or \code{NULL}; the starting \verb{NMF} rank for
#' automatic selection. \code{NULL} uses \code{min(n - 1L, 10L)} where \eqn{n} is
#' the number of subgroups.
#' @param zeta_threshold numeric \eqn{> 0}; rank-selection cutoff on \eqn{\zeta}
#' (see \sQuote{Details}); smaller values give fewer curves. Defaults to \code{1},
#' matching the \verb{BPC} reference.
#' @param null_class logical; if \code{TRUE} (the default) subgroups whose
#' winner-take-all loading falls below \eqn{1/(2\sqrt{n})} are left unassigned
#' (\code{NA}); if \code{FALSE} every subgroup is forced into a curve.
#' @param nmf_max_iters,nmf_tol passed to \code{\link{naive_nmf}}.
#' @param verbose logical; whether to report progress.
#'
#' @returns A named list of class \code{ravetools_bpc}:
#' \describe{
#' \item{\code{curves}}{Numeric matrix, time \eqn{\times} number of \verb{BPCs};
#'   column \eqn{q} is the basis profile curve \eqn{B_q(t)}, the first
#'   linear-kernel \verb{PCA} component of its member trials (unit-norm,
#'   sign-oriented to a positive mean projection).}
#' \item{\code{time}}{Numeric vector, the (windowed) time axis for \code{curves};
#'   the sample index when \code{time} was not supplied.}
#' \item{\code{group_labels}}{The unique stimulation subgroup labels, in the order
#'   used by \code{clusters} and \code{xi}.}
#' \item{\code{clusters}}{Integer vector, the \verb{BPC} index assigned to each
#'   subgroup; \code{NA} for subgroups left in the null class.}
#' \item{\code{excluded_groups}}{The labels of subgroups not represented by any
#'   \verb{BPC}.}
#' \item{\code{weights}}{A \code{data.frame} with one row per (\verb{BPC},
#'   subgroup) membership: \code{bpc}, \code{group}, \code{n_trials},
#'   \code{weight} (mean residual-normalized \eqn{\alpha}) and \code{p_value}.}
#' \item{\code{alpha}}{Numeric matrix, trials \eqn{\times} number of \verb{BPCs};
#'   the per-trial projection \eqn{\alpha} of every trial onto each basis curve.}
#' \item{\code{xi}}{Numeric matrix, the \eqn{n \times n} significance matrix
#'   \eqn{\Xi}.}
#' \item{\code{nmf}}{The \code{\link{naive_nmf}} result for the selected rank, with
#'   \code{H} (raw loadings) and \code{H0} (winner-take-all, thresholded).}
#' \item{\code{n_bpc}}{Integer, the number of basis curves found.}
#' \item{\code{time_window}}{The effective analysis window used (or \code{NULL}).}
#' }
#'
#' @references
#' The \verb{BPC} method is described in \doi{10.1371/journal.pcbi.1008710}, with
#' a reference Python implementation at
#' \url{https://github.com/MultimodalNeuroimagingLab/bpc_jupyter}.
#'
#' @seealso \code{\link{crp}}, \code{\link{crp_cluster}}, \code{\link{naive_nmf}}
#'
#' @examples
#' # Three response shapes, several stimulation groups per shape.
#' \donttest{
#' set.seed(1)
#' n_time <- 300L
#' tt <- seq(-0.2, 1, length.out = n_time)
#' shapes <- list(
#'   exp(-((tt - 0.08) / 0.03)^2) - 0.5 * exp(-((tt - 0.18) / 0.04)^2),
#'   exp(-((tt - 0.38) / 0.03)^2) - 0.5 * exp(-((tt - 0.50) / 0.04)^2),
#'   exp(-((tt - 0.70) / 0.04)^2)
#' )
#'
#' # 3 stimulation groups per shape, 12 trials each
#' V <- NULL
#' groups <- NULL
#' g <- 0L
#' for (s in seq_along(shapes)) {
#'   for (rep in seq_len(3L)) {
#'     g <- g + 1L
#'     trials <- outer(shapes[[s]], runif(12L, 0.5, 1.5)) +
#'       matrix(rnorm(n_time * 12L, sd = 0.2), n_time, 12L)
#'     V <- cbind(V, trials)
#'     groups <- c(groups, rep(g, 12L))
#'   }
#' }
#'
#' res <- bpc(V, groups, time = tt, time_window = c(0, 1), verbose = TRUE)
#' res$n_bpc
#' res$clusters
#' plot(res)
#' }
#'
#' @export
bpc <- function(
  x,
  groups,
  time = NULL,
  time_window = NULL,
  n_bpc = NULL,
  initial_rank = NULL,
  zeta_threshold = 1,
  null_class = TRUE,
  nmf_max_iters = 1e4,
  nmf_tol = c(1e-4, 1e-8),
  verbose = TRUE
) {

  # ---- validate inputs ----------------------------------------------------
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric matrix with shape time x trials.")
  }
  if (ncol(x) < 2L) {
    stop("`x` must have at least 2 trials (columns).")
  }
  if (length(groups) != ncol(x)) {
    stop("`groups` must have length `ncol(x)` (one label per trial/column).")
  }

  # ---- crop rows to the analysis window (optional) ------------------------
  if (!is.null(time)) {
    time <- as.numeric(time)
    if (length(time) != nrow(x)) {
      stop("`time` must have length `nrow(x)`.")
    }
  }
  if (!is.null(time_window)) {
    if (is.null(time)) {
      stop("`time_window` requires `time` to be supplied.")
    }
    time_window <- as.numeric(time_window)
    if (length(time_window) == 1L) {
      time_window <- c(0, time_window)
    } else {
      time_window <- time_window[c(1L, 2L)]
    }
    rng <- range(time)
    time_window[is.na(time_window)] <- rng[is.na(time_window)]
    if (is.unsorted(time_window)) {
      stop("`time_window` must be a numeric c(lo, hi) with lo < hi.")
    }
    row_sel <- time >= time_window[1L] & time <= time_window[2L]
  } else {
    row_sel <- rep(TRUE, nrow(x))
  }

  V <- x[row_sel, , drop = FALSE]
  axis <- if (is.null(time)) seq_len(nrow(x))[row_sel] else time[row_sel]

  # drop rows with any non-finite trial value
  keep <- rowSums(!is.finite(V)) == 0L
  V <- V[keep, , drop = FALSE]
  axis <- axis[keep]
  if (nrow(V) < 2L) {
    stop("Fewer than 2 finite timepoints remain after windowing.")
  }

  # ---- subgroup index lists -----------------------------------------------
  group_labels <- sort(unique(groups))
  n <- length(group_labels)
  if (n < 2L) {
    stop("`groups` must contain at least 2 distinct subgroups.")
  }
  idx_list <- lapply(group_labels, function(g) which(groups == g))
  group_sizes <- lengths(idx_list)
  if (any(group_sizes < 2L)) {
    small <- group_labels[group_sizes < 2L]
    stop("Every subgroup needs at least 2 trials; offending subgroup(s): ",
         paste(small, collapse = ", "), ".")
  }

  # ---- normalized internal projections P = V0^T V -------------------------
  norms <- sqrt(colSums(V * V))
  norms[norms == 0] <- 1
  V0 <- t(t(V) / norms)
  V0[!is.finite(V0)] <- 0
  P <- crossprod(V0, V)               # trials x trials

  # ---- significance matrix Xi ---------------------------------------------
  xi <- bpc_significance_matrix(P, idx_list)

  # ---- NMF on the non-negative, rescaled Xi -------------------------------
  t0 <- xi
  t0[t0 < 0 | !is.finite(t0)] <- 0
  mx <- max(t0)
  if (mx > 0) t0 <- t0 / mx

  if (is.null(initial_rank)) {
    initial_rank <- min(n - 1L, 10L)
  }
  initial_rank <- max(1L, min(as.integer(initial_rank), n - 1L))
  if (!is.null(n_bpc)) {
    n_bpc <- max(1L, min(as.integer(n_bpc), n - 1L))
  }

  sel <- crp_cluster_nmf_select(
    t0,
    n_clusters = n_bpc,
    initial_rank = initial_rank,
    zeta_threshold = zeta_threshold,
    nmf_max_iters = nmf_max_iters,
    nmf_tol = nmf_tol,
    verbose = verbose
  )
  nmf <- sel$nmf

  # ---- winner-take-all assignment + null class ----------------------------
  # Subgroups join the BPC with the largest row-normalized NMF loading; with
  # `null_class` a loading below 1/(2*sqrt(n)) - half the uniform loading
  # 1/sqrt(n) - leaves the subgroup unassigned (NA, the null class).
  Hn <- nmf$H
  rn <- sqrt(rowSums(Hn * Hn))
  rn[rn == 0] <- 1
  Hn <- Hn / rn
  raw_cluster <- apply(Hn, 2L, which.max)
  if (isTRUE(null_class)) {
    winner <- Hn[cbind(raw_cluster, seq_len(n))]
    raw_cluster[winner <= 1 / (2 * sqrt(n))] <- NA_integer_
  }
  present <- sort(unique(raw_cluster))
  clusters <- match(raw_cluster, present)
  n_found <- length(present)

  # thresholded winner-take-all loadings (rows = BPC, cols = subgroup)
  H0 <- matrix(0, nrow = nrow(Hn), ncol = n)
  assigned <- !is.na(raw_cluster)
  if (any(assigned)) {
    H0[cbind(raw_cluster[assigned], which(assigned))] <-
      Hn[cbind(raw_cluster[assigned], which(assigned))]
  }

  # ---- basis profile curves -----------------------------------------------
  basis_curves <- bpc_basis(V, clusters, idx_list, n_found)

  # ---- per-trial projection weights and significance ----------------------
  stat <- bpc_curve_statistics(V, basis_curves, clusters, idx_list, group_labels)

  excluded_groups <- group_labels[is.na(clusters)]

  rownames(basis_curves) <- NULL
  dimnames(xi) <- list(as.character(group_labels), as.character(group_labels))

  nmf$H0 <- H0

  structure(
    class = "ravetools_bpc",
    list(
      curves = basis_curves,
      time = axis,
      group_labels = group_labels,
      clusters = clusters,
      excluded_groups = excluded_groups,
      weights = stat$weights,
      alpha = stat$alpha,
      xi = xi,
      nmf = nmf,
      n_bpc = n_found,
      time_window = if (is.null(time_window)) NULL else time_window
    )
  )
}

#' @title Plot Basis Profile Curve results
#' @description S3 plot method for \code{ravetools_bpc} objects returned by
#' \code{\link{bpc}}. Draws a three-panel figure: the basis profile curves, the
#' subgroup-by-subgroup significance matrix \eqn{\Xi} ordered by assignment, and
#' the per-subgroup projection weights.
#' @param x a \code{ravetools_bpc} object.
#' @param ... ignored.
#' @return Invisibly returns \code{x}.
#' @export
plot.ravetools_bpc <- function(x, ...) {
  Q <- x$n_bpc
  cols <- grDevices::hcl.colors(max(Q, 2L), palette = "Dark 3")

  op <- graphics::par(mfrow = c(1L, 3L))
  on.exit(graphics::par(op), add = TRUE)

  # ---- Panel 1: basis profile curves --------------------------------------
  graphics::matplot(
    x$time, x$curves, type = "l", lty = 1L, lwd = 2L, col = cols[seq_len(Q)],
    las = 1L, xlab = "Time (s)", ylab = "Basis amplitude (a.u.)",
    main = "Basis profile curves"
  )
  graphics::abline(h = 0, col = "#cccccc", lty = 3L)
  graphics::legend(
    "topright", bty = "n", cex = 0.8, lwd = 2L, col = cols[seq_len(Q)],
    legend = sprintf("BPC %d (n=%d)", seq_len(Q),
                     tabulate(x$clusters, nbins = Q))
  )

  # ---- Panel 2: significance matrix ordered by assignment -----------------
  ord <- order(x$clusters, na.last = TRUE)
  S <- x$xi[ord, ord, drop = FALSE]
  graphics::image(
    seq_len(nrow(S)), seq_len(ncol(S)), S,
    col = grDevices::hcl.colors(64L, palette = "YlOrRd", rev = TRUE),
    xlab = "Subgroup (ordered by BPC)", ylab = "Subgroup",
    main = expression("Significance matrix " * Xi), las = 1L
  )
  brk <- cumsum(tabulate(x$clusters, nbins = Q))
  brk <- brk[brk > 0 & brk < nrow(S)] + 0.5
  graphics::abline(v = brk, h = brk, col = "#202020", lty = 1L, lwd = 1L)

  # ---- Panel 3: per-subgroup projection weights ---------------------------
  w <- x$weights
  if (nrow(w)) {
    w <- w[order(w$bpc), , drop = FALSE]
    graphics::barplot(
      w$weight, col = cols[w$bpc], border = NA, las = 1L,
      names.arg = w$group, cex.names = 0.7,
      ylab = "Projection weight", main = "Projection weights"
    )
  } else {
    graphics::plot.new()
    graphics::title(main = "Projection weights")
  }

  invisible(x)
}


# --------------------------------------------------------------------------
# Internal helpers (not exported)
# --------------------------------------------------------------------------

# Significance matrix Xi (BPC `nativeNormalized`). `P` is the trials x trials
# projection matrix (normalized x native); `idx_list` gives the trial columns of
# each subgroup. Entry (k, l) is the one-sample t-statistic versus zero of the
# off-diagonal within-group block (k == l) or the full cross block (k != l). The
# resulting n x n matrix is generally asymmetric.
bpc_significance_matrix <- function(P, idx_list) {
  n <- length(idx_list)
  xi <- matrix(0, nrow = n, ncol = n)
  for (k in seq_len(n)) {
    ik <- idx_list[[k]]
    for (l in seq_len(n)) {
      il <- idx_list[[l]]
      if (k == l) {
        block <- P[ik, ik, drop = FALSE]
        b <- block[row(block) != col(block)]   # all off-diagonal elements
      } else {
        b <- as.numeric(P[ik, il, drop = FALSE])
      }
      b <- b[is.finite(b)]
      if (length(b) < 2L) next
      s <- stats::sd(b)
      if (!is.finite(s) || s == 0) next
      xi[k, l] <- mean(b) / (s / sqrt(length(b)))
    }
  }
  xi
}

# Per-BPC basis profile curve: the first linear-kernel PCA component of all
# trials in the BPC's member subgroups, sign-oriented to a positive mean
# projection. Reuses `crp_kt_pca`; mirrors `crp_cluster_basis`.
bpc_basis <- function(V, clusters, idx_list, n_found) {
  Tn <- nrow(V)
  basis_curves <- matrix(0, nrow = Tn, ncol = n_found)
  for (q in seq_len(n_found)) {
    members <- which(clusters == q)
    cols <- unlist(idx_list[members], use.names = FALSE)
    Vq <- V[, cols, drop = FALSE]
    pca <- crp_kt_pca(Vq)
    Bq <- pca$E[, 1L]
    if (mean(crossprod(Bq, Vq)) < 0) Bq <- -Bq
    basis_curves[, q] <- Bq
  }
  basis_curves
}

# Per-trial coefficients and projection weights (BPC `curvesStatistics` +
# `projectionWeights`). For each basis curve B: alpha = B^T V (one value per
# trial), residual ep = V - B alpha^T, and the per-trial residual power is the
# diagonal of ep^T ep. For each member subgroup the residual-normalized
# coefficients alpha / sqrt(ep2) give a mean projection weight and a one-sample
# t-test p-value.
bpc_curve_statistics <- function(V, basis_curves, clusters, idx_list, group_labels) {
  n_found <- ncol(basis_curves)
  K <- ncol(V)
  alpha <- matrix(0, nrow = K, ncol = n_found)
  rows <- list()
  for (q in seq_len(n_found)) {
    B <- basis_curves[, q]
    al <- as.numeric(crossprod(B, V))            # length K
    alpha[, q] <- al
    ep <- V - outer(B, al)
    ep2 <- colSums(ep * ep)                      # per-trial residual power
    members <- which(clusters == q)
    for (g in members) {
      idx <- idx_list[[g]]
      ratio <- al[idx] / sqrt(ep2[idx])
      ratio <- ratio[is.finite(ratio)]
      pval <- NA_real_
      if (length(ratio) >= 2L) {
        tt <- tryCatch(stats::t.test(ratio, mu = 0), error = function(e) NULL)
        if (!is.null(tt)) pval <- tt$p.value
      }
      rows[[length(rows) + 1L]] <- data.frame(
        bpc = q,
        group = group_labels[g],
        n_trials = length(idx),
        weight = if (length(ratio)) mean(ratio) else NA_real_,
        p_value = pval,
        stringsAsFactors = FALSE
      )
    }
  }
  weights <- if (length(rows)) {
    do.call(rbind, rows)
  } else {
    data.frame(
      bpc = integer(0), group = group_labels[0], n_trials = integer(0),
      weight = numeric(0), p_value = numeric(0), stringsAsFactors = FALSE
    )
  }
  list(weights = weights, alpha = alpha)
}
