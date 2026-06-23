#' @title Cluster electrodes by their canonical \verb{CRP} response shape
#' @description
#' Groups recording electrodes by the \emph{shape} of their canonical evoked
#' response. Each electrode is summarized by \code{\link{crp}} into one
#' amplitude-normalized canonical shape; an electrode-by-electrode similarity
#' matrix is factorized with \code{\link{naive_nmf}} to find clusters that share
#' a response shape, and a representative basis profile curve is extracted per
#' cluster. This applies the basis-profile-curve (\verb{BPC}) approach
#' \emph{across electrodes} rather than across stimulation sites; see
#' \sQuote{References}.
#'
#' @details
#' \enumerate{
#' \item \strong{Common window.} Electrodes may have different time axes
#'   (\code{\link{crp}} drops \code{NA} samples); the overlapping time domain is
#'   used, \code{time_window} is clipped into it, and each electrode is subset on
#'   the fly with no interpolation.
#' \item \strong{Similarity.} With \code{paired = TRUE} each entry is the
#'   one-sample t-statistic, across an electrode pair's common trials, of the
#'   per-trial \verb{CCEP} cross-projection of the \eqn{1/\alpha}-rescaled
#'   responses; negatives are zeroed, the matrix made symmetric and scaled to a
#'   maximum of one. With \code{paired = FALSE} it is the cosine cross-projection
#'   of the \code{C_full} curves.
#' \item \strong{Rank selection.} \code{\link{naive_nmf}} factorizes the
#'   similarity at rank \eqn{Q} (from \code{initial_rank}); each rank is re-run a
#'   few times and the lowest-error fit kept. \eqn{Q} is reduced while
#'   \eqn{\zeta} - the sum of the upper off-diagonal of the row-normalized
#'   \eqn{HH^\top} - exceeds \code{zeta_threshold}.
#' \item \strong{Assignment.} Each electrode takes its winner-take-all cluster
#'   over the normalized \verb{NMF} loadings; with \code{null_class}, loadings
#'   below \eqn{1/(2\sqrt{N})} are left unassigned.
#' \item \strong{Basis curves.} Per cluster, the first linear kernel \verb{PCA}
#'   component of its members' \code{C_full} curves (as in \code{\link{crp}}).
#' }
#'
#' @param crp_list a named list of \code{ravetools_crp} objects (one per
#' electrode) from \code{\link{crp}}; list names are used as electrode labels.
#' Time axes need not match - the common overlapping domain is used (an error is
#' raised if the domains do not overlap). See \sQuote{Details}.
#' @param paired logical; \code{TRUE} (the default) builds the trial-level
#' \verb{BPC} t-statistic similarity and requires the single trials in
#' \code{.data$V}; \code{FALSE} uses the cosine similarity of the \code{C_full}
#' curves. See \sQuote{Details}.
#' @param time_window numeric \code{c(lo, hi)} analysis window, clipped into the
#' common domain; an \code{NA} extends to that domain edge. Defaults to
#' \code{c(0, NA)} (stimulus onset to the end of the domain).
#' @param n_clusters integer or \code{NULL}; if given, fixes the number of
#' clusters and skips the automatic selection (which otherwise uses
#' \code{zeta_threshold}).
#' @param initial_rank integer or \code{NULL}; the starting \verb{NMF} rank for
#' automatic selection. \code{NULL} uses \code{min(length(crp_list) - 1L, 10L)}.
#' @param zeta_threshold numeric \eqn{> 0}; rank-selection cutoff on \eqn{\zeta}
#' (see \sQuote{Details}); smaller values give fewer clusters. Defaults to
#' \code{1}, matching the \verb{BPC} reference.
#' @param null_class logical; if \code{TRUE} (the default) electrodes whose
#' winner-take-all loading falls below \eqn{1/(2\sqrt{N})} are left unassigned
#' (\code{NA}); if \code{FALSE} every electrode is forced into a cluster. See
#' \sQuote{Details}.
#' @param nmf_max_iters,nmf_tol passed to \code{\link{naive_nmf}}.
#' @param verbose logical; whether to report progress.
#'
#' @returns A named list of class \code{ravetools_crp_cluster}:
#' \describe{
#' \item{\code{clusters}}{Integer vector, the cluster index assigned to each
#'   electrode (winner-take-all over the row-normalized \verb{NMF} loadings).
#'   When \code{null_class = TRUE}, electrodes whose top loading falls below
#'   \eqn{1/(2\sqrt{N})} are left unassigned (\code{NA}).}
#' \item{\code{basis_curves}}{Numeric matrix, time \eqn{\times} number of
#'   clusters; column \eqn{q} is the basis profile curve \eqn{B_q(t)} for
#'   cluster \eqn{q}, the first linear-kernel \verb{PCA} component of its member
#'   \code{C_full} curves, sign-oriented to the cluster mean.}
#' \item{\code{basis_times}}{Numeric vector, the time axis for
#'   \code{basis_curves} (the common overlapping time axis, restricted to
#'   \code{time_window}).}
#' \item{\code{similarity}}{Numeric matrix, the electrode-by-electrode
#'   similarity \eqn{\Xi} that was factorized (non-negative, scaled to a maximum
#'   of one).}
#' \item{\code{nmf}}{The \code{\link{naive_nmf}} result for the selected rank.}
#' \item{\code{n_clusters}}{Integer, the number of clusters found.}
#' \item{\code{domain}}{Numeric \code{c(lo, hi)}, the overlapping time domain
#'   shared by all electrodes.}
#' \item{\code{paired},\code{time_window}}{The settings used; \code{time_window}
#'   is the effective window after clipping into \code{domain}.}
#' }
#'
#' @references
#' The \verb{BPC} method is described in \doi{10.1371/journal.pcbi.1008710}; the
#' underlying \verb{CRP} method in \doi{10.1371/journal.pcbi.1011105}.
#'
#' @seealso \code{\link{crp}}, \code{\link{naive_nmf}}
#'
#' @examples
#' # Four response shapes; shapes 3 and 4 differ only in amplitude, so they
#' # cluster together once shapes are amplitude-normalized.
#' \donttest{
#' set.seed(1)
#' n_time <- 300L
#' tt <- seq(-0.2, 1, length.out = n_time)
#' shapes <- list(
#'   exp(-((tt - 0.08) / 0.03)^2) - 0.5 * exp(-((tt - 0.18) / 0.04)^2),
#'   exp(-((tt - 0.38) / 0.03)^2) - 0.5 * exp(-((tt - 0.50) / 0.04)^2),
#'   exp(-((tt - 0.70) / 0.04)^2),
#'   exp(-((tt - 0.70) / 0.04)^2) * 2
#' )
#'
#' # 4 electrodes per shape, each parameterized with crp()
#' crp_list <- list()
#' for (g in seq_along(shapes)) {
#'   for (e in seq_len(4L)) {
#'     V <- outer(shapes[[g]], runif(15L, 0.5, 1.5)) +
#'       matrix(rnorm(n_time * 15L, sd = 0.15), n_time, 15L)
#'     crp_list[[sprintf("elec_%d_%d", g, e)]] <-
#'       crp(V, tt, remove_artifacts = FALSE)
#'   }
#' }
#'
#' res <- crp_cluster(crp_list, verbose = TRUE)
#' res$n_clusters
#' table(res$clusters)
#' plot(res)
#' }
#'
#' @export
crp_cluster <- function(
  crp_list,
  paired = TRUE,
  time_window = c(0, NA),
  n_clusters = NULL,
  initial_rank = NULL,
  zeta_threshold = 1,
  null_class = TRUE,
  nmf_max_iters = 1e4,
  nmf_tol = c(1e-4, 1e-8),
  verbose = TRUE
) {

  # DIPSAUS DEBUG START
  # pipeline <- ravepipeline::pipeline(
  #   pipeline_name = "voltage_explorer",
  #   paths = file.path(
  #     rstudioapi::getActiveProject(),
  #     "../rave-pipelines/modules/"
  #   ),
  #   temporary = TRUE
  # )
  # crp_results <- pipeline$run("crp_results")
  # crp_list <- unname(lapply(crp_results, function(x) { x[[1]] }))
  # paired = TRUE
  # time_window = c(0, 0.5)
  # normalize = TRUE
  # n_clusters = NULL
  # initial_rank = NULL
  # zeta_threshold = 0.5
  # nmf_max_iters = 1e4
  # nmf_tol = c(1e-4, 1e-8)
  # verbose = TRUE


  paired <- isTRUE(as.logical(paired))

  # ---- validate inputs ----------------------------------------------------
  if (!is.list(crp_list) || length(crp_list) < 2L) {
    stop("`crp_list` must be a list of at least 2 `ravetools_crp` objects.")
  }
  is_crp <- vapply(crp_list, inherits, logical(1L), what = "ravetools_crp")
  if (!all(is_crp)) {
    stop("Every element of `crp_list` must be a `ravetools_crp` object ",
         "returned by `crp()`.")
  }
  labels <- names(crp_list)
  if (is.null(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) {
    labels <- sprintf("electrode_%d", seq_along(crp_list))
  }
  names(crp_list) <- labels
  N <- length(crp_list)

  # ---- common (overlapping) time domain across electrodes -----------------
  # crp() drops non-finite (NA) rows, so electrodes may carry different time
  # axes. We only require their domains to overlap; subsetting is done on the
  # fly (by time value, per electrode) rather than assuming shared row indices.
  domain <- crp_cluster_domain(crp_list)

  # ---- comparison window: validate, then clip into the common domain ------
  time_window <- as.numeric(time_window)
  if (!length(time_window)) {
    window <- domain
  } else {
    if (length(time_window) == 1) {
      time_window <- c(0, time_window)
    } else {
      time_window <- time_window[c(1, 2)]
    }
    time_window[is.na(time_window)] <- domain[is.na(time_window)]
    time_window[time_window < domain[[1]]] <- domain[[1]]
    time_window[time_window > domain[[2]]] <- domain[[2]]
    if (is.unsorted(time_window)) {
      stop(
        sprintf(
          "`time_window` must be a numeric c(lo, hi) with lo < hi, and must overlap with the common time domain [%.4g, %.4g].",
          domain[1L], domain[2L])
      )
    }
    window <- time_window
  }

  # ---- electrode-by-electrode similarity ----------------------------------
  # paired = TRUE: trial-level cross-projection t-statistic over pairwise common
  # trials; paired = FALSE: cosine of the C_full curves. Returned directly.
  similarity <- crp_cluster_similarity_matrix(crp_list, window, paired)

  # ---- NMF clustering with optional automatic rank selection ---------------
  if (is.null(initial_rank)) {
    initial_rank <- min(N - 1L, 10L)
  }
  initial_rank <- max(1L, min(as.integer(initial_rank), N - 1L))
  if (!is.null(n_clusters)) {
    n_clusters <- max(1L, min(as.integer(n_clusters), N - 1L))
  }

  sel <- crp_cluster_nmf_select(
    similarity,
    n_clusters = n_clusters,
    initial_rank = initial_rank,
    zeta_threshold = zeta_threshold,
    nmf_max_iters = nmf_max_iters,
    nmf_tol = nmf_tol,
    verbose = verbose
  )
  nmf <- sel$nmf

  # ---- winner-take-all assignment, renumbered to present clusters ----------
  # Each electrode joins its top cluster on the row-normalized NMF loadings. When
  # `null_class = TRUE` (the BPC reference behavior) a membership is kept only if
  # that normalized loading exceeds 1/(2*sqrt(N)) - half the value a uniform
  # loading 1/sqrt(N) would give - else the electrode is left unassigned (NA, the
  # "null" class). When FALSE every electrode is forced into its winner cluster.
  Hn <- nmf$H
  rn <- sqrt(rowSums(Hn * Hn))
  rn[rn == 0] <- 1
  Hn <- Hn / rn
  raw_cluster <- apply(Hn, 2L, which.max)
  if (isTRUE(null_class)) {
    winner <- Hn[cbind(raw_cluster, seq_len(N))]
    raw_cluster[winner <= 1 / (2 * sqrt(N))] <- NA_integer_
  }
  present <- sort(unique(raw_cluster))
  clusters <- match(raw_cluster, present)
  n_found <- length(present)

  # ---- per-cluster basis curves + per-electrode metrics --------------------
  # Subsets each electrode's C_full on the fly to `window` using its own time.
  fm <- crp_cluster_feature_matrix(crp_list, window, TRUE)
  feature_mat <- fm$features
  basis_times <- fm$times
  basis <- crp_cluster_basis(feature_mat, clusters, n_found)
  basis_curves <- basis$basis_curves
  # metrics <- crp_cluster_metrics(feature_mat, clusters, basis_curves)

  rownames(basis_curves) <- NULL

  structure(
    class = "ravetools_crp_cluster",
    list(
      clusters = clusters,
      basis_curves = basis_curves,
      basis_times = basis_times,
      similarity = similarity,
      # metrics = metrics,
      nmf = nmf,
      n_clusters = n_found,
      paired = paired,
      domain = domain,
      time_window = window
    )
  )
}

#' @title Plot electrode clustering results
#' @description S3 plot method for \code{ravetools_crp_cluster} objects returned
#' by \code{\link{crp_cluster}}. Draws a two-panel figure: the per-cluster basis
#' profile curves, and the electrode-by-electrode similarity matrix ordered by
#' cluster.
#' @param x a \code{ravetools_crp_cluster} object.
#' @param ... ignored.
#' @return Invisibly returns \code{x}.
#' @export
plot.ravetools_crp_cluster <- function(x, ...) {
  Q <- x$n_clusters
  cols <- grDevices::hcl.colors(max(Q, 2L), palette = "Dark 3")

  op <- graphics::par(mfrow = c(1L, 2L))
  on.exit(graphics::par(op), add = TRUE)

  # ---- Panel 1: basis profile curves --------------------------------------
  bt <- x$basis_times
  bc <- x$basis_curves
  graphics::matplot(
    bt, bc, type = "l", lty = 1L, lwd = 2L, col = cols[seq_len(Q)],
    las = 1L, xlab = "Time (s)", ylab = "Basis amplitude (a.u.)",
    main = "Basis profile curves"
  )
  graphics::abline(h = 0, col = "#cccccc", lty = 3L)
  graphics::legend(
    "topright", bty = "n", cex = 0.8, lwd = 2L, col = cols[seq_len(Q)],
    legend = sprintf("cluster %d (n=%d)", seq_len(Q),
                     tabulate(x$clusters, nbins = Q))
  )

  # ---- Panel 2: similarity matrix ordered by cluster ----------------------
  ord <- order(x$clusters)
  S <- x$similarity[ord, ord, drop = FALSE]
  graphics::image(
    seq_len(nrow(S)), seq_len(ncol(S)), S,
    col = grDevices::hcl.colors(64L, palette = "YlOrRd", rev = TRUE),
    xlab = "Electrode (ordered by cluster)", ylab = "Electrode",
    main = "Similarity matrix", las = 1L
  )
  # cluster block boundaries
  brk <- cumsum(tabulate(x$clusters, nbins = Q))
  brk <- brk[-length(brk)] + 0.5
  graphics::abline(v = brk, h = brk, col = "#202020", lty = 1L, lwd = 1L)

  invisible(x)
}


# --------------------------------------------------------------------------
# Internal helpers (not exported)
# --------------------------------------------------------------------------

# Overlapping time domain (range) shared by all electrodes. Different crp objects
# may carry different time axes (crp() drops non-finite rows), so the common
# domain is [max of the per-electrode minima, min of the per-electrode maxima].
crp_cluster_domain <- function(crp_list) {
  rng <- vapply(crp_list, function(obj) {
    range(obj$parameters$params_times_full)
  }, numeric(2L))
  lo <- max(rng[1L, ])
  hi <- min(rng[2L, ])
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
    stop("The `crp` objects have no overlapping time domain.")
  }
  c(lo, hi)
}

# Assemble the time x electrode feature matrix from each electrode's C_full,
# subsetting on the fly to the (domain-clipped) `window` using each electrode's
# own time axis. Because the axes may be heterogeneous, the selected timepoints
# are checked for agreement across electrodes (no interpolation). Returns the
# feature matrix and its shared time axis. Optionally unit-normalized so
# electrodes cluster by waveform shape rather than amplitude.
crp_cluster_feature_matrix <- function(crp_list, window, normalize) {
  picked <- lapply(crp_list, function(obj) {
    tt <- obj$parameters$params_times_full
    sel <- tt >= window[1L] & tt <= window[2L]
    list(times = tt[sel], values = obj$parameters$C_full[sel])
  })

  ref_times <- picked[[1L]]$times
  if (length(ref_times) < 2L) {
    stop("Fewer than 2 timepoints fall within the common time domain.")
  }
  for (i in seq_along(picked)) {
    if (length(picked[[i]]$times) != length(ref_times) ||
        !isTRUE(all.equal(picked[[i]]$times, ref_times, tolerance = 1e-8))) {
      stop("Electrodes have misaligned timepoints within the common domain; ",
           "`C_full` cannot be stacked without interpolation (electrode '",
           names(crp_list)[i], "' differs).")
    }
  }

  feature_mat <- vapply(picked, function(p) {
    p$values
  }, numeric(length(ref_times)))
  feature_mat <- as.matrix(feature_mat)
  if (isTRUE(normalize)) {
    norms <- sqrt(colSums(feature_mat * feature_mat))
    norms[norms == 0] <- 1
    feature_mat <- t(t(feature_mat) / norms)
  }
  feature_mat[!is.finite(feature_mat)] <- 0
  list(features = feature_mat, times = ref_times)
}

# Electrode-by-electrode similarity from the cross-projection of normalized
# C_full curves (BPC's Xi). Negatives are zeroed and the matrix is scaled to a
# maximum of one, matching the BPC significance-matrix convention.
crp_cluster_cross_matrix <- function(feature_mat) {
  norms <- sqrt(colSums(feature_mat * feature_mat))
  norms[norms == 0] <- 1
  F0 <- t(t(feature_mat) / norms)
  P <- crossprod(F0)            # cosine similarity, electrodes x electrodes
  P <- (P + t(P)) / 2           # enforce symmetry
  P[P < 0] <- 0
  mx <- max(P)
  if (mx > 0) P <- P / mx
  P
}

# Electrode-by-electrode similarity matrix (BPC significance matrix), returned
# directly. When `paired` is TRUE it is the one-sample t-statistic, across each
# electrode pair's *common* trials, of the per-trial CCEP cross-projection
# between the alpha-rescaled single-trial responses. When FALSE it falls back to
# the cosine cross-projection of the C_full curves (no trials needed). Both are
# made non-negative, symmetric and scaled to a maximum of one for `naive_nmf`.
crp_cluster_similarity_matrix <- function(crp_list, window, paired) {
  if (!isTRUE(paired)) {
    fm <- crp_cluster_feature_matrix(crp_list, window, normalize = TRUE)
    return(crp_cluster_cross_matrix(fm$features))
  }

  N <- length(crp_list)

  # alpha-rescaled, window-subset single trials per electrode (~ canonical shape
  # per trial, so amplitude variance is removed and the t-statistic sharpens)
  similarity <- lapply(seq_len(N), function(i) {
    obj_i <- crp_list[[i]]
    # al <- obj$parameters$al
    # al[!is.finite(al) | al == 0] <- 1
    tt_i <- obj_i$parameters$params_times
    # tt_i <- obj_i$.data$time
    sel_i <- tt_i >= window[1L] & tt_i <= window[2L]
    bad_trials_i <- obj_i$bad_trials

    tp_i <- sum(sel_i)
    window_i <- range(tt_i[sel_i])
    trial_i <- seq_len(ncol(obj_i$parameters$V_tR) + length(bad_trials_i))
    # trial_i <- seq_len(ncol(obj_i$.data$V) + length(bad_trials_i))
    if (length(bad_trials_i)) {
      trial_i <- trial_i[-bad_trials_i]
    }


    X_i <- vapply(seq_len(N), function(j) {

      # if (i == j) { return(0.0) }

      obj_j <- crp_list[[j]]

      bad_trials_j <- obj_j$bad_trials
      trial_j <- seq_len(ncol(obj_j$parameters$V_tR) + length(bad_trials_j))
      # trial_j <- seq_len(ncol(obj_j$.data$V) + length(bad_trials_j))
      if (length(bad_trials_j)) {
        trial_j <- trial_j[-bad_trials_j]
      }

      common_trials <- intersect(trial_i, trial_j)

      if (length(common_trials) <= 1) { return(0.0) }

      # crp() requires increasing time, paired=TRUE already assumed time consistency
      tt_j <- obj_j$parameters$params_times
      # tt_j <- obj_j$.data$time
      sel_j <- tt_j >= window_i[1L] & tt_j <= window_i[2L]

      tp_j <- sum(sel_j)
      if (!tp_j) { return(0.0) }

      window_j <- range(tt_j[sel_j])
      sel_i <- tt_i >= window_j[1L] & tt_i <= window_j[2L]
      tp_i <- sum(sel_i)
      if (!tp_i) { return(0.0) }

      # Assuming sum(sel_i) == sum(sel_j)
      if (tp_i != tp_j) {
        stop(sprintf("Number of time-point mismatch between object %d (%d points) vs %d (%d points)",
                     i, tp_i, j, tp_j))
      }

      trial_sel_j <- trial_j %in% common_trials
      vi <- obj_i$parameters$V_tR[sel_i, trial_i %in% common_trials, drop = FALSE]
      vj <- obj_j$parameters$V_tR[sel_j, trial_sel_j, drop = FALSE]
      # vi <- obj_i$.data$V[sel_i, trial_i %in% common_trials, drop = FALSE]
      # vj <- obj_j$.data$V[sel_j, trial_sel_j, drop = FALSE]

      # Calculate rescales
      # for vi, normalize to norm_2 = 1
      norm_i <- sqrt(colSums(vi ^ 2))
      norm_i[norm_i == 0] <- 1
      ui <- t(vi) / norm_i

      # For vj, rescale by `al`
      al_j <- obj_j$parameters$al
      al_j[al_j == 0] <- 1
      uj <- t(vj) / al_j[trial_sel_j]

      S <- rowSums(ui * uj)
      S <- S[is.finite(S)]
      if (length(S) <= 1) { return(0.0) }

      sd_s <- stats::sd(S)
      if (!is.finite(sd_s) || sd_s == 0) { return(0.0) }

      # standardized S
      mean(S) / (sd_s / sqrt(length(S)))
    }, 0.0)
    X_i
  })

  similarity <- t(simplify2array(unname(similarity), higher = TRUE))
  similarity[similarity < 0 | !is.finite(similarity)] <- 0
  # diag(similarity) <- 0

  # off-diagonal structure drives clustering: scale those to [0, 1] and set a
  # unit self-similarity (the diagonal t-statistic of trial norms is huge and
  # would otherwise dominate the max-scaling).
  m <- max(similarity)
  if (m > 0) {
    similarity <- similarity / m
  }

  # diag(similarity) <- 1

  similarity
}

# Run NMF on the similarity matrix, optionally reducing the rank automatically.
# When `n_clusters` is NULL the rank starts at `initial_rank` and drops by one
# while zeta - the sum of the upper-half off-diagonal elements of (row-normalized)
# HH^T - exceeds `zeta_threshold` (the BPC degeneracy criterion; see
# `crp_cluster_zeta`).
crp_cluster_nmf_select <- function(
  similarity, n_clusters, initial_rank, zeta_threshold,
  nmf_max_iters, nmf_tol, verbose
) {
  # NNMF is sensitive to its random initialization, so - following the BPC
  # reference (which uses 20) - re-run it a few times per rank and keep the
  # lowest-reconstruction-error factorization. 5 is sufficient here.
  n_restarts <- 5L
  fit_rank <- function(k) {
    fits <- lapply(seq_len(n_restarts), function(.) {
      naive_nmf(similarity, k = k, tol = nmf_tol,
                max_iters = nmf_max_iters, verbose = FALSE)
    })
    err <- vapply(fits, function(f) f$error[[1L]], numeric(1L))
    fits[[which.min(err)]]
  }

  if (!is.null(n_clusters)) {
    return(list(nmf = fit_rank(n_clusters), zeta = NA_real_, rank = n_clusters))
  }

  k <- initial_rank
  nmf <- fit_rank(k)
  zeta <- crp_cluster_zeta(nmf$H)
  while (k > 1L && is.finite(zeta) && zeta > zeta_threshold) {
    if (isTRUE(verbose)) {
      message(sprintf("rank %d: HH^T off-diagonal sum = %.3f > %.2f, reducing",
                      k, zeta, zeta_threshold))
    }
    k <- k - 1L
    nmf <- fit_rank(k)
    zeta <- crp_cluster_zeta(nmf$H)
  }
  if (isTRUE(verbose)) {
    message(sprintf("selected rank %d (HH^T off-diagonal sum = %.3f)", k, zeta))
  }
  list(nmf = nmf, zeta = zeta, rank = k)
}

# NMF degeneracy, per the BPC reference (Miller et al.): the rows of H (the NNMF
# components' loadings over electrodes) are unit-normalized, then zeta is the SUM
# of the upper-half off-diagonal elements of HH^T - the total shared structure
# between components. (The reference notebook normalizes H rows and sums
# triu(HH^T, 1); the paper also mentions a max variant, but the released code
# uses the sum.) A degenerate factorization has high overlap; once the rank
# matches the true clusters the components become near-disjoint and the sum drops.
crp_cluster_zeta <- function(H) {
  if (nrow(H) < 2L) return(0)
  rn <- sqrt(rowSums(H * H))
  rn[rn == 0] <- 1
  Hn <- H / rn                            # unit-normalize each component (row)
  HHt <- tcrossprod(Hn)                   # normalized HH^T
  sum(HHt[upper.tri(HHt)])                # BPC: sum of upper off-diagonals
}

# Per-cluster basis profile curve: the first linear-kernel PCA component of the
# member C_full curves, sign-oriented to the cluster mean. Reuses `crp_kt_pca`.
crp_cluster_basis <- function(feature_mat, clusters, n_found) {
  Tn <- nrow(feature_mat)
  basis_curves <- matrix(0, nrow = Tn, ncol = n_found)
  for (q in seq_len(n_found)) {
    members <- which(clusters == q)
    Xq <- feature_mat[, members, drop = FALSE]
    pca <- crp_kt_pca(Xq)
    Bq <- pca$E[, 1L]
    if (sum(Bq * rowMeans(Xq)) < 0) Bq <- -Bq
    basis_curves[, q] <- Bq
  }
  list(basis_curves = basis_curves)
}

# Per-electrode parameterization against its assigned basis curve: alpha
# (projection), signal-to-noise and explained variance, defined exactly as in
# `crp_core`.
crp_cluster_metrics <- function(feature_mat, clusters, basis_curves) {
  N <- ncol(feature_mat)
  alpha <- numeric(N)
  Vsnr <- numeric(N)
  expl_var <- numeric(N)
  for (i in seq_len(N)) {
    f <- feature_mat[, i]
    B <- basis_curves[, clusters[i]]
    a <- sum(B * f)
    ep <- f - a * B
    epep <- sum(ep * ep)
    ff <- sum(f * f)
    alpha[i] <- a
    Vsnr[i] <- if (epep > 0) a / sqrt(epep) else NA_real_
    expl_var[i] <- if (ff > 0) 1 - epep / ff else NA_real_
  }
  data.frame(
    cluster = unname(clusters),
    alpha = alpha,
    Vsnr = Vsnr,
    expl_var = expl_var,
    stringsAsFactors = FALSE
  )
}
