#' @title Canonical Response Parameterization (\verb{CRP})
#' @description
#' Parameterizes single-trial evoked responses (e.g. cortico-cortical evoked
#' potentials, \verb{CCEPs}) using the Canonical Response Parameterization
#' method (see 'Citation'). The function estimates the response
#' duration \eqn{\tau_R}, the time after stimulus at which the evoked
#' response has its most consistent, shared structure across trials.
#' The estimator is obtained from the time course of cross-trial projection
#' magnitudes, extracts the canonical response shape \eqn{C(t)} via a linear
#' kernel-trick PCA on the trial matrix truncated at \eqn{\tau_R}, and
#' reports per-trial weights, residuals, signal-to-noise, explained variance
#' and extraction-significance statistics.
#'
#' This is an R translation of \code{CRP_method.m} (and the surrounding
#' artifact-rejection / duration-uncertainty logic in
#' \code{CRP_illustration.m}) from the upstream MATLAB reference
#' implementation; see \sQuote{References}.
#'
#' @param x numeric matrix of single-trial evoked voltages with shape
#' \code{time x trials} (i.e. rows are timepoints, columns are trials);
#' the matrix orientation matches the variable \code{V} / \code{data} in
#' the MATLAB reference. At least two trials are required.
#' @param time numeric vector of length \code{nrow(x)} giving the
#' stimulus-aligned time (in seconds) of each row of \code{x}; must be
#' monotonically increasing and span \code{[t_start, t_end]}.
#' @param t_start,t_end numeric scalars, post-stimulation start and end times
#' (in seconds) defining the analysis window. Defaults match the MATLAB
#' illustration (\code{0.015 s} to \code{1 s}).
#' @param remove_artifacts logical; if \code{TRUE} (the default), an initial
#' \verb{CRP} pass is run to identify outlier/artifact trials, which are then
#' dropped before the final pass. See \sQuote{Details}.
#' @param artifact_interval character, one of \code{"full"} (the default,
#' matching the active option in the MATLAB illustration) or \code{"tR"};
#' selects whether per-trial outlier statistics are computed on the
#' projection magnitudes for the full window or only at the response
#' duration \eqn{\tau_R}.
#' @param artifact_p_threshold numeric, p-value threshold below which a
#' trial is flagged as artifact (provided its mean projection is also
#' below the cohort mean); defaults to \code{1e-5}.
#' @param threshold_quantile numeric in \code{(0, 1)}; the fraction of the
#' peak mean projection magnitude used to derive the duration-uncertainty
#' bounds \code{tau_R_lower} and \code{tau_R_upper}. Defaults to
#' \code{0.98} as in the manuscript.
#' @param time_step integer, sampling step (in samples) used when sweeping
#' candidate response duration; defaults to \code{5L}, matching the
#' MATLAB \code{t_step}. Larger values are faster but smooth the
#' projection profile.
#'
#' @returns A named list with the following elements:
#' \describe{
#' \item{\code{parameters}}{A list of single-trial parameterizations
#' (\code{crp_parms} in MATLAB):
#' \describe{
#'   \item{\code{C}}{Numeric vector of length \eqn{T_R} (timepoints up to
#'     \eqn{\tau_R}), the canonical response shape \eqn{C(t)}: the first
#'     eigenvector of the linear kernel PCA on \code{V_tR}, unit-normalized
#'     (\eqn{\|C\| = 1}). The matching time axis is in \code{params_times}.}
#'   \item{\code{al}}{Numeric vector of length \eqn{K} (number of trials),
#'     the per-trial alpha coefficient \eqn{\alpha_k = C^\top V_k}: scalar
#'     projection of trial \eqn{k} onto \eqn{C(t)}. Larger magnitude means
#'     the trial resembles the canonical shape more strongly; sign reflects
#'     polarity relative to \eqn{C}.}
#'   \item{\code{al_p}}{Numeric vector of length \eqn{K}, alpha-prime
#'     \eqn{\alpha_k / \sqrt{T_R}}: \code{al} rescaled to remove the
#'     duration dependence from the unit-norm convention on \eqn{C}.
#'     Expressed in \eqn{\mu V} and comparable across electrodes or
#'     conditions with different \eqn{\tau_R}.}
#'   \item{\code{ep}}{Numeric matrix of shape \eqn{T_R \times K}, the
#'     per-trial residual \eqn{\epsilon_k(t) = V_k(t) - \alpha_k C(t)}
#'     after the shared component is removed. Access trial \eqn{k} via
#'     \code{ep[, k]}.}
#'   \item{\code{epep_root}}{Numeric vector of length \eqn{K},
#'     \eqn{\|\epsilon_k\| = \sqrt{\epsilon_k^\top \epsilon_k}}: L2 norm
#'     of the residual per trial. Smaller values indicate the canonical
#'     shape describes that trial more faithfully.}
#'   \item{\code{Vsnr}}{Numeric vector of length \eqn{K}, per-trial
#'     signal-to-noise \eqn{\alpha_k / \|\epsilon_k\|}. Values \eqn{> 1}
#'     indicate the canonical component is larger than the residual.}
#'   \item{\code{expl_var}}{Numeric vector of length \eqn{K}, per-trial
#'     explained variance \eqn{1 - \|\epsilon_k\|^2 / \|V_k\|^2}: fraction
#'     of each trial's energy accounted for by \eqn{\alpha_k C(t)}.
#'     Ranges in \eqn{[0, 1]}.}
#'   \item{\code{tR}}{Numeric scalar, response duration \eqn{\tau_R} in
#'     seconds: the time at which mean cross-trial projection magnitude is
#'     maximized.}
#'   \item{\code{params_times}}{Numeric vector of length \eqn{T_R}, time
#'     axis for \code{C}, \code{V_tR}, \code{ep}, and
#'     \code{avg_trace_tR}.}
#'   \item{\code{V_tR}}{Numeric matrix \eqn{T_R \times K}, trial matrix
#'     truncated to \eqn{\tau_R} - the data actually decomposed.}
#'   \item{\code{avg_trace_tR}}{Numeric vector of length \eqn{T_R}, simple
#'     trial average truncated to \eqn{\tau_R}.}
#' }}
#' \item{\code{projections}}{A list of projection-stage outputs
#' (\code{crp_projs} in MATLAB):
#' \describe{
#'   \item{\code{proj_tpts}}{Numeric vector, candidate duration time points
#'     (seconds) at which projection magnitudes were evaluated.}
#'   \item{\code{S_all}}{Numeric matrix; rows are non-redundant off-diagonal
#'     trial-pair projections, columns correspond to \code{proj_tpts}.
#'     Units: \eqn{\mu V \cdot s^{1/2}}.}
#'   \item{\code{mean_proj_profile}}{Numeric vector, mean of \code{S_all}
#'     across trial pairs at each candidate duration, the profile whose
#'     maximum defines \eqn{\tau_R}.}
#'   \item{\code{var_proj_profile}}{Numeric vector, variance of
#'     \code{S_all} across trial pairs at each candidate duration.}
#'   \item{\code{tR_index}}{Integer, column index into \code{S_all} and
#'     \code{proj_tpts} corresponding to \eqn{\tau_R}.}
#'   \item{\code{avg_trace_input}}{Numeric vector, simple trial average
#'     over the full analysis window (not truncated to \eqn{\tau_R}).}
#'   \item{\code{stat_indices}}{Integer vector, row indices of \code{S_all}
#'     used for the significance t-tests, constructed so each trial-pair
#'     comparison appears at most once.}
#'   \item{\code{t_value_tR}, \code{p_value_tR}}{t-statistic and one-sided
#'     p-value (H1: mean projection \eqn{> 0}) at \eqn{\tau_R}. Primary
#'     extraction-significance test reported in the manuscript.}
#'   \item{\code{t_value_full}, \code{p_value_full}}{Same test at the full
#'     analysis-window duration.}
#' }}
#' \item{\code{bad_trials}}{Integer vector of column indices into the
#' original \code{x} flagged and removed as artifacts; \code{integer(0)}
#' when none removed or when \code{remove_artifacts = FALSE}.}
#' \item{\code{tau_R}}{Numeric scalar, estimated response duration
#' \eqn{\tau_R} in seconds (convenience copy of \code{parameters$tR}).}
#' \item{\code{tau_R_lower}, \code{tau_R_upper}}{Numeric scalars, lower
#' and upper threshold-crossing times (seconds) bracketing \eqn{\tau_R}
#' at the \code{threshold_quantile} fraction of the peak mean projection
#' magnitude.}
#' \item{\code{t_start}, \code{t_end}}{The analysis window used.}
#' \item{\code{sample_rate}}{Numeric, sampling rate inferred from
#' \code{time}.}
#' }
#'
#' @details
#' Briefly, the algorithm proceeds in three stages:
#' \enumerate{
#' \item For a sweep of candidate durations \eqn{k}, compute pairwise
#' L2-normalized cross-projection magnitudes between trials truncated to
#' \eqn{[0, k]}. The duration that maximizes the mean projection magnitude
#' is taken as the response duration \eqn{\tau_R}.
#' \item Apply linear kernel-trick PCA to the trial matrix truncated to
#' \eqn{\tau_R}; the first principal component is the canonical response
#' shape \eqn{C(t)}.
#' \item Project \eqn{C(t)} into each trial to obtain per-trial weights
#' \eqn{\alpha_k}; the residual \eqn{\epsilon_k = V_k - \alpha_k C}
#' summarizes trial-by-trial deviation from the canonical shape.
#' }
#' Significance is assessed by a one-sided t-test on the off-diagonal
#' projection magnitudes against zero, restricted to a non-overlapping
#' subset of comparison pairs to avoid double-counting.
#'
#' When \code{remove_artifacts = TRUE}, the function performs an initial
#' \verb{CRP} pass and runs an unpaired t-test for each trial comparing the
#' projections it participates in against all other off-diagonal
#' projections. Trials with \code{p < artifact_p_threshold} \emph{and}
#' mean projection below the cohort mean are dropped, and \verb{CRP} is re-run.
#'
#' @references
#' The \verb{CRP} algorithm is described in \doi{10.1371/journal.pcbi.1011105},
#' with a reference MATLAB implementation at
#' \url{https://github.com/kaijmiller/crp_scripts}. See
#' \code{citation("ravetools")} for the full bibliographic entry.
#'
#' @examples
#' set.seed(42)
#'
#' # Synthetic CCEP-like data: shared canonical shape with per-trial scaling
#' n_time   <- 500L
#' n_trials <- 20L
#' tt <- seq(-0.05, 1, length.out = n_time)
#' canonical <- exp(-((tt - 0.10) / 0.05)^2) -
#'              0.5 * exp(-((tt - 0.30) / 0.10)^2)
#' V <- (outer(canonical, runif(n_trials, 0.5, 1.5)) +
#'        matrix(rnorm(n_time * n_trials, sd = 0.3), n_time, n_trials)) * 2
#'
#' res <- crp(V, tt)
#'
#' op <- par(mfrow = c(1, 3), mar = c(4.5, 4, 3, 1))
#' on.exit({ par(op) })
#'
#' # ---- Panel 1: all trials (full window) + mean + C(t) overlay ----------
#' parms <- res$parameters
#' matplot(tt, V, type = "l", lty = 1,
#'         col = "#80808060", xlab = "Time (s)",
#'         ylab = expression(mu * V),
#'         main = expression("Canonical shape " * C(t)))
#' # scale C(t) to the amplitude of the mean trace for overlay;
#' # C(t) ends at tau_R so the line is cut off there naturally
#' C_scaled <- parms$C * max(abs(rowMeans(V))) / max(abs(parms$C))
#' lines(parms$params_times, C_scaled, col = "#FFFF0080", lwd = 3)
#'
#' # Mean
#' lines(tt, rowMeans(V), col = "black", lwd = 1)
#' legend("topright", c("mean", "C(t) scaled"),
#'        col = c("black", "#FFFF00"), lty = c(1, 2), lwd = 2,
#'        bty = "n", cex = 0.8)
#'
#' # ---- Panel 2: per-trial alpha-prime weights -----------------------------
#' barplot(sort(parms$al_p), col = "steelblue", border = NA, las = 1,
#'         xlab = "Trial (sorted)",
#'         ylab = expression(alpha * "'" ~ (mu * V)),
#'         main = expression("Per-trial " * alpha * "' (alpha-prime)"))
#' abline(h = c(0, mean(parms$al_p)), lty = 2)
#'
#' # ---- Panel 3: mean projection profile with tau_R bounds ----------------
#' proj <- res$projections
#' plot(proj$proj_tpts, proj$mean_proj_profile, type = "l", lwd = 2,
#'      xlab = "Candidate duration (s)",
#'      ylab = expression(bar(S) ~ (mu * V %.% s^{0.5})),
#'      main = expression("Projection profile & " * tau[R]),
#'      las = 1)
#' abline(v = c(res$tau_R_lower, res$tau_R, res$tau_R_upper),
#'        col = c("cyan3", "orange2", "red"),
#'        lty = c(2, 1, 2), lwd = 2)
#' legend("topright",
#'        legend = expression(tau[lb], tau[R], tau[ub]),
#'        col = c("cyan3", "orange2", "red"),
#'        lty = c(2, 1, 2), lwd = 2, bty = "n")
#'
#' par(op)
#'
#' @export
crp <- function(
  x, time,
  t_start = 0.015, t_end = 1,
  remove_artifacts = TRUE,
  artifact_interval = c("full", "tR"),
  artifact_p_threshold = 1e-5,
  threshold_quantile = 0.98,
  time_step = 5L
) {

  # ---- input validation ----------------------------------------------------
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric matrix with shape time x trials.")
  }
  if (ncol(x) < 2L) {
    stop("`x` must have at least 2 trials (columns).")
  }
  if (!is.numeric(time) || length(time) != nrow(x)) {
    stop("`time` must be a numeric vector with length equal to nrow(x).")
  }
  if (any(diff(time) <= 0)) {
    stop("`time` must be strictly monotonically increasing.")
  }
  if (time[1L] > t_start) {
    stop("The first element of `time` is greater than `t_start`.")
  }
  if (time[length(time)] < t_end) {
    stop("The last element of `time` is less than `t_end`.")
  }
  artifact_interval <- match.arg(artifact_interval)
  if (!is.numeric(threshold_quantile) || length(threshold_quantile) != 1L ||
      threshold_quantile <= 0 || threshold_quantile >= 1) {
    stop("`threshold_quantile` must be a single number in (0, 1).")
  }
  time_step <- as.integer(time_step)
  if (length(time_step) != 1L || is.na(time_step) || time_step < 1L) {
    stop("`time_step` must be a positive integer scalar.")
  }

  # ---- subset to analysis window ------------------------------------------
  tpts <- which(time > t_start & time <= t_end)
  if (length(tpts) < 11L) {
    stop("Analysis window contains too few samples (need >= 11).")
  }
  V <- x[tpts, , drop = FALSE]
  t_win <- time[tpts]
  srate <- 1 / mean(diff(time))

  # ---- optional artifact removal ------------------------------------------
  bad_trials <- integer(0)
  if (isTRUE(remove_artifacts)) {
    init <- .crp_core(V, t_win, srate, time_step)
    mp <- .crp_trial_tests(init$crp_projs, artifact_interval, ncol(V))
    flag <- mp$p < artifact_p_threshold & mp$m < mean(mp$m)
    bad_trials <- which(flag)
    if (length(bad_trials) > 0L) {
      if (ncol(V) - length(bad_trials) < 2L) {
        stop(
          "Artifact removal would leave fewer than 2 trials; ",
          "consider raising `artifact_p_threshold` or setting ",
          "`remove_artifacts = FALSE`."
        )
      }
      V <- V[, -bad_trials, drop = FALSE]
    }
  }

  # ---- final CRP pass ------------------------------------------------------
  out <- .crp_core(V, t_win, srate, time_step)
  crp_projs <- out$crp_projs
  crp_params <- out$crp_params

  # ---- duration-uncertainty bounds (threshold_quantile of peak) ----------
  m <- crp_projs$mean_proj_profile
  mp <- max(m)
  cidx <- which.max(m)
  thr <- threshold_quantile * mp

  n <- length(m)
  up_cross <- 1L + which(m[-1L] > thr & m[-n] <= thr)
  up_cross <- up_cross[up_cross <= cidx]
  lb <- if (length(up_cross) > 0L) up_cross[length(up_cross)] else 1L

  down_cross <- which(m[-n] > thr & m[-1L] <= thr)
  down_cross <- down_cross[down_cross >= cidx]
  hb <- if (length(down_cross) > 0L) down_cross[1L] else n

  tau_R_lower <- crp_projs$proj_tpts[lb]
  tau_R_upper <- crp_projs$proj_tpts[hb]

  structure(
    class = "ravetools_crp",
    list(
      parameters = crp_params,
      projections = crp_projs,
      bad_trials = bad_trials,
      tau_R_lower = tau_R_lower,
      tau_R = crp_params$tR,
      tau_R_upper = tau_R_upper,
      t_start = t_start,
      t_end = t_end,
      sample_rate = srate,
      .data = list(V = V, time = t_win)
    )
  )
}

#' @title Plot \verb{CRP} results
#' @description
#' S3 plot method for objects of class \code{ravetools_crp} returned by
#' \code{\link{crp}}.  Produces a three-panel figure:
#' \enumerate{
#' \item Single-trial traces over the full analysis window with the mean
#'   and the scaled canonical shape \eqn{C(t)} overlaid (the shape is
#'   drawn only up to \eqn{\tau_R}, so the cut-off is itself informative).
#' \item Per-trial \eqn{\alpha'} weights sorted in ascending order.
#' \item Mean cross-trial projection profile with vertical lines marking
#'   \eqn{\tau_{lb}}, \eqn{\tau_R} and \eqn{\tau_{ub}}.
#' }
#' @param x an object of class \code{ravetools_crp} as returned by
#' \code{\link{crp}}.
#' @param ... additional graphical parameters passed to \code{\link{par}}
#' (e.g. \code{mar}, \code{cex.axis}); currently unused beyond restoring
#' the previous \code{par} state on exit.
#' @return Invisibly returns \code{x}.
#' @export
plot.ravetools_crp <- function(x, ...) {
  parms <- x$parameters
  proj  <- x$projections
  V  <- x$.data$V
  tt <- x$.data$time

  op <- graphics::par(mfrow = c(1, 3))
  on.exit(graphics::par(op), add = TRUE)

  # ---- Panel 1: all trials (full window) + mean + C(t) overlay ----------
  graphics::matplot(tt, V, type = "l", lty = 1, las = 1,
                    col = "#80808060", xlab = "Time (s)",
                    ylab = expression(mu * V),
                    main = expression("Canonical shape " * C(t)))

  # C(t) is drawn only up to tau_R; scale amplitude to match mean trace
  C_scaled <- parms$C * max(abs(rowMeans(V))) / max(abs(parms$C))

  graphics::abline(v = x$tau_R, col = "red", lwd = 1)
  graphics::text(
    x = x$tau_R, y = min(C_scaled),
    labels = bquote(tau[R] * "=" * .(sprintf("%.3s", x$tau_R))),
    col = "red", pos = 4, offset = 1
  )

  graphics::lines(parms$params_times, C_scaled, col = "#FFFF0080", lwd = 3)
  graphics::lines(tt, rowMeans(V), col = "black", lwd = 1)
  graphics::legend("topright", c("mean", "C(t) scaled"),
                   col = c("black", "#FFFF00"), lty = c(1, 2), lwd = 3,
                   bty = "n", cex = 0.8)

  # ---- Panel 2: per-trial alpha-prime weights -----------------------------
  graphics::barplot(
    sort(parms$al_p),
    col = "steelblue",
    border = NA,
    las = 1,
    xlab = "Trial (sorted)",
    ylab = expression(alpha * "'" ~ (mu * V)),
    main = expression("Per-trial " * alpha * "' (alpha-prime)")
  )
  graphics::abline(h = c(0, mean(parms$al_p)), lty = c(1, 2))

  # ---- Panel 3: mean projection profile with tau_R bounds ----------------
  graphics::plot(
    proj$proj_tpts, proj$mean_proj_profile, type = "l", lwd = 2,
    xlab = "Candidate duration (s)",
    ylab = expression(bar(S) ~ (mu * V %.% s^{0.5})),
    main = expression("Projection profile & " * tau[R]),
    las = 1
  )
  graphics::abline(
    v = c(x$tau_R_lower, x$tau_R, x$tau_R_upper),
    col = c("grey", "red", "grey"),
    lty = c(2, 1, 2), lwd = 2)
  graphics::legend(
    "topright",
    legend = c(
     bquote(tau[lb] * "=" * .(sprintf("%.3fs", x$tau_R_lower))),
     bquote(tau[R] * "=" * .(sprintf("%.3fs", x$tau_R))),
     bquote(tau[ub] * "=" * .(sprintf("%.3fs", x$tau_R)))
    ),
    col = c("grey", "red", "grey"),
     lty = c(2, 1, 2), lwd = 2, bty = "n")

  invisible(x)
}


# --------------------------------------------------------------------------
# Internal helpers (not exported)
# --------------------------------------------------------------------------

# Single-trial cross-projection magnitudes (MATLAB `ccep_proj`).
# V is time x trials; returns a numeric vector of off-diagonal projections.
.crp_ccep_proj <- function(V) {
  norms <- sqrt(colSums(V * V))
  # remove zero norms since 0/1=0
  norms[norms == 0] <- 1
  V0 <- t(t(V) / norms)
  V0[!is.finite(V0)] <- 0
  P <- crossprod(V0, V) # trials x trials
  diag(P) <- NA_real_
  S0 <- as.numeric(P)
  S0[!is.na(S0)]
}

# Linear kernel-trick PCA (MATLAB `kt_pca`).
# Returns columns of E as unit-normalized eigenvectors of X (descending),
# and S as corresponding singular values (sqrt of eigenvalues).
.crp_kt_pca <- function(X) {
  XtX <- crossprod(X)
  eg <- eigen(XtX, symmetric = TRUE)
  # eigen() already returns values in descending order
  evals <- pmax(eg$values, 0)
  S <- sqrt(evals)
  ES <- X %*% eg$vectors
  # divide each column by S (avoid divide-by-zero)
  S_safe <- ifelse(S > 0, S, 1)
  E <- ES / rep(S_safe, each = nrow(X))
  E[, S == 0] <- 0
  list(E = E, S = S)
}

# Non-overlapping pair indices for significance testing
# (port of MATLAB `get_stat_indices`).
.crp_stat_indices <- function(N) {
  # base indices: every other element of 1:(N^2 - N)
  stat_indices <- seq.int(1L, N * N - N, by = 2L)
  if (N %% 2L == 1L) {
    # odd N: offset every even-indexed column block by 1
    b <- integer(length(stat_indices))
    half <- (N - 1L) %/% 2L
    for (k in seq_len(N)) {
      if (k %% 2L == 0L) {
        idx <- ((k - 1L) * half + 1L):(k * half)
        b[idx] <- 1L
      }
    }
    stat_indices <- stat_indices + b
  }
  stat_indices
}

# Per-trial outlier statistics (port of `trial_tests` from CRP_illustration.m).
.crp_trial_tests <- function(crp_projs, interval, num_trials) {
  S_test <- if (identical(interval, "tR")) {
    crp_projs$S_all[, crp_projs$tR_index]
  } else {
    crp_projs$S_all[, ncol(crp_projs$S_all)]
  }

  N <- num_trials
  L <- length(S_test) # = N^2 - N

  stat_indices <- .crp_stat_indices(N)
  excl_indices <- setdiff(seq_len(L), stat_indices)

  m <- numeric(N)
  p <- numeric(N)

  for (q in seq_len(N)) {
    # projections of normalized other trials into this trial, before this
    # cluster's own normalized projections (only present for q > 1)
    if (q > 1L) {
      t_indices1 <- seq.int(q - 1L, q * (N - 1L), by = N - 1L)
    } else {
      t_indices1 <- integer(0)
    }
    # projections of normalized other trials into this trial, after this
    # cluster's own normalized projections
    t_indices2 <- if (q * N <= L) {
      seq.int(q * N, L, by = N - 1L)
    } else {
      integer(0)
    }
    # projections of this trial into other trials
    t_indices3 <- (q - 1L) * (N - 1L) + seq_len(N - 1L)

    g_in <- unique(c(t_indices1, t_indices2, t_indices3))
    g_out <- setdiff(seq_len(L), c(g_in, excl_indices))

    # unpaired t-test (Welch by default in R; MATLAB ttest2 is pooled, but
    # the threshold is heuristic so this is acceptable).
    a <- S_test[g_in]
    b <- S_test[g_out]
    if (length(a) < 2L || length(b) < 2L ||
        stats::var(a) == 0 && stats::var(b) == 0) {
      p[q] <- NA_real_
    } else {
      tt <- tryCatch(
        stats::t.test(a, b),
        error = function(e) NULL
      )
      p[q] <- if (is.null(tt)) NA_real_ else tt$p.value
    }
    m[q] <- mean(S_test[t_indices3])
  }

  list(m = m, p = p)
}

# Core CRP computation (MATLAB `CRP_method`).
.crp_core <- function(V, t_win, srate, time_step) {
  Tn <- nrow(V)
  K <- ncol(V)
  proj_tpts <- seq.int(10L, Tn, by = time_step)
  if (length(proj_tpts) < 1L) {
    stop("Window too short for the requested `time_step`.")
  }

  npairs <- K * K - K
  S_all <- matrix(0, nrow = npairs, ncol = length(proj_tpts))
  m <- numeric(length(proj_tpts))
  v2 <- numeric(length(proj_tpts))

  # For each cumulative time point, calculate the projection
  for (i in seq_along(proj_tpts)) {
    k <- proj_tpts[i]
    S <- .crp_ccep_proj(V[seq_len(k), , drop = FALSE]) / sqrt(srate)
    S_all[, i] <- S
    m[i] <- mean(S)
    v2[i] <- stats::var(S)
  }

  # tt is the time point where the trial consistency reaches the maximum
  tt <- which.max(m)

  V_tR <- V[seq_len(proj_tpts[tt]), , drop = FALSE]
  pca <- .crp_kt_pca(V_tR)
  C <- pca$E[, 1L]
  # Eigenvector sign is arbitrary; orient C to align with the mean trace
  if (sum(C * rowMeans(V_tR)) < 0) { C <- -C }

  al <- as.numeric(crossprod(C, V_tR))
  ep <- V_tR - C %*% rbind(al)

  stat_indices <- .crp_stat_indices(K)

  s_tR <- S_all[stat_indices, tt]
  s_full <- S_all[stat_indices, ncol(S_all)]

  t_value_tR <- mean(s_tR) / (stats::sd(s_tR) / sqrt(length(s_tR)))
  p_value_tR <- stats::t.test(s_tR, mu = 0, alternative = "greater")$p.value

  t_value_full <- mean(s_full) / (stats::sd(s_full) / sqrt(length(s_full)))
  p_value_full <- stats::t.test(s_full, mu = 0, alternative = "greater")$p.value

  ep2 <- colSums(ep * ep)
  v2_tR <- colSums(V_tR * V_tR)

  crp_projs <- list(
    proj_tpts = t_win[proj_tpts],
    S_all = S_all,
    mean_proj_profile = m,
    var_proj_profile = v2,
    tR_index = tt,
    avg_trace_input = rowMeans(V),
    stat_indices = stat_indices,
    t_value_tR = t_value_tR,
    p_value_tR = p_value_tR,
    t_value_full = t_value_full,
    p_value_full = p_value_full
  )

  crp_params <- list(
    V_tR = V_tR,
    al = al,
    C = C,
    ep = ep,
    tR = t_win[proj_tpts[tt]],
    params_times = t_win[seq_len(proj_tpts[tt])],
    avg_trace_tR = rowMeans(V_tR),
    al_p = al / sqrt(length(C)),
    epep_root = sqrt(ep2),
    Vsnr = al / sqrt(ep2),
    expl_var = 1 - ep2 / v2_tR
  )

  list(crp_projs = crp_projs, crp_params = crp_params)
}
