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
#' @param detect_onset logical; if \code{TRUE}, estimate the response
#' \emph{onset} \eqn{\tau_{onset}} via a reverse cumulative-projection scan
#' (see \sQuote{Details}). The forward canonical shape and per-trial weights
#' are left unchanged; only the \emph{reported} canonical shape \code{C} is
#' restricted to the active response support \eqn{[\tau_{onset}, \tau_R]}.
#' The onset may fall before \code{t_start} (down to \code{onset_search_start}).
#' Defaults to \code{FALSE}.
#' @param onset_search_start numeric scalar or \code{NULL}; the earliest time
#' (in seconds) the onset scan may reach when \code{detect_onset = TRUE}.
#' \code{NULL} (the default) uses the first loaded time point, so the scan can
#' look across the entire signal before \eqn{\tau_R}. Raise it to keep the scan
#' away from an early stimulation artifact. Clamped into
#' \code{[min(time), t_end]}.
#'
#' @returns A named list with the following elements:
#' \describe{
#' \item{\code{parameters}}{A list of single-trial parameterizations
#' (\code{crp_parms} in MATLAB):
#' \describe{
#'   \item{\code{C}}{Numeric vector, the reported canonical response shape
#'     \eqn{C(t)}, taken as the slice of \code{C_full} over its active support
#'     (so it shares identical values with \code{C_full}). By default this is
#'     \eqn{[t_{start}, \tau_R]}, where it equals the first eigenvector of the
#'     linear kernel PCA on \code{V_tR} (oriented to the mean trace, unit-norm,
#'     length \eqn{T_R}). When \code{detect_onset = TRUE} it is restricted to
#'     \eqn{[\tau_{onset}, \tau_R]} (so its norm is then \eqn{\le 1}). The
#'     matching time axis is \code{params_times}.}
#'   \item{\code{C_full}}{Numeric vector spanning the \emph{entire} loaded
#'     time range (every row of \code{time}, including any baseline at
#'     \eqn{t < 0}). Obtained by applying the trial-space \code{loading} to the
#'     full retained data, so \code{C_full} coincides exactly with the
#'     untrimmed forward \code{C} on \eqn{[t_{start}, \tau_R]} and extrapolates
#'     the shape everywhere else. Time axis is \code{params_times_full}.
#'     Outside \eqn{[t_{start}, \tau_R]} cross-trial consistency is not
#'     optimized, so those portions are more variable.}
#'   \item{\code{loading}}{Numeric vector of length \eqn{K} (retained
#'     trials), the trial-space loading \eqn{g = V_{tR}^\top C = s_1 v_1}
#'     recovered from the forward decomposition. Applying it to a
#'     time \eqn{\times} trials matrix and dividing by \eqn{\|g\|^2}
#'     reconstructs the canonical shape over that time range; this is how
#'     \code{C_full} is formed.}
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
#'     after the shared component is removed, computed over the forward window
#'     \eqn{[t_{start}, \tau_R]} against the untrimmed forward canonical shape.
#'     Access trial \eqn{k} via \code{ep[, k]}. (When
#'     \code{detect_onset = TRUE} the reported \code{C} may be shorter than
#'     \eqn{T_R}; \code{ep} always stays on the full forward window.)}
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
#'   \item{\code{params_times}}{Numeric vector, time axis for the reported
#'     \code{C}; length matches \code{C} (\eqn{T_R}, or the onset-trimmed
#'     support when \code{detect_onset = TRUE}).}
#'   \item{\code{params_times_full}}{Numeric vector, time axis for
#'     \code{C_full}; equals the full loaded \code{time} (every row,
#'     including the baseline at \eqn{t < 0}).}
#'   \item{\code{V_tR}}{Numeric matrix \eqn{T_R \times K}, trial matrix
#'     truncated to \eqn{\tau_R} - the data actually decomposed. Together
#'     with \code{ep} and \code{avg_trace_tR} it always spans the full
#'     forward window \eqn{[t_{start}, \tau_R]}, independently of
#'     \code{detect_onset}.}
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
#'   \item{\code{tR_sample}}{Integer, row (sample) index of \eqn{\tau_R}
#'     within the windowed data; \code{V[seq_len(tR_sample), ]} is the
#'     truncated trial matrix \code{V_tR}.}
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
#' \item{\code{tau_onset}, \code{tau_onset_lower},
#' \code{tau_onset_upper}}{Numeric scalars, the estimated response onset
#' time (seconds) and its lower/upper threshold-crossing bounds, from the
#' reverse projection scan; all \code{NA} unless \code{detect_onset = TRUE}.}
#' \item{\code{onset}}{List with the reverse-scan profile (\code{onset_tpts}
#' and \code{mean_proj_profile}, mapped onto original time); \code{NULL}
#' unless \code{detect_onset = TRUE} (also \code{NULL} when the search window
#' had too few samples).}
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
#' When \code{detect_onset = TRUE}, a complementary reverse pass estimates the
#' response \emph{onset}. The retained trials between \code{onset_search_start}
#' and \eqn{\tau_R} are time-reversed, and the same cumulative cross-projection
#' profile is computed growing backward from \eqn{\tau_R}. The backward
#' duration that maximizes cross-trial consistency marks \eqn{\tau_{onset}},
#' the time at which trials begin to share structure - useful when the response
#' is delayed by an unknown latency. Because the scan can extend before
#' \code{t_start}, the onset may be earlier than the analysis window. This pass
#' only locates a time: it re-uses the forward loading, so the canonical shape
#' and per-trial weights are unchanged (the reported \code{C} is merely sliced
#' to \eqn{[\tau_{onset}, \tau_R]}).
#'
#' Time points whose data are not finite (any \code{NA}, \code{NaN} or
#' \code{Inf} across trials) are dropped before analysis. \code{t_start},
#' \code{t_end} and \code{onset_search_start} are clamped into the available
#' time range rather than triggering an error.
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
#' tt <- seq(-0.5, 1, length.out = n_time)
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
#' # ---- Onset detection (reverse scan) on the same data -------------------
#' res_onset <- crp(V, tt, detect_onset = TRUE)
#' c(tau_onset = res_onset$tau_onset, tau_R = res_onset$tau_R)
#'
#' @export
crp <- function(
  x, time,
  t_start = 0.015, t_end = 1,
  remove_artifacts = TRUE,
  artifact_interval = c("full", "tR"),
  artifact_p_threshold = 1e-5,
  threshold_quantile = 0.98,
  time_step = 5L,
  detect_onset = FALSE,
  onset_search_start = NULL
) {

  # ---- input validation ----------------------------------------------------
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric matrix with shape time x trials.")
  }
  if (ncol(x) < 2L) {
    stop("`x` must have at least 2 trials (columns).")
  }
  time <- as.numeric(time)
  if (anyNA(time) || length(time) != nrow(x)) {
    stop("`time` must be a non-NA numeric vector with length equal to nrow(x).")
  }

  # ---- drop time points with non-finite data (NA / NaN / Inf) -------------
  finite_rows <- is.finite(rowSums(x))
  if (!all(finite_rows)) {
    x <- x[finite_rows, , drop = FALSE]
    time <- time[finite_rows]
  }
  if (length(time) < 11L) {
    stop("`x` has too few finite time points (need >= 11 after removing ",
         "non-finite rows).")
  }
  if (any(diff(time) <= 0)) {
    stop("`time` must be strictly monotonically increasing.")
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
  detect_onset <- isTRUE(as.logical(detect_onset))

  # ---- clamp the analysis window into the available time range ------------
  # t_start / t_end only bound where the response may *start* / *end*; samples
  # outside the loaded `time` cannot constrain them, so clamp instead of
  # erroring (e.g. when non-finite rows trimmed the edges or t_end exceeds the
  # last sample).
  t_start <- max(t_start, min(time))
  t_end <- min(t_end, max(time))
  if (t_start >= t_end) {
    stop("`t_start` must be less than `t_end`.")
  }

  # onset search floor: how far before t_start the reverse scan may look
  if (is.null(onset_search_start)) {
    onset_search_start <- min(time)
  } else {
    onset_search_start <- as.numeric(onset_search_start)
    if (length(onset_search_start) != 1L || !is.finite(onset_search_start)) {
      stop("`onset_search_start` must be a single finite number or NULL.")
    }
    onset_search_start <- min(max(onset_search_start, min(time)), t_end)
  }

  # ---- subset to analysis window ------------------------------------------
  tpts <- which(time >= t_start & time <= t_end)
  if (length(tpts) < 11L) {
    stop("Analysis window contains too few samples (need >= 11).")
  }
  V <- x[tpts, , drop = FALSE]
  t_win <- time[tpts]
  srate <- 1 / mean(diff(time))

  # ---- optional artifact removal ------------------------------------------
  bad_trials <- integer(0)
  if (isTRUE(remove_artifacts)) {
    init <- crp_core(V, t_win, srate, time_step)
    mp <- crp_trial_tests(init$crp_projs, artifact_interval, ncol(V))
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
  retained_cols <- setdiff(seq_len(ncol(x)), bad_trials)

  # ---- final CRP pass ------------------------------------------------------
  out <- crp_core(V, t_win, srate, time_step)
  crp_projs <- out$crp_projs
  crp_params <- out$crp_params

  # ---- full-range canonical shape -----------------------------------------
  # Apply the trial-space loading to the entire loaded time range (every row of
  # the retained trials, including any baseline at t < 0). By construction
  # C_full equals the forward C exactly on [t_start, tau_R] and extrapolates the
  # shape everywhere else, so C_full, the reported C and the onset-trimmed C all
  # share identical values on their common support.
  x_ret <- x[, retained_cols, drop = FALSE]
  g <- crp_params$loading
  gg <- sum(g * g)
  C_full <- if (gg > 0) as.numeric(x_ret %*% g) / gg else numeric(nrow(x_ret))
  crp_params$C_full <- C_full
  crp_params$params_times_full <- time

  # ---- duration-uncertainty bounds (threshold_quantile of peak) ----------
  bnd <- crp_profile_bounds(crp_projs$mean_proj_profile, threshold_quantile)
  tau_R_lower <- crp_projs$proj_tpts[bnd$lb]
  tau_R_upper <- crp_projs$proj_tpts[bnd$hb]

  # ---- optional response-onset detection (reverse projection scan) -------
  # Re-uses the forward loading: the scan only *locates* the onset time and
  # never recomputes C or alpha. It may reach before t_start (down to
  # onset_search_start) so a delayed response with unknown latency can be found.
  tau_onset <- NA_real_
  tau_onset_lower <- NA_real_
  tau_onset_upper <- NA_real_
  onset <- NULL
  if (isTRUE(detect_onset)) {
    idx_tR_full <- tpts[crp_projs$tR_sample]
    idx_lo_full <- which(time >= onset_search_start)[1L]
    if (idx_tR_full - idx_lo_full + 1L >= 11L) {
      onset <- crp_onset(
        x_ret[idx_lo_full:idx_tR_full, , drop = FALSE],
        time[idx_lo_full:idx_tR_full],
        srate, time_step, threshold_quantile
      )
      tau_onset <- onset$tau_onset
      tau_onset_lower <- onset$tau_onset_lower
      tau_onset_upper <- onset$tau_onset_upper
    }
  }

  # ---- reported canonical shape (slice of C_full) -------------------------
  # Always derived from C_full so C_full, the reported C and the onset-trimmed C
  # share identical values on their common support. Defaults to
  # [t_start, tau_R]; restricted to [tau_onset, tau_R] when an onset was found.
  report_lo <- if (is.finite(tau_onset)) tau_onset else t_start
  sel <- crp_params$params_times_full >= report_lo &
    crp_params$params_times_full <= crp_params$tR
  crp_params$C <- crp_params$C_full[sel]
  crp_params$params_times <- crp_params$params_times_full[sel]

  structure(
    class = "ravetools_crp",
    list(
      parameters = crp_params,
      projections = crp_projs,
      bad_trials = bad_trials,
      tau_R_lower = tau_R_lower,
      tau_R = crp_params$tR,
      tau_R_upper = tau_R_upper,
      tau_onset = tau_onset,
      tau_onset_lower = tau_onset_lower,
      tau_onset_upper = tau_onset_upper,
      onset = onset,
      t_start = t_start,
      t_end = t_end,
      sample_rate = srate,
      .data = list(
        V = x_ret,
        time = time,
        time_range = range(time),
        window = c(t_start, t_end)
      )
    )
  )
}

#' @title Plot \verb{CRP} results
#' @description
#' S3 plot method for objects of class \code{ravetools_crp} returned by
#' \code{\link{crp}}.  Produces a three-panel figure:
#' \enumerate{
#' \item Single-trial traces over the entire loaded time range with the mean,
#'   the canonical shape \eqn{C(t)} on its active support (solid) and the
#'   full-range extension \eqn{C_{full}(t)} (dashed) overlaid; the analysis
#'   window edges and, when present, \eqn{\tau_{onset}} are marked.
#' \item Per-trial \eqn{\alpha'} weights sorted in ascending order.
#' \item Mean cross-trial projection profile with vertical lines marking
#'   \eqn{\tau_{lb}}, \eqn{\tau_R} and \eqn{\tau_{ub}}, plus the reverse
#'   onset profile and \eqn{\tau_{onset}} when \code{detect_onset} was used.
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
  window <- x$.data$window
  avg_trace <- if (is.matrix(V) && length(V) > 0) {
    rowMeans(V)
  } else {
    parms$C_full
  }

  op <- graphics::par(mfrow = c(1, 3))
  on.exit(graphics::par(op), add = TRUE)

  # ---- Panel 1: all trials (entire time range) + mean + C(t) overlay -----
  # Scale C_full (and the reported C, a slice of it) by the SAME factor so they
  # coincide on shared support; the factor is derived from the in-window region
  # for a stable overlay even though C_full is drawn across the whole range.
  inwin <- tt >= window[1L] & tt <= window[2L]
  cden <- max(abs(parms$C_full[inwin]), na.rm = TRUE)
  sf <- if (cden > 0) max(abs(avg_trace[inwin]), na.rm = TRUE) / cden else 1
  C_full_scaled <- parms$C_full * sf
  C_scaled <- parms$C * sf
  has_onset <- is.finite(x$tau_onset)

  ylim <- range(c(V, C_full_scaled, avg_trace), na.rm = TRUE)

  plot(range(tt), ylim, type = "n", las = 1,
       xlab = "Time (s)", ylab = expression(mu * V),
       main = expression("Canonical shape " * C(t)))

  if (is.matrix(V) && length(V) > 0) {
    # We might trim V to save space, hence V might not be available
    graphics::matlines(tt, V, lty = 1, col = "#80808060")
  }

  # analysis-window edges (where the forward tau_R search operated)
  graphics::abline(v = window, col = "#cccccc", lty = 3, lwd = 1)

  graphics::abline(v = x$tau_R, col = "red", lwd = 1)
  graphics::text(
    x = x$tau_R, y = min(C_full_scaled),
    labels = bquote(tau[R] * "=" * .(sprintf("%.3fs", x$tau_R))),
    col = "red", pos = 4, offset = 0.5
  )

  # full-range canonical shape (dashed = extrapolation outside active support)
  graphics::lines(parms$params_times_full, C_full_scaled,
                  col = "#FFB300", lty = 3, lwd = 3)
  # canonical shape on the active response support (solid)
  graphics::lines(parms$params_times, C_scaled, col = "#FFD500", lwd = 3)
  graphics::lines(tt, avg_trace, col = "black", lwd = 1)

  # onset marker (only when detect_onset was requested)
  if (has_onset) {
    graphics::abline(v = x$tau_onset, col = "forestgreen", lwd = 1, lty = 2)
    graphics::text(
      x = x$tau_onset, y = max(C_full_scaled),
      labels = bquote(tau[onset] * "=" * .(sprintf("%.3fs", x$tau_onset))),
      col = "forestgreen", pos = 4, offset = 0.5
    )
  }

  legend_lab <- c("mean", "C(t)", "C_full(t)")
  legend_col <- c("black", "#FFD500", "#FFB300")
  legend_lty <- c(1, 1, 3)
  graphics::legend("topright", legend = legend_lab,
                   col = legend_col, lty = legend_lty, lwd = 2,
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

  # ---- Panel 3: mean projection profile with tau_R (and onset) ----------
  onset <- x$onset
  if (has_onset && is.list(onset) &&
      length(onset$onset_tpts) > 0L) {
    has_onset_prof <- TRUE

    yr <- range(c(proj$mean_proj_profile, onset$mean_proj_profile), na.rm = TRUE)
    xr <- range(c(proj$proj_tpts, onset$onset_tpts), na.rm = TRUE)
  } else {
    has_onset_prof <- FALSE

    yr <- range(proj$mean_proj_profile, na.rm = TRUE)
    xr <- range(proj$proj_tpts, na.rm = TRUE)
  }

  graphics::plot(xr, yr, type = "n", xlab = "Candidate duration (s)", ylim = yr,
                 ylab = expression(bar(S) ~ (mu * V %.% s^{0.5})),
                 main = expression("Projection profile & " * tau[R]), las = 1)

  graphics::lines(proj$proj_tpts, proj$mean_proj_profile, lwd = 2)
  if (has_onset_prof) {
    # reverse (onset) projection profile, mapped back onto original time
    graphics::lines(onset$onset_tpts, onset$mean_proj_profile,
                    col = "gray", lwd = 2)
  }
  graphics::abline(v = x$tau_R, col = "#FF0000", lty = 1, lwd = 2)
  if (has_onset) {
    graphics::abline(v = x$tau_onset, col = "forestgreen", lwd = 2)
    leg <- c(
      bquote(tau[onset] * "=" * .(sprintf("%.3fs", x$tau_onset))),
      bquote(tau[R] * "=" * .(sprintf("%.3fs", x$tau_R)))
    )
    graphics::legend(
      "topright",
      legend = as.expression(leg),
      col = c("forestgreen", "red"),
      lwd = 2,
      bty = "n"
    )
  } else {
    leg <- c(bquote(tau[R] * "=" * .(sprintf("%.3fs", x$tau_R))))
    graphics::legend(
      "topright", legend = as.expression(leg),
      col = "red", lwd = 2, bty = "n")
  }

  invisible(x)
}


# --------------------------------------------------------------------------
# Internal helpers (not exported)
# --------------------------------------------------------------------------

# Single-trial cross-projection magnitudes (MATLAB `ccep_proj`).
# V is time x trials; returns a numeric vector of off-diagonal projections.
crp_ccep_proj <- function(V) {
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
crp_kt_pca <- function(X) {
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
crp_stat_indices <- function(N) {
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
crp_trial_tests <- function(crp_projs, interval, num_trials) {
  S_test <- if (identical(interval, "tR")) {
    crp_projs$S_all[, crp_projs$tR_index]
  } else {
    crp_projs$S_all[, ncol(crp_projs$S_all)]
  }

  N <- num_trials
  L <- length(S_test) # = N^2 - N

  stat_indices <- crp_stat_indices(N)
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
        error = function(e) { NULL }
      )
      p[q] <- if (is.null(tt)) NA_real_ else tt$p.value
    }
    m[q] <- mean(S_test[t_indices3])
  }

  list(m = m, p = p)
}

# Cumulative cross-trial projection profile over candidate truncation lengths.
# Sweeps a truncation length `k` over the rows (time) of `V`; for each `k` it
# computes the off-diagonal cross-trial projection magnitudes and their mean and
# variance. Shared by the forward CRP pass and the reverse onset scan (which
# passes a time-reversed `V`). Indices are returned in the sample domain so the
# caller can map them back onto its own time axis.
crp_proj_profile <- function(V, srate, time_step) {
  Tn <- nrow(V)
  K <- ncol(V)
  proj_idx <- seq.int(10L, Tn, by = time_step)
  if (length(proj_idx) < 1L) {
    stop("Window too short for the requested `time_step`.")
  }

  npairs <- K * K - K
  S_all <- matrix(0, nrow = npairs, ncol = length(proj_idx))
  m <- numeric(length(proj_idx))
  v2 <- numeric(length(proj_idx))

  # For each cumulative time point, calculate the projection
  for (i in seq_along(proj_idx)) {
    k <- proj_idx[i]
    S <- crp_ccep_proj(V[seq_len(k), , drop = FALSE]) / sqrt(srate)
    S_all[, i] <- S
    m[i] <- mean(S)
    v2[i] <- stats::var(S)
  }

  list(
    proj_idx = proj_idx,
    S_all = S_all,
    mean_proj_profile = m,
    var_proj_profile = v2
  )
}

# Core CRP computation (MATLAB `CRP_method`).
crp_core <- function(V, t_win, srate, time_step) {
  K <- ncol(V)

  # Cumulative cross-trial projection profile (shared with the onset scan)
  prof <- crp_proj_profile(V, srate, time_step)
  proj_tpts <- prof$proj_idx
  S_all <- prof$S_all
  m <- prof$mean_proj_profile
  v2 <- prof$var_proj_profile

  # tt is the time point where the trial consistency reaches the maximum
  tt <- which.max(m)
  tR_sample <- proj_tpts[tt]

  V_tR <- V[seq_len(tR_sample), , drop = FALSE]
  pca <- crp_kt_pca(V_tR)
  C <- pca$E[, 1L]
  # Eigenvector sign is arbitrary; orient C to align with the mean trace
  if (sum(C * rowMeans(V_tR)) < 0) { C <- -C }

  al <- as.numeric(crossprod(C, V_tR))
  ep <- V_tR - C %*% rbind(al)

  # Trial-space loading recovered from the unit-norm, sign-oriented C and V_tR
  # (g = s1 * v1). Applying it to any time x trials matrix and dividing by
  # ||g||^2 reconstructs the canonical shape there; crp() uses it to build the
  # full-range C_full. Kept in trial space (length K) so it generalizes beyond
  # the analysis window.
  loading <- as.numeric(crossprod(V_tR, C))   # length K (retained trials)

  stat_indices <- crp_stat_indices(K)

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
    tR_sample = tR_sample,
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
    loading = loading,
    ep = ep,
    tR = t_win[tR_sample],
    params_times = t_win[seq_len(tR_sample)],
    avg_trace_tR = rowMeans(V_tR),
    al_p = al / sqrt(length(C)),
    epep_root = sqrt(ep2),
    Vsnr = al / sqrt(ep2),
    expl_var = 1 - ep2 / v2_tR
  )

  list(crp_projs = crp_projs, crp_params = crp_params)
}

# Threshold-crossing bounds around the peak of a projection profile.
# Returns the peak index and the lower/upper indices where `m` last (going up)
# and first (coming down) crosses `threshold_quantile * max(m)` on either side
# of the peak. Shared by the forward tau_R bounds and the reverse onset bounds.
crp_profile_bounds <- function(m, threshold_quantile) {
  n <- length(m)
  cidx <- which.max(m)
  thr <- threshold_quantile * m[cidx]

  up_cross <- 1L + which(m[-1L] > thr & m[-n] <= thr)
  up_cross <- up_cross[up_cross <= cidx]
  lb <- if (length(up_cross) > 0L) up_cross[length(up_cross)] else 1L

  down_cross <- which(m[-n] > thr & m[-1L] <= thr)
  down_cross <- down_cross[down_cross >= cidx]
  hb <- if (length(down_cross) > 0L) down_cross[1L] else n

  list(peak = cidx, lb = lb, hb = hb)
}

# Response-onset estimation by reverse cumulative projection.
# `Xo` is the retained trial matrix over [onset_search_start, tau_R] in forward
# time order, with matching times `to` (so the last row is tau_R). The trials
# are time-reversed, so the cumulative truncation grows backward from tau_R; the
# same cross-trial projection profile is computed and the backward duration that
# maximizes cross-trial consistency marks the response onset tau_onset, with the
# threshold crossings giving its uncertainty bounds. The scan can reach before
# t_start. It only locates a time - it never recomputes the canonical shape C or
# the weights alpha, which are taken from the forward pass.
crp_onset <- function(Xo, to, srate, time_step, threshold_quantile) {
  n <- nrow(Xo)

  # too few samples to run a backward sweep
  if (n < 11L) {
    return(list(
      tau_onset = NA_real_,
      tau_onset_lower = NA_real_,
      tau_onset_upper = NA_real_,
      onset_tpts = numeric(0),
      mean_proj_profile = numeric(0)
    ))
  }

  V_rev <- Xo[n:1L, , drop = FALSE]
  prof <- crp_proj_profile(V_rev, srate, time_step)
  m <- prof$mean_proj_profile
  bnd <- crp_profile_bounds(m, threshold_quantile)

  # map backward sample extents to forward-time rows of `Xo`/`to`: a backward
  # extent of `k` samples from tau_R (row n) reaches forward row (n - k + 1). A
  # larger extent reaches an earlier time, so the earliest (lower) onset bound
  # uses the upper extent hb and the latest (upper) bound uses the lower extent.
  onset_tpts <- to[n - prof$proj_idx + 1L]
  ord <- order(onset_tpts)

  list(
    tau_onset = to[n - prof$proj_idx[bnd$peak] + 1L],
    tau_onset_lower = to[n - prof$proj_idx[bnd$hb] + 1L],
    tau_onset_upper = to[n - prof$proj_idx[bnd$lb] + 1L],
    onset_tpts = onset_tpts[ord],
    mean_proj_profile = m[ord]
  )
}
