#' @title Common Average Re-referencing by Least Anti-Correlation (CARLA)
#' @description
#' Selects an optimal subset of channels to use as the common average reference
#' (CAR) for cortico-cortical evoked potential (\verb{CCEP}) data, following the
#' \verb{CARLA} (see 'Reference' and 'Citation'). Channels are ranked in
#' increasing order of their cross-trial covariance (or variance when only one
#' trial is available); subsets are then iteratively grown and the size that
#' yields the least anti-correlation between the candidate reference and the
#' remaining unreferenced channels is selected as optimal.
#'
#' @param x numeric array of shape \code{channels x time x trials}; if a
#' matrix is provided, it is treated as \code{channels x time} (single trial).
#' The signal must already be cropped to the responsive time window
#' (typically \code{0.01 s} to \code{0.3 s} post-stimulus).
#' @param nboot integer, number of bootstrapped trial resamplings used to
#' estimate the optimization statistic; defaults to \code{100}. Ignored when
#' a single trial is supplied.
#' @param sensitive logical; if \code{TRUE} (and more than one trial is
#' supplied), the more sensitive cutoff is used, corresponding to the channel
#' count just before the first statistically significant decrease in the mean
#' anti-correlation curve. The default \code{FALSE} returns the global maximum
#' of that curve.
#' @param min_size integer, minimum subset size considered when
#' \code{sensitive = TRUE}; defaults to
#' \code{max(2, ceiling(0.1 * nchan_good))}, where \code{nchan_good} is
#' the number of usable channels after the bad-channel mask has been
#' applied (see Details). Floored at 2 because subset size 1 is never
#' evaluated.
#' @param absolute_rank logical; if \code{FALSE} (the default), the
#' per-channel ranking statistic is the signed mean of the off-diagonal
#' trial-to-trial covariances, as in the original CARLA manuscript. This
#' assumes a phase-locked evoked response, so trial-pair covariances on
#' responsive channels are positive and a clean signed mean separates them
#' from non-responsive channels. Set to \code{TRUE} to use the mean of the
#' \strong{absolute} covariances instead, which is more robust to evoked
#' responses whose polarity flips across trials (e.g. alternating-polarity
#' stimulation, biphasic \verb{CCEPs} with jitter), at the cost of an upward
#' bias on the non-responsive floor. Ignored when only one trial is supplied
#' (the per-channel variance is used in that case).
#' @param virtual_reference logical; if \code{TRUE}, runs the modified
#' CARLA before the iterative subset
#' evaluation: rank channels once, designate the channel with the median
#' rank as a "virtual reference", subtract its signal from every channel,
#' then re-rank on the subtracted data. The virtual channel itself is
#' pinned to the very end of the new ordering so it is never picked into
#' the CAR until the final subset size. This is intended for data where
#' the recording reference is contaminated by stimulation artifact or
#' evoked activity: a contaminated reference makes every channel look
#' artificially similar and biases the covariance/variance ranking statistic,
#' but subtracting a mid-rank proxy of that contamination unbiases the
#' ranking so genuinely responsive channels rise to the top. Default
#' \code{FALSE} reproduces the original CARLA algorithm bit-for-bit.
#'
#' @returns A list with the following elements:
#' \describe{
#' \item{\code{channels}}{integer vector, sorted indices (1-based) of the
#' channels chosen to construct the common average reference.}
#' \item{\code{car}}{numeric matrix of shape \code{time x trials} (or numeric
#' vector when only one trial is supplied) holding the common average signal
#' computed from \code{channels}. Subtract this from each channel of the
#' original signal to obtain the re-referenced data.}
#' \item{\code{order}}{integer vector, indices of the \strong{good}
#' channels sorted in increasing order of the ranking statistic. Bad
#' channels (zero variance / all-\code{NA}; see Details) are excluded.}
#' \item{\code{vars}}{numeric vector of length \code{nchan} containing the
#' per-channel ranking statistic (mean cross-trial covariance, or variance
#' for a single trial). Bad channels keep their raw value (zero or
#' \code{NA}) so callers can audit the mask.}
#' \item{\code{n_optimum}}{integer, the optimal subset size selected
#' (indexes into \code{order}).}
#' \item{\code{zmin_mean}}{numeric matrix of shape
#' \code{length(order) x nboot} (or a length-\code{length(order)} vector
#' for a single trial / \code{nboot = 1}) holding, for each subset size
#' and bootstrap, the mean Fisher z-transformed correlation of the most
#' globally anti-correlated unreferenced channel against the candidate CAR.
#' Row 1 is always \code{NA} (subset size 1 is not evaluated).}
#' \item{\code{bad_channels}}{integer vector of channel indices that were
#' excluded from the analysis because their ranking statistic was zero or
#' \code{NA} (flat / dead channels, plus the virtual channel itself when
#' \code{virtual_reference = TRUE}).}
#' \item{\code{virtual_channel}}{integer, the index (1-based) of the
#' channel used as the virtual reference when \code{virtual_reference =
#' TRUE}; \code{NA_integer_} otherwise.}
#' \item{\code{vars1}}{numeric vector of the first-pass ranking statistic
#' when \code{virtual_reference = TRUE} (the post-subtraction statistic is
#' returned in \code{vars}); \code{NULL} otherwise.}
#' }
#'
#' @details
#' The function is a faithful port of the core \code{CARLA.m} routine; it does
#' not perform notch filtering, time-window cropping,
#' or grouping by stimulation site. Those steps belong to the surrounding
#' preprocess pipeline (see the example below).
#'
#' For each candidate subset size \eqn{n = 2, \ldots, N}{n = 2, ..., N}, the
#' candidate reference is computed as the channel-wise mean of the \eqn{n}
#' lowest-ranked channels. Each of those \eqn{n} channels is then correlated,
#' in its unreferenced form, against every channel of the candidate
#' re-referenced subset. The resulting Pearson correlations are
#' Fisher z-transformed; the row corresponding to the most globally
#' anti-correlated channel is recorded as \code{zmin}. The optimal \eqn{n}
#' is the one that maximizes (i.e. makes least negative) the mean of
#' \code{zmin}.
#'
#' \strong{Bad-channel mask.} Channels whose ranking statistic is exactly
#' zero or \code{NA} are flagged as "bad" and excluded both from the CAR
#' candidate pool and from the target-channel correlation rows. This
#' catches flat / dead channels (constant signal, no usable variance) and
#', when \code{virtual_reference = TRUE}, automatically removes the virtual
#' channel itself, since subtracting it from itself yields a zero trace.
#' All references to channel indices in the returned \code{order},
#' \code{n_optimum}, and \code{zmin_mean} are with respect to the
#' \strong{good} channels only; \code{vars} is full-length so callers can
#' inspect the raw scores.
#'
#' @references
#' The CARLA algorithm (\code{virtual_reference = FALSE}) is described in
#' \doi{10.1016/j.jneumeth.2024.110153}; the modified CARLA precursor
#' (\code{virtual_reference = TRUE}) is described in
#' \doi{10.1016/j.jneumeth.2025.110461}, with a reference implementation
#' at \url{https://github.com/hharveygit/SPES_reference_contam}. See
#' \code{citation("ravetools")} for the full bibliographic entries of
#' both manuscripts.
#'
#' @examples
#'
#' # ---- Simulate a small CCEP-like dataset --------------------------------
#' # 16 channels, 12 trials, sampled at 1 kHz, 0.5 s peri-stimulus epoch.
#' # Channels 1:4 are "responsive" (carry an evoked potential); the rest are
#' # noise-only and should make up the optimal CAR.
#' srate <- 1000
#' tt <- seq(-0.1, 0.4 - 1 / srate, by = 1 / srate)   # time, seconds
#' nchan <- 16
#' ntrial <- 30
#' resp_ch <- 1:4
#' noise_ch <- 4:7
#'
#' # Evoked potential template: damped sinusoid starting at t = 0
#' ep <- ifelse(tt >= 0,
#'              80 * exp(-tt / 0.05) * sin(2 * pi * 12 * tt),
#'              0)
#'
#' # channels x time x trials
#' x_full <- array(rnorm(nchan * length(tt) * ntrial, sd = 5),
#'                 dim = c(nchan, length(tt), ntrial))
#' for (ch in resp_ch) {
#'   for (k in seq_len(ntrial)) {
#'     x_full[ch, , k] <- x_full[ch, , k] + ep * runif(1, -0.8, 1.2)
#'   }
#' }
#' for (ch in noise_ch) {
#'   for (k in seq_len(ntrial)) {
#'     tmp <- x_full[ch, , k]
#'     x_full[ch, , k] <- tmp + sign(ch %% 2 - 0.5) * 5 *
#'       runif(length(tmp), 0.8, 1.2)
#'   }
#' }
#' # Add artifacts common to all channels and trials
#' artifacts <- 6 * sin(2 * pi * 60 * tt) + 7 * sin(2 * pi * 24 * tt)
#' x_full <- sweep(x_full, 2L, artifacts, "+")
#'
#' # ---- 1. Notch filter line noise (per channel, per trial) ---------------
#' # The CARLA paper notch-filters before ranking; the re-reference itself
#' # is applied to the original (unfiltered) signal.
#' x_clean <- x_full
#' for (ch in seq_len(nchan)) {
#'   for (k in seq_len(ntrial)) {
#'     x_clean[ch, , k] <- notch_filter(
#'       x_full[ch, , k], sample_rate = srate,
#'       lb = c(59, 119, 179), ub = c(61, 121, 181)
#'     )
#'   }
#' }
#'
#' # ---- 2. Crop to the responsive window (0.01 s to 0.3 s post-stim) ------
#' resp_idx <- which(tt > 0.0 & tt <= 0.3)
#' x_resp <- x_clean[, resp_idx, , drop = FALSE]
#'
#' # ---- 3. Run CARLA to pick reference channels ---------------------------
#' fit <- carla(x_resp, sensitive = TRUE, absolute_rank = TRUE,
#'              virtual_reference = TRUE)
#' fit$channels        # selected reference channels (should exclude 1:4)
#' fit$n_optimum       # number of channels in the optimal CAR
#'
#' # ---- 4. Re-reference the ORIGINAL (unfiltered) signal ------------------
#' # mean or median, your choice!
#' car_full <- apply(x_full[fit$channels, , , drop = FALSE], c(2, 3), mean)
#'
#' # old-style: using all channels for CAR
#' car_old <- apply(x_full, c(2, 3), mean)
#'
#' x_reref <- sweep(x_full, c(2, 3), car_full, "-")
#' x_compare <- sweep(x_full, c(2, 3), car_old, "-")
#'
#' # ---- 5. Inspect: evoked potential is preserved on responsive channels --
#' op <- graphics::par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
#' ravetools::plot_signals(
#'   signals = x_full[, , 1],
#'   sample_rate = srate,
#'   main = "Trial 1 - (raw)")
#'
#' ravetools::plot_signals(
#'   signals = x_clean[, , 1],
#'   sample_rate = srate,
#'   main = "Notch-filtered")
#'
#' ravetools::plot_signals(
#'   signals = x_reref[, , 1],
#'   sample_rate = srate,
#'   main = sprintf("CARLA-ref (n=%d)", length(fit$channels)))
#'
#' ravetools::plot_signals(
#'   signals = x_compare[, , 1],
#'   sample_rate = srate,
#'   main = "Conventional CAR for comparison")
#'
#' col <- adjustcolor(seq_len(nchan))
#' col[resp_ch] <- adjustcolor(col[resp_ch], alpha.f = 0.2)
#'
#' graphics::matplot(tt, t(rowMeans(x_full, dims = 2L)),
#'                   type = "l", lty = 1, xlab = "Time (s)", ylab = "uV",
#'                   main = "Trial-averaged (raw", col = col)
#'
#' graphics::matplot(tt, t(rowMeans(x_clean, dims = 2L)),
#'                   type = "l", lty = 1, xlab = "Time (s)", ylab = "uV",
#'                   main = "Notch-filtered", col = col)
#'
#' graphics::matplot(tt, t(rowMeans(x_reref, dims = 2L)),
#'                   type = "l", lty = 1, xlab = "Time (s)", ylab = "uV",
#'                   main = "CARLA-referenced", col = col)
#'
#' graphics::matplot(tt, t(rowMeans(x_compare, dims = 2L)),
#'                   type = "l", lty = 1, xlab = "Time (s)", ylab = "uV",
#'                   main = "conventional CAR-referenced", col = col)
#'
#' graphics::par(op)
#'
#' @export
carla <- function(x, nboot = 100L, sensitive = FALSE, min_size = NULL,
                  absolute_rank = FALSE, virtual_reference = FALSE) {

  stopifnot(is.numeric(x))
  if (is.matrix(x)) {
    dim(x) <- c(dim(x), 1L)
  }
  d <- dim(x)
  if (length(d) != 3L) {
    stop("`x` must be a numeric matrix (channels x time) or 3-D array ",
         "(channels x time x trials).")
  }
  n_ch <- d[[1L]] # Number of channels
  n_t  <- d[[2L]] # Number of time points
  n_tr <- d[[3L]] # Number of trials
  if (n_ch < 2L) {
    stop("`x` must have at least 2 channels.")
  }
  if (n_t < 2L) {
    stop("`x` must have at least 2 time points in the response window.")
  }

  nboot <- as.integer(nboot)
  if (is.na(nboot) || nboot < 1L) {
    stop("`nboot` must be a positive integer.")
  }
  if (n_tr == 1L) {
    nboot_eff <- 1L
  } else {
    nboot_eff <- nboot
  }

  # ---- 1. Per-channel ranking statistic --------------------------------------
  # In the original manuscript there is no `abs()`: the assumption is a
  # phase-locked evoked response, so trial-pair covariances are positive
  # and their signed mean is a clean score. That signed-mean form is the
  # default here (`absolute_rank = FALSE`). Setting `absolute_rank = TRUE`
  # takes the magnitude instead, which is robust to polarity-flipping
  # evoked responses (e.g. alternating-polarity stimulation, biphasic
  # CCEPs) where signed covariances can cancel.
  cov_summary <- if (isTRUE(absolute_rank)) {
    function(v) mean(abs(v), na.rm = TRUE)
  } else {
    function(v) mean(v, na.rm = TRUE)
  }
  rank_fn <- function(signal) {
    if (n_tr == 1L) {
      apply(signal, 1L, function(row) stats::var(as.numeric(row)))
    } else {
      vapply(seq_len(n_ch), function(ii) {
        slice <- matrix(signal[ii, , ], nrow = n_t, ncol = n_tr)  # time x trials
        cv <- stats::cov(slice)                                   # trials x trials
        cov_summary(cv[upper.tri(cv)])
      }, numeric(1L))
    }
  }

  vars <- rank_fn(x)

  # Signal model (denote i: channel, j: trial, t: time):
  #   signal_ij(t) =
  #       r_i(t)        <- evoked response, shared across trials of channel i
  #                        (zero or near-zero on non-responsive channels)
  #     + a(t)          <- global artifact / common-mode signal we want to
  #                        remove via re-referencing (line noise, stim leak,
  #                        ground/ref drift); shared across channels
  #     + e_ij(t)       <- zero-mean noise, independent across trials and
  #                        (approximately) across channels
  #
  # The ranking statistic `vars[i]` estimates the across-trial covariance of
  # channel i on the response window. Taking the expectation:
  #   Cov(signal_i,j1, signal_i,j2)  ~=  Var(r_i) + Var(a)
  # because r_i and a are constant across trials while e_ij averages out.
  # Therefore:
  #   - Non-responsive channels  -> vars[i] ~= Var(a)            (small, common floor)
  #   - Responsive channels      -> vars[i] ~= Var(r_i) + Var(a) (larger)
  # Sorting `vars` in increasing order puts the least-responsive channels
  # first; growing the CAR subset from the top of `ord` keeps the evoked
  # response out of the reference for as long as possible.

  # ---- 1b. Optional virtual-reference pre-step (modified CARLA) -------------
  # Huang et al. (2025): when the recording reference is contaminated by
  # stim/evoked activity, every channel inherits a copy of the
  # contamination, which makes the ranking statistic above biased --
  # responsive and non-responsive channels look more similar than they are.
  # The fix is to pick the channel whose first-pass rank is in the middle
  # (so it is unlikely to be either pure noise or a strongly responsive
  # channel, but should still carry the contamination), subtract it from
  # every channel, and re-rank on the subtracted signal. The virtual
  # channel itself becomes identically zero after subtraction, so we
  # remove it from the new order and append it at the very end: it is
  # never picked into the CAR until the final subset size and never
  # influences any of the correlations downstream.
  if (virtual_reference) {
    vars1 <- vars
    # MATLAB `round(nChs/2)` rounds halves AWAY from zero; R's base
    # `round` is banker's (half-to-even). Use `floor(n/2 + 0.5)` to match
    # the reference implementation exactly.
    virtual_channel <- order(vars1)[as.integer(floor(n_ch / 2 + 0.5))]
    # Force n_t x n_tr matrix even when n_tr == 1 (drop=TRUE would
    # collapse to a length-n_t vector and break sweep over c(2, 3)).
    vref <- matrix(x[virtual_channel, , ], nrow = n_t, ncol = n_tr)
    x_work <- sweep(x, c(2L, 3L), vref, "-")
    # The virtual channel is now identically zero in `x_work`, so the
    # bad-channel mask below will (correctly) drop it from `ord`. There is
    # no need for the explicit `c(ord[ord != vc], vc)` reordering trick.
    vars <- rank_fn(x_work)
  } else {
    vars1 <- NULL
    virtual_channel <- NA_integer_
    x_work <- x
  }

  # ---- 1c. Bad-channel mask --------------------------------------------------
  # Flat / dead channels (zero variance) and channels that produced an NA
  # ranking statistic (e.g. all-NA input rows) cannot meaningfully appear in
  # the CAR - their correlation rows would be NaN - and including them
  # introduces noise into the bootstrap statistic. Drop them from `ord` and
  # shrink the iterative loop accordingly.
  bad_mask <- is.na(vars) | vars == 0
  bad_channels <- which(bad_mask)
  ord <- order(vars)
  ord <- ord[!bad_mask[ord]]
  n_good <- length(ord)
  if (n_good < 2L) {
    stop("`carla` needs at least 2 good channels (found ", n_good,
         " after dropping flat / NA channels).")
  }

  # ---- 2. Iterative subset evaluation ----------------------------------------
  # zmin <- array(NA_real_, dim = c(n_ch, n_ch, nboot_eff))
  #
  # for (ii in seq.int(2L, n_ch)) {
  #   sel <- ord[seq_len(ii)]
  #   sub <- x[sel, , , drop = FALSE]                # ii x n_t x n_tr
  #   car_sub <- colMeans(sub)                       # n_t x n_tr (drops dim 1)
  #   if (n_tr == 1L) {
  #     dim(car_sub) <- c(n_t, 1L)
  #   }
  #   sub_reref <- sweep(sub, c(2L, 3L), car_sub, "-")
  #
  #   if (n_tr == 1L) {
  #     Useg   <- t(matrix(sub[, , 1L],       nrow = ii, ncol = n_t))
  #     Useg_r <- t(matrix(sub_reref[, , 1L], nrow = ii, ncol = n_t))
  #     r <- stats::cor(Useg, Useg_r)
  #     diag(r) <- NA_real_
  #     z <- atanh(r)
  #     kk <- which.min(rowMeans(z, na.rm = TRUE))
  #     zmin[sel, ii, 1L] <- z[kk, ]
  #     next
  #   }
  #
  #   for (b in seq_len(nboot_eff)) {
  #     inds <- sample.int(n_tr, n_tr, replace = TRUE)
  #     Useg_m <- rowMeans(sub[, , inds, drop = FALSE],       dims = 2L)  # ii x n_t
  #     Uref_m <- rowMeans(sub_reref[, , inds, drop = FALSE], dims = 2L)
  #     r <- stats::cor(t(Useg_m), t(Uref_m))
  #     diag(r) <- NA_real_
  #     z <- atanh(r)
  #     kk <- which.min(rowMeans(z, na.rm = TRUE))
  #     zmin[sel, ii, b] <- z[kk, ]
  #   }
  # }

  # Per-`ii` callback returns the length-`nboot_eff` vector of bootstrap
  # zmin_mean values for that subset size: i.e. for each bootstrap b, the
  # mean over target channels of `z[kk, ]` for the worst-contaminator row
  # `kk`. Collapsing the target-channel dim here (rather than keeping a
  # full `n_ch x nboot_eff` slab and averaging later) avoids allocating a
  # `n_ch x nboot_eff x (n_ch - 1)` array we never read.
  zmin_per_ii <- vapply(seq.int(2L, n_good), function(ii) {
    # evaluate as if the referencing channels are the first ii least important channels
    sel <- ord[seq_len(ii)]

    sub <- x_work[sel, , , drop = FALSE]           # ii x n_t x n_tr
    car_sub <- colMeans(sub)                       # n_t x n_tr (drops dim 1)

    if (n_tr == 1L) {
      dim(car_sub) <- c(n_t, 1L)
    }

    # Compute the reference based on these ii channels
    sub_reref <- sweep(sub, c(2L, 3L), car_sub, "-")

    if (nboot_eff == 1L) {
      # Single-pass: either a single trial was supplied, or the user asked
      # for `nboot = 1` with multiple trials. In either case there is no
      # bootstrap; we just collapse trials by simple mean (a no-op when
      # n_tr == 1) and run one correlation pass.
      Useg_m <- rowMeans(sub,       dims = 2L)     # ii x n_t
      Uref_m <- rowMeans(sub_reref, dims = 2L)
      r <- stats::cor(t(Useg_m), t(Uref_m))
      diag(r) <- NA_real_
      z <- atanh(r)
      kk <- which.min(rowMeans(z, na.rm = TRUE))

      # length nboot_eff (= 1)
      mean(z[kk, ], na.rm = TRUE)
    } else {

      # Repeat nboot times: ii x nboot_eff
      sub_zmin <- replicate(nboot_eff, {

        # Bootstrap trials
        inds <- sample.int(n_tr, n_tr, replace = TRUE)

        # Resample trials
        Useg_m <- rowMeans(sub[, , inds, drop = FALSE],       dims = 2L)  # ii x n_t
        Uref_m <- rowMeans(sub_reref[, , inds, drop = FALSE], dims = 2L)

        # Compute the correlation between raw vs referenced signals
        # Signal model (denote i: channel, j: trial, t: time):
        #   signal_ij(t) =
        #       r_i(t)        <- evoked response, shared across trials of channel i
        #                        (zero or near-zero on non-responsive channels)
        #     + a(t)          <- global artifact / common-mode signal we want to
        #                        remove via re-referencing (line noise, stim leak,
        #                        ground/ref drift); shared across channels
        #     + e_ij(t)       <- zero-mean noise, independent across trials and
        #                        (approximately) across channels
        #
        # Bootstrap-averaging across `n_tr` trials suppresses e_ij, so the
        # rows of `Useg_m` and `Uref_m` are approximately
        #   Useg_m[k, t] ~= r_{sel[k]}(t) + a(t)
        #   Uref_m[k, t] ~= r_{sel[k]}(t) - mean_{m in sel} r_m(t)
        #                   (the shared a(t) cancels in the re-reference)
        #
        # Then `stats::cor(t(Useg_m), t(Uref_m))` ~=
        #   r[k, l] = Cor_t( r_{sel[k]} + a,
        #                    r_{sel[l]} - mean_{m in sel} r_m )
        # an `ii x ii` matrix whose row k summarizes how much the unreferenced
        # channel `sel[k]` co-varies (or anti-co-varies) with every other
        # re-referenced channel in the candidate subset.
        #
        # Interpretation of off-diagonal entry r[k, l]  (k != l):
        #   - Large POSITIVE: channel `sel[k]` shares a strong waveform with
        #     channel `sel[l]` AFTER `sel[l]` has been re-referenced. Because
        #     a(t) has been removed from `Uref_m`, the only structure left to
        #     drive a positive correlation is the evoked response r_{sel[l]}
        #     (or, more precisely, r_{sel[l]} minus the subset's mean response),
        #     which means `sel[k]` itself carries a similar response. This is
        #     the signature of two responsive channels being mistakenly placed
        #     into the CAR subset.
        #   - Large NEGATIVE: channel `sel[k]` looks like the negative of the
        #     re-referenced waveform of `sel[l]`. This happens when the CAR
        #     constructed from the subset already contains a sizeable copy of
        #     r_{sel[k]}: subtracting that CAR from `sel[l]` flips a piece of
        #     `sel[k]`'s response into `sel[l]`'s residual with the opposite
        #     sign, so the unreferenced `sel[k]` and the re-referenced `sel[l]`
        #     are anticorrelated. Anticorrelation is therefore the diagnostic
        #     that some channel in the subset is leaking response structure
        #     into the common average.
        #   - Near ZERO: the two channels share neither a residual response
        #     nor a CAR-mediated cross-talk; this is the desired regime for
        #     channels that legitimately belong in the reference (their only
        #     shared content was a(t), which has been cancelled).
        #
        # Within each row k, taking `mean_l r[k, l]` (excluding the diagonal)
        # collapses these l-wise relationships into a single anticorrelation
        # score for channel `sel[k]`; `kk` is then the channel with the most
        # negative such score, i.e. the worst contaminator under the current
        # candidate subset.
        #
        # The self-self diagonal
        # is dropped (set to NA) because it is dominated by the trivial
        # contribution of `sel[k]` to its own re-reference and would bias the
        # row mean toward zero.
        #
        # Fisher z-transform stabilizes the variance of the correlation
        # estimator before averaging. The row with the most negative mean
        # (`kk`) identifies the channel whose presence in the subset induces
        # the strongest global anticorrelation in the rest of the subset --
        # i.e. the channel whose responses are most contaminating the
        # candidate reference. Recording z[kk, ] in `zmin` lets the outer
        # loop track, as `ii` grows, when the next added channel starts
        # injecting evoked-response structure into the CAR; that transition
        # is what `n_optimum` ultimately picks up on.

        r <- stats::cor(t(Useg_m), t(Uref_m))
        diag(r) <- NA_real_
        z <- atanh(r)
        kk <- which.min(rowMeans(z, na.rm = TRUE))
        z[kk, ]
      })
      # length nboot_eff: mean over target channels (sel rows) of the
      # worst-contaminator z-row, per bootstrap.
      colMeans(sub_zmin, na.rm = TRUE)
    }
  }, FUN.VALUE = numeric(nboot_eff))

  # `vapply` with a numeric(nboot_eff) FUN.VALUE returns shape
  #   (boot, ii_index)   where ii_index runs 1..(n_ch - 1)
  # corresponding to subset sizes 2..n_ch. Transpose to the (ii, boot)
  # layout used downstream and prepend an all-NA row so the ii axis spans
  # 1..n_ch (size 1 is meaningless and never selected by `which.max`).
  if (nboot_eff == 1L) {
    # vapply drops to a length-(n_ch - 1) vector when FUN.VALUE has length 1
    zmin_mean <- matrix(c(NA_real_, zmin_per_ii), ncol = 1L)
  } else {
    zmin_mean <- rbind(NA_real_, t(zmin_per_ii), deparse.level = 0L)
  }

  # ---- 3. Pick the optimum subset size ---------------------------------------
  if (sensitive && nboot_eff > 1L) {
    if (is.null(min_size)) {
      # MATLAB: `nMin = ceil(0.1 * nChs)`. We follow that on the GOOD-channel
      # count, but floor at 2 because subset size 1 is never evaluated (row
      # 1 of `zmin_mean` is structurally NA), so a 1-floor would short-
      # circuit the sensitive search to `n_optimum = 1`.
      n_min <- max(2L, as.integer(ceiling(0.1 * n_good)))
    } else {
      n_min <- max(2L, as.integer(min_size))
    }
    zmm_x_trs <- rowMeans(zmin_mean, na.rm = TRUE)              # length n_good
    ii <- n_min
    n_optimum <- n_good
    while (ii <= n_good) {
      if (ii == n_good) {
        n_optimum <- ii
        break
      }
      if (!is.na(zmm_x_trs[ii + 1L]) &&
          zmm_x_trs[ii + 1L] > zmm_x_trs[ii]) {
        ii <- ii + 1L
        next
      }
      if (ii > 1L) {
        zmm_x_trs[seq_len(ii - 1L)] <- NA_real_
      }
      next_greater <- which(!is.na(zmm_x_trs) & zmm_x_trs > zmm_x_trs[ii])[1L]
      if (is.na(next_greater)) {
        n_optimum <- ii
        break
      }
      zmm_curr   <- zmin_mean[ii, ]
      idx        <- which.min(zmm_x_trs[seq_len(next_greater - 1L)])
      zmm_trough <- zmin_mean[idx, ]
      diffs <- as.numeric(outer(zmm_trough, zmm_curr, "-"))
      # left-tailed 95% CI on the bootstrap difference (trough - peak)
      conf_upper <- stats::quantile(diffs, 0.95, names = FALSE, type = 7L,
                                    na.rm = TRUE)
      if (!is.na(conf_upper) && conf_upper < 0) {
        n_optimum <- ii
        break
      }
      ii <- next_greater
    }
    if (n_optimum == n_min) {
      warning("CARLA optimum was detected at the minimum subset size floor.")
    }
  } else {
    if (n_tr == 1L) {
      zmm <- as.numeric(zmin_mean)
    } else {
      zmm <- rowMeans(zmin_mean, na.rm = TRUE)
    }
    n_optimum <- which.max(zmm)
  }

  channels <- sort(ord[seq_len(n_optimum)])
  # car <- apply(x[channels, , , drop = FALSE], c(2L, 3L), agg_fun)
  if (n_tr == 1L) {
    # car <- as.numeric(car)
    zmin_mean <- as.numeric(zmin_mean)
  }

  list(
    channels        = channels,
    # car             = car,
    order           = ord,
    vars            = vars,
    n_optimum       = as.integer(n_optimum),
    zmin_mean       = zmin_mean,
    bad_channels    = as.integer(bad_channels),
    virtual_channel = as.integer(virtual_channel),
    vars1           = vars1
  )
}
