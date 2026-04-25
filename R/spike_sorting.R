# detect: The detection method ('neg', 'pos', 'both'). Default is 'neg'.
# segment_duration: Duration of each segment in seconds. Default is 5 minutes.
# stdmin: Threshold factor for detection. Default is 5.
# stdmax: Threshold factor for artifact removal. Default is 50.
# w_pre: Pre-event window size.
# w_post: Post-event window size.
# min_ref_per: Minimum refractory period in milliseconds. Default is 1.5.

lillie.test <- function(x) {
  DNAME <- deparse(substitute(x))
  x <- sort(x[stats::complete.cases(x)])
  n <- length(x)
  if (n < 5)
    stop("sample size must be greater than 4")
  p <- stats::pnorm((x - mean(x)) / sd(x))
  Dplus <- max(seq(1:n) / n - p)
  Dminus <- max(p - (seq(1:n) - 1) / n)
  K <- max(Dplus, Dminus)
  if (n <= 100) {
    Kd <- K
    nd <- n
  } else {
    Kd <- K * ((n / 100)^0.49)
    nd <- 100
  }
  pvalue <- exp(-7.01256 * Kd^2 * (nd + 2.78019) + 2.99587 *
                  Kd * sqrt(nd + 2.78019) - 0.122119 + 0.974598 / sqrt(nd) +
                  1.67997 / nd)
  if (pvalue > 0.1) {
    KK <- (sqrt(n) - 0.01 + 0.85 / sqrt(n)) * K
    if (KK <= 0.302) {
      pvalue <- 1
    } else if (KK <= 0.5) {
      pvalue <- 2.76773 - 19.828315 * KK + 80.709644 *
        KK^2 - 138.55152 * KK^3 + 81.218052 * KK^4
    } else if (KK <= 0.9) {
      pvalue <- -4.901232 + 40.662806 * KK - 97.490286 *
        KK^2 + 94.029866 * KK^3 - 32.355711 * KK^4
    } else if (KK <= 1.31) {
      pvalue <- 6.198765 - 19.558097 * KK + 23.186922 *
        KK^2 - 12.234627 * KK^3 + 2.423045 * KK^4
    } else {
      pvalue <- 0
    }
  }
  RVAL <- list(statistic = c(D = K), p.value = pvalue, method = "Lilliefors (Kolmogorov-Smirnov) normality test",
               data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}

spike_detection <- function(
    x,
    sample_rate,
    x_bp2 = NULL,
    x_bp4 = NULL,
    detect_method = c("negative", "positive", "both"),
    segment_duration = 5,
    std_min = 5,
    std_max = 50,
    w_pre = 20,
    w_post = 44,
    min_refractory = 1.5
) {

  detect_method <- match.arg(detect_method)

  spike_signal <- x
  nyquist <- sample_rate / 2

  if (is.null(x_bp2)) {
    ellip2 <- ravetools::ellip(
      n = 2,
      Rp = 0.1,
      Rs = 40,
      w = c(300, 3000) / nyquist,
      type = "pass"
    )
    signal_bp2 <- ravetools::filtfilt(b = ellip2$b, a = ellip2$a, x = spike_signal)
  } else {
    signal_bp2 <- x_bp2
  }

  if (is.null(x_bp4)) {
    ellip4 <- ravetools::ellip(
      n = 4,
      Rp = 0.1,
      Rs = 40,
      w = c(300, 3000) / nyquist,
      type = "pass"
    )
    signal_bp4 <- ravetools::filtfilt(b = ellip4$b, a = ellip4$a, x = spike_signal)
  } else {
    signal_bp4 <- x_bp4
  }

  ref <- ceiling(min_refractory * sample_rate / 1000)
  sample_ref <- floor(ref / 2)

  # We process data by segments
  total_n_timepoints <- length(spike_signal)
  segment_n_timepoints <- sample_rate * segment_duration * 60
  n_segments <- ceiling(total_n_timepoints / segment_n_timepoints)
  pre_to_post <- seq(-w_pre, w_post)


  spike_timepoints <- lapply(seq_len(n_segments), function(ii_seg) {
    timepoint_start <- max((ii_seg - 1) *  segment_n_timepoints + 1 - w_pre, 1)
    timepoint_end <- min(ii_seg * segment_n_timepoints + w_pre, total_n_timepoints)
    time_sel <- seq(timepoint_start, timepoint_end)

    slice_orig <- spike_signal[time_sel]
    slice_bp2 <- signal_bp2[time_sel]
    slice_bp4 <- signal_bp4[time_sel]

    threshold <- std_min * median(abs(slice_bp4)) / 0.6745
    threhsold_max <- threshold * std_min

    switch(
      detect_method,
      "negative" = {
        slice_threshold <- -slice_bp4
      },
      "positive" = {
        slice_threshold <- slice_bp4
        xaux <- which(slice_bp4 > threshold)
      },
      {
        # both
        slice_threshold <- abs(slice_bp4)
      }
    )

    xaux <- which(slice_threshold > threshold)
    xaux <- xaux[xaux > w_pre & xaux < (length(time_sel) - w_post)]

    tmp_env <- new.env(parent = emptyenv())
    tmp_env$xaux0 <- 0
    tmp_env$index <- NULL

    lapply(seq_along(xaux), function(ii) {

      xaux_ii <- xaux[[ii]]

      if (xaux_ii >= tmp_env$xaux0 + ref) {
        iaux <- which.min(slice_bp2[(xaux_ii - 1) + seq_len(sample_ref)])

        if ( max(abs(slice_bp2[pre_to_post + xaux_ii])) < threhsold_max ) {
          tmp_env$xaux0 <- xaux_ii + iaux - 1
          tmp_env$index <- c(tmp_env$index, tmp_env$xaux0)
        }
      }
      return()
    })


    spike_timepoints <- time_sel[tmp_env$index]

    spike_timepoints
  })

  spike_timepoints <- unlist(spike_timepoints)

  unique(spike_timepoints)
}



extract_waveforms_for_channel <- function(
    spike_timepoints, x_bp2, w_pre = 20, w_post = 44, int_factor = 5,
    detect_method = c("negative", "positive", "both")) {

  detect_method <- match.arg(detect_method)

  signal_bp2 <- x_bp2

  total_spikes <- length(spike_timepoints)

  len_waveform <- w_pre + w_post

  # Range to extract for each spike (with padding for interpolation)
  extra <- 2
  window_len <- len_waveform + 2 * extra
  s <- seq_len(window_len) - 1
  ints <- seq(0, window_len - 1, by = 1 / int_factor)

  # time-index offset
  pre_to_post <- seq(-w_pre - extra, w_post + extra - 1)

  # subset and super-sample
  interpolated_spikes <- sapply(seq_len(total_spikes), function(ii) {

    idx <- pre_to_post + spike_timepoints[[ii]]
    if (ii == 1) {
      idx[idx < 1] <- 1
    } else if (ii == total_spikes) {
      idx[idx > total_spikes] <- total_spikes
    }
    spike <- signal_bp2[idx]

    # Spline interpolation (like scipy's splrep/splev)
    spline_fit <- stats::splinefun(s, spike, method = "natural")
    spline_fit(ints)
  }, simplify = TRUE, USE.NAMES = FALSE)
  interpolated_spikes <- t(interpolated_spikes)

  # Find alignment within window
  align_win <- seq((w_pre + extra - 1) * int_factor + 1, (w_pre + extra + 1) * int_factor + 1)

  switch(
    detect_method,
    "positive" = {
      iaux <- apply(interpolated_spikes[, align_win, drop = FALSE], 1, which.max)
    },
    "negative" = {
      iaux <- apply(interpolated_spikes[, align_win, drop = FALSE], 1, which.min)
    },
    {
      iaux <- apply(abs(interpolated_spikes[, align_win, drop = FALSE]), 1, which.max)
    }
  )

  iaux <- iaux + (w_pre + extra - 1) * int_factor
  pre_to_post_interp <- seq(- w_pre * int_factor + int_factor, w_post * int_factor,
                            by = int_factor)

  spike_waveforms <- sapply(seq_along(iaux), function(ii) {
    interpolated_spikes[ii, iaux[[ii]] + pre_to_post_interp]
  }, simplify = TRUE, USE.NAMES = FALSE)
  spike_waveforms <- t(spike_waveforms)

  # validate
  # matplot(t(spike_waveforms[sample(total_spikes, 1000, ), ]),
  #         type = "l", lty = 1, col = "gray80")
  spike_waveforms
}


haar_feature_extraction <- function(
    spike_waveforms,
    level = 4,
    ls = 64,
    max_inputs = ceiling(0.75 * 64),
    min_inputs = 10,
    nd = 10
) {
  # spikes: matrix or list of spike waveforms, each row or element is a waveform

  nspk <- nrow(spike_waveforms)
  spike_mat <- spike_waveforms

  cc <- matrix(0, nrow = nspk, ncol = ls)

  cc <- t(apply(spike_mat, 1L, function(spike) {
    # Haar wavelet decomposition
    w <- waveslim::dwt(spike, wf = "haar", n.levels = level, boundary = "periodic")
    coeffs <- c(unlist(w[paste0("d", seq_len(level))]), w[[paste0("s", level)]])
    unname(coeffs[seq_len(ls)])
  }))

  cc[is.na(cc)] <- 0
  # ls <- ncol(cc)

  ks <- sapply(seq_len(ls), function(ii) {

    coef_ii <- cc[, ii]
    m <- mean(coef_ii)
    thr_dist <- stats::sd(coef_ii) * 3
    thr_dist_min <- m - thr_dist
    thr_dist_max <- m + thr_dist
    aux <- cc[coef_ii > thr_dist_min & coef_ii < thr_dist_max, ii]

    stats <- 0
    if (length(aux) > 10) {
      stats <- lillie.test(aux)$statistic
    }

    unname(stats)

  })

  sorted_indices <- order(ks)
  # ind <- sorted_indices
  A <- ks[sorted_indices][(length(ks) - max_inputs + 1):length(ks)]
  ncoeff <- length(A)
  maxA <- max(A)

  if (nd > 1 && ncoeff >= nd) {
    d <- ((A[nd:length(A)] - A[1:(length(A) - nd + 1)]) / maxA) * (ncoeff / nd)
    all_above1 <- which(d >= 1)
  } else {
    d <- numeric(0)
    all_above1 <- integer(0)
  }

  if (length(all_above1) >= 2) {
    aux2 <- diff(all_above1)
    # temp_bla <- stats::filter(c(0, aux2, 0), rep(1/3, 3), sides = 2, circular = TRUE)
    temp_bla <- ravetools::convolve_signal(c(0, aux2, 0), rep(1 / 3, 3))
    temp_bla[1] <- aux2[1]
    temp_bla[length(temp_bla)] <- aux2[length(aux2)]
    thr_knee_diff <- all_above1[which(temp_bla[-1] == 1)[1]] + (nd / 2)
    inputs <- max_inputs - thr_knee_diff
  } else {
    inputs <- min_inputs
  }

  if (inputs > max_inputs) {
    inputs <- max_inputs
  } else if (inputs < min_inputs) {
    inputs <- min_inputs
  }

  coeff <- sorted_indices[(length(sorted_indices) - inputs + 1):length(sorted_indices)]

  return(cc[, coeff])
}

