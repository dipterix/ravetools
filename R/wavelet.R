#' @title 'Morlet' wavelet transform (Discrete)
#' @name wavelet
#' @description Transform analog voltage signals with 'Morlet'
#' wavelets: complex wavelet kernels with \eqn{\pi/2} phase
#' differences.
#' @param data numerical vector such as analog voltage signals
#' @param freqs frequency in which \code{data} will be projected on
#' @param srate sample rate, number of time points per second
#' @param wave_num desired number of cycles in wavelet kernels to
#' balance the precision in time and amplitude (control the
#' smoothness); positive integers are strongly suggested
#' @param precision the precision of computation; choices are
#' \code{'float'} (default) and \code{'double'}.
#' @param trend choices are \code{'constant'}: center the signal at zero;
#' \code{'linear'}: remove the linear trend; \code{'none'} do nothing
#' @param ... further passed to \code{\link{detrend}};
#' @returns \code{wavelet_kernels} returns wavelet kernels to be
#' used for wavelet function; \code{morlet_wavelet} returns a file-based array
#' if \code{precision} is \code{'float'}, or a list of real and imaginary
#' arrays if \code{precision} is \code{'double'}
#' @param frequency_range frequency range to calculate, default is 2 to 200
#' @param cycle_range number of cycles corresponding to \code{frequency_range}.
#' For default frequency range (2 - 200), the default \code{cycle_range} is
#' 3 to 20. That is, 3 wavelet kernel cycles at 2 Hertz, and 20 cycles at 200
#' Hertz.
#' @param signature signature to calculate kernel path to save, internally used
#' @param segment_length optional positive integer; when provided, long signals
#' are processed in overlapping segments of this length (in samples) using
#' batched \code{\link{mvfftw_c2c}} convolutions instead of a single full-length
#' FFT. This dramatically reduces peak memory and FFT cost for long recordings
#' (e.g. multi-hour). Must be strictly greater than the longest wavelet kernel
#' length (the kernel at the lowest frequency); otherwise an error is raised.
#' Default is \code{NULL}, which uses the legacy single-shot path. When
#' \code{segment_length >= length(data)} the function silently falls back to
#' the legacy path. Results match the legacy path on the interior of the
#' signal up to floating-point error; the first/last
#' \code{ceiling(max_kernel_len/2)} samples may differ near the global
#' boundaries (both implementations have boundary artifacts there).
#'
#' @examples
#'
#' \donttest{
#'
#' # generate sine waves
#' time <- seq(0, 3, by = 0.01)
#' x <- sin(time * 20*pi) + exp(-time^2) * cos(time * 10*pi)
#'
#' plot(time, x, type = 'l')
#'
#' # freq from 1 - 15 Hz; wavelet using float precision
#' freq <- seq(1, 15, 0.2)
#' coef <- morlet_wavelet(x, freq, 100, c(2,3))
#'
#' # to get coefficients in complex number from 1-10 time points
#' coef[1:10, ]
#'
#' # power
#' power <- Mod(coef[])^2
#'
#' # Power peaks at 5Hz and 10Hz at early stages
#' # After 1.0 second, 5Hz component fade away
#' image(power, x = time, y = freq, ylab = "frequency")
#'
#' # wavelet using double precision
#' coef2 <- morlet_wavelet(x, freq, 100, c(2,3), precision = "double")
#' power2 <- (coef2$real[])^2 + (coef2$imag[])^2
#'
#' image(power2, x = time, y = freq, ylab = "frequency")
#'
#' # The maximum relative change of power with different precisions
#' max(abs(power/power2 - 1))
#'
#' # display kernels
#' freq <- seq(1, 15, 1)
#' kern <- wavelet_kernels(freq, 100, c(2,3))
#' print(kern)
#'
#' plot(kern)
#'
#' }
#'
NULL

#' @rdname wavelet
#' @export
wavelet_kernels <- function(freqs, srate, wave_num) {
  # calculate wavelet cycles for each frequencies
  if (length(wave_num) != length(freqs)) {
    # calculate wavelet cycles for each frequencies
    ratio <- (log(max(wave_num)) - log(min(wave_num))) / (log(max(freqs)) - log(min(freqs)))
    wavelet_cycles <- round(exp((log(freqs) - log(min(freqs))) * ratio + log(min(wave_num))))
  } else {
    wavelet_cycles <- wave_num
  }
  f_l <- length(freqs)


  # wavelet window calc - each columns of final wave is a wavelet kernel (after fft)
  # sts = wavelet_cycles / (2 * pi * freqs)
  # wavelet_wins = cbind(-3 * sts, 3 * sts)

  fft_waves <- lapply(seq_len(f_l), function(ii) {
    fq <- freqs[ii]
    cycles <- wavelet_cycles[ii]
    # standard error
    st <- cycles / (2 * pi * fq)

    # calculate window size
    wavelet_win <- seq(-3 * st, 3 * st, by = 1 / srate)

    # half of window length
    w_l_half <- (length(wavelet_win) - 1) / 2

    # wavelet 1: calc sinus in complex domain
    tmp_sine <- exp((0 + 1i) * 2 * pi * fq / srate * (-w_l_half:w_l_half))

    # Gaussian normalization part
    A <- 1 / sqrt(st * sqrt(pi))

    # wavelet 2: calc gaussian wrappers
    tmp_gaus_win <- A * exp(-wavelet_win^2 / (2 * (cycles / (2 * pi * fq))^2))

    # wave kernel
    tmp_wavelet <- tmp_sine * tmp_gaus_win

    tmp_wavelet
  })

  structure(list(
    kernels = fft_waves,
    wavelet_cycles = wavelet_cycles,
    sample_rate = srate,
    frequencies = freqs
  ), class = "ravetools-wavelet-kernels")
}

#' @export
`print.ravetools-wavelet-kernels` <- function(x, plot = FALSE, ...) {
  cat("Discrete wavelet kernels\n")
  cat("  number of kernels/frequencies:", length(x$kernels), "\n")
  cat(sprintf("  frequency range: %.2f Hz - %.2f Hz\n", min(x$frequencies), max(x$frequencies)))
  cat(sprintf("  number of cycles: %.2f - %.2f\n", min(x$wavelet_cycles), max(x$wavelet_cycles)))
  invisible(x)
}

#' @export
`plot.ravetools-wavelet-kernels` <- function(
  x, cex = 1.2, cex.lab = cex * 1.2, cex.main = cex * 1.33,
  cex.axis = cex, mai = c(0.8, 0.5, 0.4, 0.1), ...) {

  fft_waves <- x$kernels
  srate <- x$sample_rate
  freqs <- x$frequencies
  wavelet_cycles <- x$wavelet_cycles
  max_l <- as.integer(max(sapply(fft_waves, length)) + 0.1 * srate)
  s <- sapply(fft_waves, function(s) {
    l <- (max_l - length(s))
    pre <- floor(l / 2)
    post <- ceiling(l / 2)
    c(rep(NA, pre), s, rep(NA, post))
  })
  s_re <- Re(s)
  s_im <- Im(s)
  ind <- exp(seq(log(min(freqs)), log(max(freqs)), length.out = 10))
  ind <- sapply(ind, function(i) {
    which.min(abs(i - freqs))[[1]]
  })
  ind <- unique(sort(ind))
  gap <- seq_along(ind) * 1.8 * max(abs(s_re), na.rm = TRUE)
  tmp_re <- t(t(s_re[, ind]) + gap)
  tmp_im <- t(t(s_im[, ind]) + gap)
  tmp <- rbind(tmp_re)
  x_all <- seq_len(max_l) / srate
  x_re <- x_all
  x_im <- x_re + max(x_all)
  # grid::grid.newpage()
  lay <- rbind(c(1, 1), c(2, 3))

  graphics::layout(mat = lay)

  old_mai <- graphics::par("mai")
  graphics::par(mai = mai)
  on.exit({
    graphics::par(mai = old_mai)
  }, add = TRUE)

  # old_mar <- par('mar')
  # on.exit({
  #   par(mar = old_mar)
  # })
  # par(mar = c(5.1, 4.1, 4.1, 2.1) * (cex / 2 + 0.5))
  graphics::matplot(y = tmp_re, x = x_re, type = "l", col = "red",
                    xlim = c(0, max(x_im)), ylim = c(min(tmp_re, na.rm = TRUE), max(gap) + 1.5 * min(gap)),
                    lty = 1, cex.lab = cex.lab, cex.main = cex.main, xlab = "Wavelet Length (seconds)", cex.axis = cex.axis,
                    ylab = "Frequency (Hz)", main = "Wavelet Kernels (Real & Imaginary)", yaxt = "n", xaxt = "n")

  graphics::matlines(y = tmp_im, x = x_im, type = "l", col = "red", lty = 1)

  n_halftickers <- 7
  x_actual <- c(x_re, x_im)
  x_label <- c(x_all, x_all) - mean(x_all)
  xind <- seq(1, length(x_re), length.out = n_halftickers)
  xind <- c(xind, xind + length(x_re))
  xind <- as.integer(xind[-n_halftickers])
  x_label <- x_label[xind]
  x_label[n_halftickers] <- abs(x_label[n_halftickers])
  x_label <- sprintf("%.2f", x_label)
  x_label[n_halftickers] <- paste0("\u00B1", x_label[n_halftickers])
  x_text <- stats::median(x_actual)

  graphics::axis(1, at = x_actual[xind], x_label, cex.axis = cex.axis)
  graphics::axis(2, at = gap, freqs[ind], cex.axis = cex.axis, las = 1)
  graphics::abline(h = gap, col = "grey80", lty = 2)
  leading_mod <- sapply(ind, function(ii) {
    x <- s[, ii]
    cycles <- wavelet_cycles[ii]
    x <- x[!is.na(x)] #Mod(x[1]) / max(Mod(x)) * 100  #= 1.111%
    c(length(x) / srate, cycles)
  })
  graphics::text(x = x_text, y = gap, "|", cex = cex)
  graphics::text(x = x_text, y = gap, sprintf("%.3f", leading_mod[1, ]), cex = cex,  pos = 2)
  graphics::text(x = x_text, y = gap, sprintf("%.2f", leading_mod[2, ]), cex = cex,  pos = 4)
  y_mini_title <- min(gap) + max(gap)
  graphics::text(x = x_text, y = y_mini_title, "|", cex = cex.lab)
  graphics::text(x = x_text, y = y_mini_title, "Wave Length", cex = cex.lab, pos = 2)
  graphics::text(x = x_text, y = y_mini_title, "# of Cycles", cex = cex.lab, pos = 4)

  # plot freq over wavelength and wave cycles
  wave_len <- sapply(fft_waves, length) / srate
  graphics::plot(freqs, wave_len, type = "l", ylab = "Wavelet Length (seconds)",
                 xlab = "Frequency (Hz)", main = "Wavelet Length | Frequency",
                 las = 1, cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis, col = "grey80")
  graphics::points(freqs, wave_len, col = "red", pch = 4)

  graphics::plot(freqs, wavelet_cycles, type = "l", ylab = "Wavelet Cycle",
                 xlab = "Frequency (Hz)", main = "Wavelet Cycle | Frequency",
                 las = 1, cex.lab = cex.lab, cex.main = cex.main, cex.axis = cex.axis, col = "grey80")
  graphics::points(freqs, wavelet_cycles, col = "red", pch = 4)

  return(invisible())
}


# Build a segmentation plan to chunk a long signal into overlapping segments.
#
# The overlap on each side equals ceil(max_kernel_len / 2) so the interior
# rows of each segment (i.e. samples at least `overlap` rows from both segment
# ends) are NOT corrupted by that segment's circular-FFT wraparound and
# therefore match the legacy single-shot result up to floating-point error.
wavelet_segment_plan <- function(data_length, segment_length, max_kernel_len) {
  data_length <- as.integer(data_length)
  segment_length <- as.integer(segment_length)
  max_kernel_len <- as.integer(max_kernel_len)
  if (segment_length <= max_kernel_len) {
    stop(sprintf(
      "`segment_length` (%d) must be greater than the longest wavelet kernel length (%d). Increase `segment_length`, raise the lowest frequency, or reduce the wavelet cycle count.",
      segment_length, max_kernel_len
    ))
  }
  overlap <- as.integer(ceiling(max_kernel_len / 2))
  step <- segment_length - 2L * overlap
  if (step < 1L) {
    stop(sprintf(
      "`segment_length` (%d) is too small relative to the wavelet kernel length (%d).",
      segment_length, max_kernel_len
    ))
  }
  n_seg <- as.integer(ceiling(data_length / step))
  out_starts <- (seq_len(n_seg) - 1L) * step + 1L
  out_ends <- pmin(out_starts + step - 1L, data_length)
  in_starts <- out_starts - overlap
  in_ends <- in_starts + segment_length - 1L
  local_starts <- rep(overlap + 1L, n_seg)
  local_ends <- overlap + (out_ends - out_starts + 1L)
  list(
    n_seg = n_seg,
    overlap = overlap,
    step = step,
    segment_length = segment_length,
    data_length = data_length,
    out_starts = out_starts,
    out_ends = out_ends,
    in_starts = in_starts,
    in_ends = in_ends,
    local_starts = local_starts,
    local_ends = local_ends
  )
}

# Build a segment_length x n_seg real matrix from a 1D signal using a plan.
# Out-of-range input rows are zero-padded.
wavelet_build_segment_matrix <- function(data, plan) {
  d_l <- length(data)

  if (storage.mode(data) != "double") {
    storage.mode(data) <- "double"
  }

  segment_length <- plan$segment_length
  seg_mat <- vapply(seq_len(plan$n_seg), function(k) {
    in_s <- plan$in_starts[k]
    in_e <- plan$in_ends[k]
    src_s <- max(in_s, 1L)
    src_e <- min(in_e, d_l)
    data_seg <- data[src_s:src_e]
    if (length(data_seg) == segment_length) {
      return(data_seg)
    } else {
      re <- double(segment_length)
      if (src_e >= src_s) {
        dst_s <- as.integer(src_s - in_s + 1L)
        dst_e <- as.integer(src_e - in_s + 1L)
        re[dst_s:dst_e] <- data_seg
      }
      return(re)
    }
  }, FUN.VALUE = double(segment_length))
  seg_mat
}

# Returns a normalized integer segment_length when segmentation should run;
# otherwise NA_integer_ (legacy single-shot path).
wavelet_normalize_segment_length <- function(segment_length, data_length) {
  if (is.null(segment_length)) { return(NA_integer_) }
  if (length(segment_length) != 1L) {
    stop("`segment_length` must be a single value or NULL")
  }
  if (is.na(segment_length)) { return(NA_integer_) }
  segment_length <- suppressWarnings(as.integer(segment_length))
  if (is.na(segment_length) || segment_length <= 0L) { return(NA_integer_) }
  if (segment_length >= data_length) { return(NA_integer_) }
  segment_length
}


## Precision-specific implementations live in:
##   R/wavelet-float.R   (morlet_wavelet_float, _segmented, kernels)
##   R/wavelet-double.R  (morlet_wavelet_double, _segmented, kernels)


#' @rdname wavelet
#' @export
morlet_wavelet <- function(data, freqs, srate, wave_num, precision = c("float", "double"),
                           trend = c("constant", "linear", "none"),
                           signature = NULL, segment_length = NULL, ...) {
  precision <- match.arg(precision)
  if (precision == "float") {
    re <- morlet_wavelet_float(data = data, freqs = freqs, srate = srate,
                               wave_num = wave_num, trend = trend,
                               signature = signature,
                               segment_length = segment_length, ...)
  } else {
    re <- morlet_wavelet_double(data = data, freqs = freqs, srate = srate,
                               wave_num = wave_num, trend = trend,
                               signature = signature,
                               segment_length = segment_length, ...)
  }
  return(re)
}

#' @rdname wavelet
#' @export
wavelet_cycles_suggest <- function(
    freqs,
    frequency_range = c(2, 200),
    cycle_range = c(3, 20)
) {
  v1 <- log(cycle_range[[1]])
  v2 <- log(cycle_range[[2]])
  cycle <- (v2 - v1) / (log(frequency_range[[2]]) - log(frequency_range[[1]])) *
    (log(freqs) - log(frequency_range[[1]])) + v1
  cycle <- round(exp(cycle))
  cycle[cycle <= 0] <- 1

  data.frame(
    Frequency = freqs,
    Cycles = cycle
  )
}


# x <- rnorm(10000)
# y2 <- morlet_wavelet(x, freqs = 2:200, srate = 2000, wave_num = c(2,20), demean = TRUE)
# y1 <- rave::wavelet(x, freqs = 2:200, srate = 2000, wave_num = c(2,20), demean = TRUE)

