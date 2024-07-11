
#' Design \code{'FIR'} filter using \code{\link{firls}}
#' @param sample_rate sampling frequency
#' @param filter_order filter order, leave \code{NA} (default) if undecided
#' @param data_size minimum length of data to apply the filter, used to
#' decide the maximum filter order. For 'FIR' filter, data length must be
#' greater than \code{3xfilter_order}
#' @param high_pass_freq high-pass frequency; default is \code{0} (no high-pass
#' filter will be applied)
#' @param high_pass_trans_freq high-pass frequency band-width; default
#' is automatically inferred from data size.
#' Frequency \code{high_pass_freq - high_pass_trans_freq} if the cutoff
#' frequency at power \code{-trans_attenuation}; see input
#' \code{trans_attenuation} or run \code{print(design_filter_fir)}
#' for illustration
#' @param low_pass_freq low-pass frequency; default is 'Nyquist' frequency (no
#' low-pass filter will be applied)
#' @param low_pass_trans_freq low-pass frequency band-width; default
#' is automatically inferred from data size.
#' Frequency \code{low_pass_freq + low_pass_trans_freq} if the cutoff
#' frequency at power \code{-trans_attenuation}; see input
#' \code{trans_attenuation} or run \code{print(design_filter_fir)} for
#' illustration
#' @param trans_attenuation allowable power attenuation (in decibel) at
#' transition frequency; default is \code{3} dB, which is half power.
#'
#' @details
#' Filter type is determined from \code{high_pass_freq} and
#' \code{low_pass_freq}. High-pass frequency is ignored if \code{high_pass_freq}
#' is less equal to 0, hence the filter is low-pass filter. When
#' \code{low_pass_freq} is greater equal to 'Nyquist' frequency, then
#' the filter is high-pass filter. When both \code{high_pass_freq} and
#' \code{low_pass_freq} are valid (positive, less than 'Nyquist'), then
#' the filter is a band-pass filter if band-pass is less than low-pass
#' frequency, otherwise the filter is band-stop.
#'
#' Although the peak amplitudes are set by \code{low_pass_freq} and
#' \code{high_pass_freq}, the transition from peak amplitude to zero require
#' a transition, which is tricky but also important to set.
#' When 'FIR' filters have too steep transition boundaries, the filter tends to
#' have ripples in peak amplitude, introducing artifacts to the final signals.
#' When the filter is too flat, components from unwanted frequencies may also
#' get aliased into the filtered signals. Ideally, the transition bandwidth
#' cannot be too steep nor too flat. In this function, users may control
#' the transition frequency bandwidths via \code{low_pass_trans_freq} and
#' \code{high_pass_trans_freq}. The power at the end of transition is defined
#' by \code{trans_attenuation}, with default value of \code{3} (i.e.
#' -3 dB, this number is automatically negated during the calculation).
#' By design, a low-pass 5 Hz filter with 1 Hz transition bandwidth results in
#' attenuating -3 dB power at 6 Hz, and -6 dB at 7 Hz (\code{5 + 1 x 2}).
#' Common choices of \code{trans_attenuation} are \code{3} (half power at
#' the transition boundary) and \code{6} (half amplitude at the transition
#' boundary).
#'
#'
#' @returns Function that generates the 'FIR' filter given filter order
#'
#' @examples
#'
#' # see illustration by `print(design_filter_fir)`
#'
#' # ---- FIR band-pass example -----------------------------
#'
#' sample_rate <- 500
#'
#' fir_generator <- design_filter_fir(
#'   sample_rate = sample_rate,
#'   high_pass_freq = 1, high_pass_trans_freq = 0.25,
#'   low_pass_freq = 50, low_pass_trans_freq = 6
#' )
#'
#' # check generator information
#' print(fir_generator)
#'
#' filter <- fir_generator()
#'
#' # -3dB cutoff (50% power) should be roughly within
#' # 0.75 (=1-0.25) ~ 56 (=50+6) Hz
#' diagnose_filter(
#'   filter$b, filter$a,
#'   fs = sample_rate,
#'   n = 1000,
#'   vlines = c(1, 50)
#' )
#'
#' # ---- High-pass filter with auto bandwith from data -----
#'
#' sample_rate <- 500
#'
#' fir_generator <- design_filter_fir(
#'   sample_rate = sample_rate,
#'   data_size = 1000,
#'   high_pass_freq = 6
#' )
#'
#' print(fir_generator)
#'
#' filter <- fir_generator()
#'
#' # -3dB cutoff should cover 6 Hz to Nyquist frequency
#' diagnose_filter(
#'   filter$b, filter$a,
#'   fs = sample_rate,
#'   n = 1000,
#'   vlines = 6,
#'   xlim = c(0, 50)
#' )
#'
#' # ---- Low-pass filter with auto bandwith from filter order -----
#'
#' sample_rate <- 500
#'
#' fir_generator <- design_filter_fir(
#'   sample_rate = sample_rate,
#'   filter_order = 200,
#'   low_pass_freq = 6
#' )
#'
#' print(fir_generator)
#'
#' filter <- fir_generator()
#'
#' # -3dB cutoff should cover 0-6 Hz
#' diagnose_filter(
#'   filter$b, filter$a,
#'   fs = sample_rate,
#'   n = 1000,
#'   vlines = 6,
#'   xlim = c(0, 50)
#' )
#'
#' # ---- Band-stop filter when data size is small -----
#' sample_rate <- 500
#'
#' # cannot have too sharp transition
#' fir_generator <- design_filter_fir(
#'   sample_rate = sample_rate,
#'   filter_order = 200,
#'   low_pass_freq = 6,
#'   high_pass_freq = 8
#' )
#'
#' print(fir_generator)
#'
#' filter <- fir_generator(scale = FALSE)
#'
#' # -3dB cutoff should cover 0-6 Hz & 8-Nyquist
#' diagnose_filter(
#'   filter$b, filter$a,
#'   fs = sample_rate,
#'   n = 1000,
#'   vlines = 6,
#'   xlim = c(0, 50)
#' )
#'
#'
#' @export
design_filter_fir <- function(
    sample_rate, filter_order = NA, data_size = NA,
    high_pass_freq = NA, high_pass_trans_freq = NA,
    low_pass_freq = NA, low_pass_trans_freq = NA,
    trans_attenuation = 3
) {

  # DIPSAUS DEBUG START
  # debug <- function(x, ...) { message(x, ...) }
  # sample_rate <- 500
  # data_size <- 1000
  # list2env(list(filter_order = NA, low_pass_freq = sample_rate/2, high_pass_freq = 4,
  #               low_pass_trans_freq = NA, high_pass_trans_freq = NA), envir=.GlobalEnv)
  # list2env(list(sample_rate = sample_rate,
  #                 high_pass_freq = 1, high_pass_trans_freq = 0.45,
  #                 low_pass_freq = 50, low_pass_trans_freq = 10), envir=.GlobalEnv)
  # trans_attenuation <- 3


  nyquist <- sample_rate / 2

  stopifnot2(
    !is.na(low_pass_freq) || !is.na(high_pass_freq),
    msg = "Please specify low-pass and/or high-pass frequencies"
  )

  if(is.na(data_size)) {
    max_nfilters <- Inf
  } else {
    max_nfilters <- max(floor((data_size - 1) / 3), 2)
  }
  if(!is.na(data_size)) {
    if(!is.na(filter_order)) {
      stopifnot2(
        filter_order <= max_nfilters,
        msg = "Filter order must be less than 1/3 of the data size"
      )
    }
  }

  r_peak <- 0
  r_pass <- abs(trans_attenuation) # amp = 0.5 if trans_attenuation=6
  r_stop <- Inf

  # calculate stop-bandwidth relative to transition bandwidth
  stopifnot2(is.finite(r_pass) && r_pass > 0, msg = "`trans_attenuation` must be finite positive number")
  factor <- 1 / (1 - 10^(-r_pass / 20))


  min_bandwidth <- 1 / 30
  if( !is.na(filter_order) ) {
    min_bandwidth <- 1 / (filter_order - 1)
  } else if(is.finite(max_nfilters)){
    min_bandwidth <- 1 / (max_nfilters - 1)
  }
  min_bandwidth <- min_bandwidth * nyquist / factor
  if(is.na(high_pass_trans_freq)) {
    high_pass_trans_freq <- min_bandwidth
  }
  if(is.na(low_pass_trans_freq)) {
    low_pass_trans_freq <- min_bandwidth
  }

  high_pass_trans_freq <- abs(high_pass_trans_freq)
  low_pass_trans_freq <- abs(low_pass_trans_freq)

  if( is.na(high_pass_freq) || high_pass_freq <= 0 ) {
    ftype <- "low"
    w_peak <- low_pass_freq / nyquist
    w_pass <- clamp((low_pass_freq + low_pass_trans_freq) / nyquist, min = 0, max = 1)
    w_stop <- clamp((low_pass_freq + low_pass_trans_freq * factor) / nyquist, min = 0, max = 1)
  } else if ( is.na(low_pass_freq) || low_pass_freq >= nyquist ) {
    ftype <- "high"
    w_peak <- high_pass_freq / nyquist
    w_pass <- clamp((high_pass_freq - high_pass_trans_freq) / nyquist, min = 0, max = 1)
    w_stop <- clamp((high_pass_freq - high_pass_trans_freq * factor) / nyquist, min = 0, max = 1)
  } else if ( isTRUE(high_pass_freq < low_pass_freq) ) {
    ftype <- "pass"
    w_peak <- c(high_pass_freq, low_pass_freq) / nyquist
    w_pass <- clamp(c(
      high_pass_freq - high_pass_trans_freq,
      low_pass_freq + low_pass_trans_freq
    ) / nyquist, min = 0, max = 1)
    w_stop <- clamp(c(
      high_pass_freq - high_pass_trans_freq * factor,
      low_pass_freq + low_pass_trans_freq * factor
    ) / nyquist, min = 0, max = 1)
  } else {
    ftype <- "stop"
    w_peak <- c(high_pass_freq, low_pass_freq) / nyquist
    w_pass <- clamp(c(
      high_pass_freq - high_pass_trans_freq,
      low_pass_freq + low_pass_trans_freq
    ) / nyquist, min = 0, max = 1)
    w_stop <- clamp(c(
      high_pass_freq - high_pass_trans_freq * factor,
      low_pass_freq + low_pass_trans_freq * factor
    ) / nyquist, min = 0, max = 1)
    if( w_stop[[1]] < w_stop[[2]] || w_pass[[1]] < w_pass[[2]] ) {
      if( high_pass_trans_freq + low_pass_trans_freq > 0 ) {
        w_stop_fac <- (high_pass_freq - low_pass_freq) / (high_pass_trans_freq + low_pass_trans_freq)
        if( w_stop_fac < 1 ) {
          w_pass <- clamp(c(
            max((high_pass_freq - high_pass_trans_freq * w_stop_fac) / nyquist, 0),
            min((low_pass_freq + low_pass_trans_freq * w_stop_fac) / nyquist, 1)
          ), min = 0, max = 1)
        }
        w_stop <- clamp(c(
          max((high_pass_freq - high_pass_trans_freq * w_stop_fac) / nyquist, 0),
          min((low_pass_freq + low_pass_trans_freq * w_stop_fac) / nyquist, 1)
        ), min = 0, max = 1)
      }
    }
  }

  params <- check_filter_params(
    type = ftype,
    w_peak = w_peak, w_pass = w_pass, w_stop = w_stop,
    r_peak = r_peak, r_pass = r_pass, r_stop = r_stop
  )
  amps0 <- 10^(-c(r_peak, r_pass, r_stop) / 20)

  # triange
  switch(
    ftype,
    "low" = {
      freqs <- c(0, w_peak, w_peak, w_pass, w_pass, w_stop, w_stop, 1)
      amps <- c(1, rep(amps0, each = 2), 0)
    },
    "high" = {
      freqs <- c(0, w_stop, w_stop, w_pass, w_pass, w_peak, w_peak, 1)
      amps <- rev(c(1, rep(amps0, each = 2), 0))
    },
    "pass" = {
      freqs <- rbind(w_stop, w_pass, w_peak)
      freqs <- unname(c(0, rep(freqs[,1], each = 2), rev(rep(freqs[,2], each = 2)), 1))
      amps <- c(0, rev(rep(amps0, each = 2)), rep(amps0, each = 2), 0)
    },
    "stop" = {
      freqs <- rbind(w_stop, w_pass, w_peak)
      freqs <- unname(c(0, rev(rep(freqs[,2], each = 2)), rep(freqs[,1], each = 2), 1))
      amps <- c(1, rep(amps0, each = 2), rev(rep(amps0, each = 2)), 1)
    }
  )

  # calculate minimal filter order
  slope <- diff(amps) / diff(freqs)
  slope <- abs(slope[is.finite(slope)])

  if(length(slope)) {
    min_nfilters <- ceiling(max(slope))
  } else {
    min_nfilters <- 2
  }

  if(is.na(filter_order)) {
    filter_order <- max_nfilters
    if( filter_order > min_nfilters ) {
      filter_order <- min_nfilters
    }
  }

  s1 <- switch(
    params$type,
    "pass" = c(
      "  The frequency response is (approximately) given by:",
      "",
      "     1-|               ----------",
      "       |             /|         | \\",
      "       |            / |         |  \\",
      "       |           /  |         |   \\",
      " |H|   |          /   |         |    \\",
      "   Rpa-|---------/----+---------+-----\\",
      "       |        /|    |         |     |\\           ",
      "       |       / |    |         |     | \\          ",
      "       |      /  |    |         |     |  \\         ",
      "     0-|-----/   |    |         |     |   \\--------|",
      "       0    Fs1 Fp1   Fh        Fl   Fp2  Fs2       Nyq",
      "",
      "  Where:",
      "",
      "      * Nyq:      (Nyquist) ",
      "            = `sample_rate` / 2",
      "      * Rpa:      (Transition power attenuation converted to amplitude) ",
      "            = 10 ^ { - `trans_attenuation` / 20 }",
      "      * Fh:       (High-pass frequency) ",
      "            = `high_pass_freq`",
      "      * Fl:       (Low-pass frequency) ",
      "            = `low_pass_freq`",
      "      * Fp1, Fp2: (Transition frequency)",
      "        Fp1 = Fh - `high_pass_trans_freq`",
      "        Fp2 = Fl + `low_pass_trans_freq`",
      "      * Fs1, Fs2: (Stop-band frequency, auto-derived)"
    ),
    "stop" = c(
      "  The frequency response is (approximately) given by:",
      "",
      "     1-|---------                  ----------|",
      "       |       | \\                / |       |",
      "   |H| |       |  \\              /  |       |",
      "   Rpa-|-------+--|\\------------/|  |       |",
      "       |       |  | \\          / |  |       |",
      "     0-|-------|  |  ----------  |  |       |",
      "       |       |  |  |        |  |  |       |",
      "       0      Fl Fp2 Fs2    Fs1 Fp1 Fh     Nyq",
      "",
      "  Where:",
      "",
      "      * Nyq:      (Nyquist) ",
      "            = `sample_rate` / 2",
      "      * Rpa:      (Transition power attenuation converted to amplitude) ",
      "            = 10 ^ { - `trans_attenuation` / 20 }",
      "      * Fh:       (High-pass frequency) ",
      "            = `high_pass_freq`",
      "      * Fl:       (Low-pass frequency) ",
      "            = `low_pass_freq`",
      "      * Fp1, Fp2: (Transition frequency)",
      "        Fp1 = Fh - `high_pass_trans_freq`",
      "        Fp2 = Fl + `low_pass_trans_freq`",
      "      * Fs1, Fs2: (Stop-band frequency, auto-derived)"
    ),
    "low" = c(
      "  The frequency response is (approximately) given by:",
      "",
      "     1-|---------              ",
      " |H|   |       | \\            ",
      "   Rpa-|-------+--\\           ",
      "       |       |  |\\          ",
      "       |       |  | \\         ",
      "     0-|-------|  |  ----------",
      "       |       |  |  |        |",
      "       0       Fl Fp Fs       Nyq",
      "",
      "  Where:",
      "",
      "      * Nyq:      (Nyquist) ",
      "            = `sample_rate` / 2",
      "      * Rpa:      (Transition power attenuation converted to amplitude) ",
      "            = 10 ^ { - `trans_attenuation` / 20 }",
      "      * Fl:       (Low-pass frequency) ",
      "            = `low_pass_freq`",
      "      * Fp:       (Transition frequency)",
      "            = Fl + `low_pass_trans_freq`",
      "      * Fs:       (Stop-band frequency, auto-derived)"
    ),
    "high" = c(
      "  The frequency response is (approximately) given by:",
      "",
      "     1-|               -----------",
      "       |             / |         |",
      "   Rpa-|------------/  |         |",
      " |H|   |           /|  |         |",
      "       |          / |  |         |",
      "     0-|----------  |  |         |",
      "       |         |  |  |         |",
      "       0        Fs  Fp Fh        Nyq",
      "",
      "  Where:",
      "",
      "      * Nyq:      (Nyquist) ",
      "            = `sample_rate` / 2",
      "      * Rpa:      (Transition power attenuation converted to amplitude) ",
      "            = 10 ^ { - `trans_attenuation` / 20 }",
      "      * Fh:       (High-pass frequency) ",
      "            = `high_pass_freq`",
      "      * Fp:       (Transition frequency)",
      "            = Fh - `high_pass_trans_freq`",
      "      * Fs:       (Stop-band frequency, auto-derived)"
    )
  )
  printable_message <- c(
    "<RAVE filter generator function>",
    sprintf("  Filter type : fir (%s)", params$type),
    sprintf("  Filter order: %.0f (suggested)", filter_order),
    sprintf("    Min order: %.0f (to achieve transition)", min_nfilters),
    sprintf("    Max order: %.0f (restricted by data size)", max_nfilters),
    "  Input parameters for `firls`: (not actual filter)",
    sprintf("    Peak frequency: [%s] x Nyquist (%s, -%.4g dB)",
            paste(sprintf("%.4g", params$w_peak), collapse = ","),
            paste(sprintf("%.1f Hz", params$w_peak * nyquist), collapse = "~"),
            params$r_peak),
    sprintf("    Transition frequency: [%s] x Nyquist (%s, -%.4g dB)",
            paste(sprintf("%.4g", params$w_pass), collapse = ","),
            paste(sprintf("%.1f Hz", params$w_pass * nyquist), collapse = "~"),
            params$r_pass),
    sprintf("    Stop-band frequency: [%s] x Nyquist (%s, -%.4g dB)",
            paste(sprintf("%.4g", params$w_stop), collapse = ","),
            paste(sprintf("%.1f Hz", params$w_stop * nyquist), collapse = "~"),
            params$r_stop),
    s1,
    "",
    "Generator Function Definition:",
    "  function(data, order = NULL, window = hamming, scale = TRUE, hilbert = FALSE) { ... }",
    "Arguments:",
    "  - `data`   : (optional) data to be filtered",
    "  - `order`  : (optional) Integer(1) filter order if not default",
    "  - `window` : (optional) Function(1) or `NULL`, function for the `FIR` window; default is `hamming`",
    "  - `scale`  : (optional) Logical(1) whether to scale the filter for unity gain; default is `TRUE`",
    "  - `hilbert`: (optional) Logical(1) use hilbert transform; default is `FALSE`.",
    "Values: ",
    "  Filtered data if `data` is provided (using `filtfilt`); otherwise will be the filter."
  )

  re <- structure(
    function(data, order = NULL, window = hamming, scale = TRUE, hilbert = FALSE) {
      if(!length(order) || !is.finite(order)) {
        order <- filter_order
      }
      has_data <- !missing(data)

      if( has_data ) {
        mo <- floor((length(data) - 1) / 3)
        if(order > mo) {
          order <- mo
        }
      }

      if( hilbert ) {
        re <- firls(N = order, freq = freqs, A = amps, ftype = "h")
      } else {
        re <- firls(N = order, freq = freqs, A = amps, ftype = "")
      }


      window_size <- length(re$b)

      if(is.function(window)) {
        wind <- window( window_size )
      } else if(length(window)) {
        wind <- window
        stopifnot2(
          length(wind) == window_size,
          msg = sprintf(
            "Expected filter window size: %d. However, the input window size is %d",
            window_size, length(wind))
        )
      } else {
        wind <- NULL
      }

      if(length(wind)) {
        b <- wind * re$b
        if(scale) {
          len <- length(b)

          switch(
            ftype,
            "low" = {
              f0 <- w_peak / 2
            },
            "high" = {
              f0 <- (1 + w_peak) / 2
            },
            "stop" = {
              # unity gain at DC
              f0 <- 0
            },
            {
              # unity gain at center of passband
              f0 <- mean(w_peak)
            }
          )
          b <- b / Mod(sum(
            exp(-1i * 2 * pi * seq(0, len - 1) * (f0 / 2)) * b
          ))

        }

        re$b <- b
      }

      if( has_data ) {
        re <- filtfilt(b = re$b, a = re$a, x = data)
      } else {
        re$check <- check_filter(
          b = re$b,
          a = re$a,
          w = c(w_peak, w_pass, w_stop),
          r_expected = rep(c(r_peak, r_pass, r_stop), each = length(w_peak))
        )
      }

      re
    },
    sample_rate = sample_rate,
    filter_type = c("fir", params$type),
    filter_order = list(suggested = filter_order, min = min_nfilters, max = max_nfilters),
    firls_args = data.frame(
      w = freqs,
      amp = amps
    ),
    checks = params,
    printable_message = printable_message,
    class = c("ravetools-design_filter_fir", "ravetools-design_filter",
              "ravetools-printable", "function")
  )

  re
}

attr(design_filter_fir, "printable_message") <- c(
  "Usage: ",
  "",
  deparse(as.call(c(list(quote(design_filter_fir)), formals(design_filter_fir)))),
  "",
  "Using band-pass filter as an example.",
  "  The frequency response is (approximately) given by:",
  "",
  "     1-|               ----------",
  "       |             /|         | \\",
  "       |            / |         |  \\",
  "       |           /  |         |   \\",
  " |H|   |          /   |         |    \\",
  "   Rpa-|---------/----+---------+-----\\",
  "       |        /|    |         |     |\\           ",
  "       |       / |    |         |     | \\          ",
  "       |      /  |    |         |     |  \\         ",
  "     0-|-----/   |    |         |     |   \\--------|",
  "       0    Fs1 Fp1   Fh        Fl   Fp2  Fs2       Nyq",
  "",
  "  Where:",
  "",
  "      * Nyq:      (Nyquist) ",
  "            = `sample_rate` / 2",
  "      * Rpa:      (Transition attenuation in amplitude) ",
  "            = 10 ^ { - `trans_attenuation` / 20 }",
  "      * Fh:       (High-pass frequency) ",
  "            = `high_pass_freq`",
  "      * Fl:       (Low-pass frequency) ",
  "            = `low_pass_freq`",
  "      * Fp1, Fp2: (Transition frequency)",
  "        Fp1 = Fh - `high_pass_trans_freq`",
  "        Fp2 = Fl + `low_pass_trans_freq`",
  "      * Fs1, Fs2: (Stop-band frequency, auto-derived)"
)

class(design_filter_fir) <- c("ravetools-printable", "function")
