
#' Design \code{'FIR'} filter using \code{\link{firls}}
#' @param method method to generate 'FIR' filter, default is using
#' \code{\link[gsignal]{kaiser}} estimate, other choices are
#' \code{\link{firls}} (with \code{hamming} window) and
#' \code{\link[gsignal]{remez}} design.
#' @param sample_rate sampling frequency
#' @param filter_order filter order, leave \code{NA} (default) if undecided
#' @param data_size minimum length of data to apply the filter, used to
#' decide the maximum filter order. For 'FIR' filter, data length must be
#' greater than \code{3xfilter_order}
#' @param high_pass_freq high-pass frequency; default is \code{NA} (no
#' high-pass filter will be applied)
#' @param high_pass_trans_freq high-pass frequency band-width; default
#' is automatically inferred from data size.
#' Frequency \code{high_pass_freq - high_pass_trans_freq} is the corner
#' of the stop-band
#' @param low_pass_freq low-pass frequency; default is \code{NA} (no
#' low-pass filter will be applied)
#' @param low_pass_trans_freq low-pass frequency band-width; default
#' is automatically inferred from data size.
#' Frequency \code{low_pass_freq + low_pass_trans_freq} is the corner
#' of the stop-band
#' @param stopband_attenuation allowable power attenuation (in decibel) at
#' transition frequency; default is \code{40} dB.
#' @param scale whether to scale the filter for unity gain
#'
#' @details
#' Filter type is determined from \code{high_pass_freq} and
#' \code{low_pass_freq}. High-pass frequency is ignored if \code{high_pass_freq}
#' is \code{NA}, hence the filter is low-pass filter. When
#' \code{low_pass_freq} is \code{NA}, then
#' the filter is high-pass filter. When both \code{high_pass_freq} and
#' \code{low_pass_freq} are valid (positive, less than 'Nyquist'), then
#' the filter is a band-pass filter if band-pass is less than low-pass
#' frequency, otherwise the filter is band-stop.
#'
#' Although the peak amplitudes are set at 1 by \code{low_pass_freq} and
#' \code{high_pass_freq}, the transition from peak amplitude to zero require
#' a transition, which is tricky but also important to set.
#' When 'FIR' filters have too steep transition boundaries, the filter tends to
#' have ripples in peak amplitude, introducing artifacts to the final signals.
#' When the filter is too flat, components from unwanted frequencies may also
#' get aliased into the filtered signals. Ideally, the transition bandwidth
#' cannot be too steep nor too flat. In this function, users may control
#' the transition frequency bandwidths via \code{low_pass_trans_freq} and
#' \code{high_pass_trans_freq}. The power at the end of transition is defined
#' by \code{stopband_attenuation}, with default value of \code{40} (i.e.
#' -40 dB, this number is automatically negated during the calculation).
#' By design, a low-pass 5 Hz filter with 1 Hz transition bandwidth results in
#' around -40 dB power at 6 Hz.
#'
#'
#' @returns 'FIR' filter in 'Arma' form.
#'
#' @examples
#'
#' # ---- Basic -----------------------------
#'
#' sample_rate <- 500
#' data_size <- 1000
#'
#' # low-pass at 5 Hz, with auto transition bandwidth
#' # from kaiser's method, with default stopband attenuation = 40 dB
#' filter <- design_filter_fir(
#'   low_pass_freq = 5,
#'   sample_rate = sample_rate,
#'   data_size = data_size
#' )
#'
#' # Passband ripple is around 0.08 dB
#' # stopband attenuation is around 40 dB
#' print(filter)
#'
#' diagnose_filter(
#'   filter$b, filter$a,
#'   fs = sample_rate,
#'   n = data_size,
#'   cutoffs = c(-3, -6, -40),
#'   vlines = 5
#' )
#'
#' # ---- Advanced ---------------------------------------------
#'
#' sample_rate <- 500
#' data_size <- 1000
#'
#' # Rejecting 3-8 Hz, with transition bandwidth 0.5 Hz at both ends
#' # Using least-square (firls) to generate FIR filter
#' # Suggesting the filter order n=160
#' filter <- design_filter_fir(
#'   low_pass_freq = 3, low_pass_trans_freq = 0.5,
#'   high_pass_freq = 8, high_pass_trans_freq = 0.5,
#'   filter_order = 160,
#'   sample_rate = sample_rate,
#'   data_size = data_size,
#'   method = "firls"
#' )
#'
#' #
#' print(filter)
#'
#' diagnose_filter(
#'   filter$b, filter$a,
#'   fs = sample_rate,
#'   n = data_size,
#'   cutoffs = c(-1, -40),
#'   vlines = c(3, 8)
#' )
#'
#'
#'
#' @export
design_filter_fir <- function(
    sample_rate, filter_order = NA, data_size = NA,
    high_pass_freq = NA, high_pass_trans_freq = NA,
    low_pass_freq = NA, low_pass_trans_freq = NA,
    stopband_attenuation = 40, scale = TRUE,
    method = c("kaiser", "firls", "remez")
) {

  method <- match.arg(method)

  # DIPSAUS DEBUG START
  # debug <- function(x, ...) { message(x, ...) }
  # sample_rate <- 200
  # data_size <- 2000
  # list2env(list(sample_rate = sample_rate,
  #                 high_pass_freq = NA, high_pass_trans_freq = NA,
  #                 low_pass_freq = 3, low_pass_trans_freq = NA), envir=.GlobalEnv)
  # list2env(list(filter_order = NA, low_pass_freq = 50, high_pass_freq = 4,
  #               low_pass_trans_freq = NA, high_pass_trans_freq = NA), envir=.GlobalEnv)
  # stopband_attenuation <- 40


  nyquist <- sample_rate / 2
  stopband_attenuation <- abs(stopband_attenuation)


  stopifnot2(
    !is.na(low_pass_freq) || !is.na(high_pass_freq),
    msg = "Please specify low-pass and/or high-pass frequencies"
  )
  stopifnot2(is.finite(stopband_attenuation) && stopband_attenuation > 0,
             msg = "`stopband_attenuation` must be a finite positive number")

  if(is.na(data_size)) {
    max_nfilters <- Inf
    # Huristic: 32 / pi / 2.285 / 333 (stopband atten 40 dB with data size 1000)
    suggested_trans_bandwidth <- 0.01
  } else {
    max_nfilters <- max(floor((data_size - 1) / 3), 2)
    if( max_nfilters %% 2 == 1 ) {
      max_nfilters <- max_nfilters - 1
    }

    if( isTRUE(filter_order > max_nfilters) ) {
      filter_order <- max_nfilters
    }
    # dw <- 2 * pi * diff(f)/fs
    # n <- max(1, ceiling((A - 8)/(2.285 * dw)))
    if( stopband_attenuation > 8 ) {
      suggested_trans_bandwidth <- (stopband_attenuation - 8) / pi / 2.285 / max_nfilters
    } else {
      suggested_trans_bandwidth <- stopband_attenuation / 22 / max_nfilters
    }
  }

  # determine the filter type
  if( is.na(high_pass_freq) ) {
    ftype <- "low"
    w_pass <- low_pass_freq / nyquist
    w_trans <- low_pass_trans_freq / nyquist
    w_trans[is.na(w_trans)] <- suggested_trans_bandwidth

    w_stop <- w_pass + abs(w_trans)

    # kaiserord
    bands <- c(w_pass, w_stop)
    mag <- c(1, 0)

    # firls
    mag2 <- c(1, 1, 0, 0)

  } else if ( is.na(low_pass_freq) ) {
    ftype <- "high"
    w_pass <- high_pass_freq / nyquist
    w_trans <- high_pass_trans_freq / nyquist
    w_trans[is.na(w_trans)] <- suggested_trans_bandwidth

    w_stop <- w_pass - abs(w_trans)

    # kaiserord
    bands <- c(w_stop, w_pass)
    mag <- c(0, 1)

    # firls
    mag2 <- c(0, 0, 1, 1)


  } else if ( isTRUE(high_pass_freq <= low_pass_freq) ) {
    ftype <- "pass"
    w_pass <- c(high_pass_freq, low_pass_freq) / nyquist
    w_trans <- c(high_pass_trans_freq, low_pass_trans_freq) / nyquist
    w_trans[is.na(w_trans)] <- suggested_trans_bandwidth

    w_stop <- clamp( w_pass + c(-1, 1) * abs(w_trans) , min = 0, max = 1)

    # kaiserord
    bands <- c(w_stop[[1]], w_pass, w_stop[[2]])
    mag <- c(0, 1, 0)

    # firls
    mag2 <- c(0, 0, 1, 1, 0, 0)

  } else {
    ftype <- "stop"
    w_pass <- c(low_pass_freq, high_pass_freq) / nyquist
    w_trans <- c(low_pass_trans_freq, high_pass_trans_freq) / nyquist
    w_trans[is.na(w_trans)] <- suggested_trans_bandwidth

    w_stop <- clamp( w_pass + c(1, -1) * abs(w_trans) , min = 0, max = 1)

    # kaiserord
    bands <- c(w_pass[[1]], w_stop, w_pass[[2]])
    mag <- c(1, 0, 1)

    # firls
    mag2 <- c(1, 1, 0, 0, 1, 1)
  }

  r_stop <- stopband_attenuation
  r_pass <- -20 * log10(1 - 10 ^ (-r_stop / 20))

  # get filter
  # c("kaiser", "ls", "remez")
  # Wc is (w_pass + w_stop) / 2
  kaisprm <- gsignal::kaiserord(bands, mag, dev = 10 ^ (-r_stop / 20))
  if(!is.na(filter_order)) {
    kaisprm$n <- filter_order
  }
  if(kaisprm$n > max_nfilters) {
    # even order for stop and high to ensure unity gain is non-zero at nyquist
    if (ftype %in% c("stop", "high") && max_nfilters %% 2 == 1) {
      max_nfilters <- max_nfilters + 1
    }
    kaisprm$n <- max_nfilters
  }

  params <- list(
    type = c("fir", ftype),
    method = method,
    order = kaisprm$n,
    scale = scale,
    passband = w_pass,
    stopband = w_stop,
    sample_rate = sample_rate
  )

  switch(
    method,
    "kaiser" = {
      filter <- fir1(
        n = kaisprm$n, w = kaisprm$Wc, type = kaisprm$type,
        window = gsignal::kaiser(kaisprm$n + 1, kaisprm$beta),
        scale = FALSE
      )
      params$cutoff <- kaisprm$Wc
      params$beta <- kaisprm$beta
      params
    },
    "firls" = {
      w <- c(0, rep(bands, each = 2), 1)
      a <- rep(mag2, each = 2)
      a <- a[-c(1, length(a))]
      filter <- firls(N = kaisprm$n, freq = w, A = a)
      wind <- hamming(length(filter$b))
      filter$b <- filter$b * wind
      params$frequency <- c(0, bands, 1)
      params$amplitude <- mag2
    },
    "remez" = {
      if(kaisprm$n < 4) {
        kaisprm$n <- 4
      }
      b <- unclass(gsignal::remez(kaisprm$n, c(0, bands, 1), mag2))
      filter <- structure(
        class = "Arma",
        list(b = b, a = 1)
      )
      params$frequency <- c(0, bands, 1)
      params$amplitude <- mag2
    }
  )
  if( scale ) {
    # scale
    switch(
      ftype,
      "low" = { f0 <- 0 },
      "high" = { f0 <- 1 },
      "pass" = {
        # unity gain at center of passband
        f0 <- mean(w_pass)
      },
      {
        if(w_stop[[1]] + w_stop[[2]] > 1) {
          # high-freqeuncy band stop
          f0 <- 0
        } else {
          f0 <- 1
        }
      }
    )
  }

  filter$parameters <- params
  filter$checks <- check_filter(b = filter$b, a = filter$a, w = c(w_pass, w_stop),
                                r_expected = -rep(c(r_pass, r_stop), each = length(w_pass)), fs = sample_rate)


  class(filter) <- c("ravetools-design_filter_fir", "ravetools-design_filter",
                     "ravetools-printable", class(filter))

  filter
}

#' @export
`format.ravetools-design_filter` <- function(x, ...) {
  params <- x$parameters
  checks <- sprintf("  %s", format(x$checks)[-c(1,2)])
  c(
    sprintf("<%s filter>", toupper(params$type[[1]])),
    sprintf("  Type : %s", params$type[[2]]),
    sprintf("  Method: %s", params$method),
    sprintf("  Order: %d", params$order),
    "  Magnitude:",
    checks
  )
}

# attr(design_filter_fir, "printable_message") <- c(
#   "Usage: ",
#   "",
#   deparse(as.call(c(list(quote(design_filter_fir)), formals(design_filter_fir)))),
#   "",
#   "Using band-pass filter as an example.",
#   "  The frequency response is (approximately) given by:",
#   "",
#   "     1-|               ----------",
#   "       |             /|         | \\",
#   "       |            / |         |  \\",
#   "       |           /  |         |   \\",
#   " |H|   |          /   |         |    \\",
#   "   Rpa-|---------/----+---------+-----\\",
#   "       |        /|    |         |     |\\           ",
#   "       |       / |    |         |     | \\          ",
#   "       |      /  |    |         |     |  \\         ",
#   "     0-|-----/   |    |         |     |   \\--------|",
#   "       0    Fs1 Fp1   Fh        Fl   Fp2  Fs2       Nyq",
#   "",
#   "  Where:",
#   "",
#   "      * Nyq:      (Nyquist) ",
#   "            = `sample_rate` / 2",
#   "      * Rpa:      (Transition attenuation in amplitude) ",
#   "            = 10 ^ { - `trans_attenuation` / 20 }",
#   "      * Fh:       (High-pass frequency) ",
#   "            = `high_pass_freq`",
#   "      * Fl:       (Low-pass frequency) ",
#   "            = `low_pass_freq`",
#   "      * Fp1, Fp2: (Transition frequency)",
#   "        Fp1 = Fh - `high_pass_trans_freq`",
#   "        Fp2 = Fl + `low_pass_trans_freq`",
#   "      * Fs1, Fs2: (Stop-band frequency, auto-derived)"
# )
#
# class(design_filter_fir) <- c("ravetools-printable", "function")
