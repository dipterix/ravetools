# # Design filters

#' Design a digital filter
#' @description
#' Provides 'FIR' and 'IIR' filter options; default is 'FIR', see also
#' \code{\link{design_filter_fir}}; for 'IIR' filters, see
#' \code{\link{design_filter_iir}}.
#'
#' @param sample_rate data sample rate
#' @param data data to be filtered, can be optional (\code{NULL})
#' @param method filter method, options are \code{"fir"} (default),
#' \code{"butter"}, \code{"cheby1"}, \code{"cheby2"}, and \code{"ellip"}
#' @param high_pass_freq,low_pass_freq high-pass or low-pass frequency,
#' see \code{\link{design_filter_fir}} or \code{\link{design_filter_iir}}
#' @param high_pass_trans_freq,low_pass_trans_freq transition bandwidths,
#' see \code{\link{design_filter_fir}} or \code{\link{design_filter_iir}}
#' @param passband_ripple allowable pass-band ripple in decibel; default is
#' \code{0.1}
#' @param stopband_attenuation minimum stop-band attenuation (in decibel) at
#' transition frequency; default is \code{40} dB.
#' @param filter_order suggested filter order; 'RAVE' may or may not adopt this
#' suggestion depending on the data and numerical feasibility
#' @param data_size used by 'FIR' filter design to determine maximum order,
#' ignored in 'IIR' filters; automatically derived from \code{data}
#' @param ... passed to filter generator functions
#' @returns If \code{data} is specified and non-empty, this function returns
#' filtered data via forward and backward \code{filtfilt}; if \code{data} is
#' \code{NULL}, then returns the generator function.
#'
#' @examples
#'
#'
#' sample_rate <- 200
#' t <- seq(0, 10, by = 1 / sample_rate)
#' x <- sin(t * 4 * pi) + sin(t * 20 * pi) +
#'   2 * sin(t * 120 * pi) + rnorm(length(t), sd = 0.4)
#'
#' # ---- Using FIR ------------------------------------------------
#'
#' # Low-pass filter
#' y1 <- design_filter(
#'   data = x,
#'   sample_rate = sample_rate,
#'   low_pass_freq = 3, low_pass_trans_freq = 0.5
#' )
#'
#' # Band-pass cheby1 filter 8-12 Hz with custom transition
#' y2 <- design_filter(
#'   data = x,
#'   method = "cheby1",
#'   sample_rate = sample_rate,
#'   low_pass_freq = 12, low_pass_trans_freq = .25,
#'   high_pass_freq = 8, high_pass_trans_freq = .25
#' )
#'
#' y3 <- design_filter(
#'   data = x,
#'   sample_rate = sample_rate,
#'   low_pass_freq = 80,
#'   high_pass_freq = 30
#' )
#'
#' oldpar <- par(mfrow = c(2, 1),
#'               mar = c(3.1, 2.1, 3.1, 0.1))
#' plot(t, x, type = 'l', xlab = "Time", ylab = "",
#'      main = "Mixture of 2, 10, and 60Hz", xlim = c(0,1))
#' # lines(t, y, col = 'red')
#' lines(t, y3, col = 'green')
#' lines(t, y2, col = 'blue')
#' lines(t, y1, col = 'red')
#' legend(
#'   "topleft", c("Input", "Low: 3Hz", "Pass 8-12Hz", "Pass 30-80Hz"),
#'   col = c(par("fg"), "red", "blue", "green"), lty = 1,
#'   cex = 0.6
#' )
#'
#' # plot pwelch
#' pwelch(x, fs = sample_rate, window = sample_rate * 2,
#'        noverlap = sample_rate, plot = 1, ylim = c(-100, 10))
#' pwelch(y1, fs = sample_rate, window = sample_rate * 2,
#'        noverlap = sample_rate, plot = 2, col = "red")
#' pwelch(y2, fs = sample_rate, window = sample_rate * 2,
#'        noverlap = sample_rate, plot = 2, col = "blue")
#' pwelch(y3, fs = sample_rate, window = sample_rate * 2,
#'        noverlap = sample_rate, plot = 2, col = "green")
#'
#'
#' # ---- Clean this demo --------------------------------------------------
#' par(oldpar)
#'
#'
#' @export
design_filter <- function(
    sample_rate,
    data = NULL,
    method = c("fir_kaiser", "firls", "fir_remez", "butter", "cheby1", "cheby2", "ellip"),
    high_pass_freq = NA, high_pass_trans_freq = NA,
    low_pass_freq = NA, low_pass_trans_freq = NA,
    passband_ripple = 0.1, stopband_attenuation = 40,
    filter_order = NA,
    ..., data_size = length(data)
) {
  method <- match.arg(method)
  if(data_size <= 0) {
    data_size <- NA
  }
  if(length(data) && is.na(data_size)) {
    data_size <- floor(length(data) - 1) / 3
  }

  if( startsWith(method, "fir") ) {
    method <- list(
      fir_kaiser = "kaiser",
      firls = "firls",
      fir_remez = "remez"
    )[[ method ]]
    filter <- design_filter_fir(
      method = method,
      sample_rate = sample_rate,
      filter_order = filter_order,
      data_size = data_size,
      high_pass_freq = high_pass_freq,
      high_pass_trans_freq = high_pass_trans_freq,
      low_pass_freq = low_pass_freq,
      low_pass_trans_freq = low_pass_trans_freq,
      stopband_attenuation = stopband_attenuation,
      ...
    )
  } else {
    filter <- design_filter_iir(
      method = method,
      sample_rate = sample_rate,
      filter_order = filter_order,
      high_pass_freq = high_pass_freq,
      high_pass_trans_freq = high_pass_trans_freq,
      low_pass_freq = low_pass_freq,
      low_pass_trans_freq = low_pass_trans_freq,
      passband_ripple = passband_ripple,
      stopband_attenuation = stopband_attenuation
    )
  }

  if(length(data)) {
    re <- filtfilt(b = filter$b, a = filter$a, x = data)
  } else {
    re <- filter
  }

  re

}

#' @title Design an 'IIR' filter
#' @param sample_rate sampling frequency
#' @param method filter method name, choices are \code{"butter"},
#' \code{"cheby1"}, \code{"cheby2"}, and \code{"ellip"}
#' @param filter_order suggested filter order. Notice filters with higher orders
#' may become numerically unstable, hence this number is only a suggested
#' number. If the filter is unstable, this function will choose a lower order;
#' leave this input \code{NA} (default) if undecided.
#' @param high_pass_freq high-pass frequency; default is \code{NA} (no high-pass
#' filter will be applied)
#' @param high_pass_trans_freq high-pass frequency band-width; default
#' is automatically inferred from filter type.
#' @param low_pass_freq low-pass frequency; default is \code{NA} (no
#' low-pass filter will be applied)
#' @param low_pass_trans_freq low-pass frequency band-width; default
#' is automatically inferred from filter type.
#' @param passband_ripple allowable pass-band ripple in decibel; default is
#' \code{0.1}
#' @param stopband_attenuation minimum stop-band attenuation (in decibel) at
#' transition frequency; default is \code{40} dB.
#' @returns A filter in 'Arma' form.
#'
#' @examples
#'
#' sample_rate <- 500
#'
#' my_diagnose <- function(
#'     filter, vlines = c(8, 12), cutoffs = c(-3, -6)) {
#'   diagnose_filter(
#'     b = filter$b,
#'     a = filter$a,
#'     fs = sample_rate,
#'     vlines = vlines,
#'     cutoffs = cutoffs
#'   )
#' }
#'
#' # ---- Default using butterworth to generate 8-12 bandpass filter ----
#'
#' # Butterworth filter with cut-off frequency
#' # 7 ~ 13 (default transition bandwidth is 1Hz) at -3 dB
#' filter <- design_filter_iir(
#'   method = "butter",
#'   low_pass_freq = 12,
#'   high_pass_freq = 8,
#'   sample_rate = 500
#' )
#'
#' filter
#'
#' my_diagnose(filter)
#'
#' ## explicit bandwidths and attenuation (sharper transition)
#'
#' # Butterworth filter with cut-off frequency
#' # passband ripple is 0.5 dB (8-12 Hz)
#' # stopband attenuation is 40 dB (5-18 Hz)
#' filter <- design_filter_iir(
#'   method = "butter",
#'   low_pass_freq = 12, low_pass_trans_freq = 6,
#'   high_pass_freq = 8, high_pass_trans_freq = 3,
#'   sample_rate = 500,
#'   passband_ripple = 0.5,
#'   stopband_attenuation = 40
#' )
#'
#' filter
#'
#' my_diagnose(filter)
#'
#' # ---- cheby1 --------------------------------
#'
#' filter <- design_filter_iir(
#'   method = "cheby1",
#'   low_pass_freq = 12,
#'   high_pass_freq = 8,
#'   sample_rate = 500
#' )
#'
#' my_diagnose(filter)
#'
#' # ---- cheby2 --------------------------------
#'
#' filter <- design_filter_iir(
#'   method = "cheby2",
#'   low_pass_freq = 12,
#'   high_pass_freq = 8,
#'   sample_rate = 500
#' )
#'
#' my_diagnose(filter)
#'
#' # ----- ellip ---------------------------------
#'
#' filter <- design_filter_iir(
#'   method = "ellip",
#'   low_pass_freq = 12,
#'   high_pass_freq = 8,
#'   sample_rate = 500
#' )
#'
#' my_diagnose(filter)
#'
#'
#'
#'
#' @export
design_filter_iir <- function(
    method = c("butter", "cheby1", "cheby2", "ellip"),
    sample_rate, filter_order = NA,
    high_pass_freq = NA, high_pass_trans_freq = NA,
    low_pass_freq = NA, low_pass_trans_freq = NA,
    passband_ripple = 0.1, stopband_attenuation = 40
) {
  method <- match.arg(method)

  # DIPSAUS DEBUG START
  # sample_rate <- 500
  # list2env(list(sample_rate = sample_rate,
  #                 high_pass_freq = 1, high_pass_trans_freq = 0.45,
  #                 low_pass_freq = 50, low_pass_trans_freq = 10), envir=.GlobalEnv)
  # list2env(list(filter_order = NA, low_pass_freq = NA, high_pass_freq = 4,
  #               low_pass_trans_freq = NA, high_pass_trans_freq = NA), envir=.GlobalEnv)
  # passband_ripple <- 0.1
  # stopband_attenuation <- 40
  # method <- "butter"


  nyquist <- sample_rate / 2
  passband_ripple <- abs(passband_ripple)
  stopband_attenuation <- abs(stopband_attenuation)

  r_pass <- passband_ripple
  r_stop <- stopband_attenuation

  stopifnot2(
    !is.na(low_pass_freq) || !is.na(high_pass_freq),
    msg = "Please specify low-pass and/or high-pass frequencies"
  )
  stopifnot2(is.finite(passband_ripple) && passband_ripple > 0,
             msg = "`passband_ripple` must be a finite positive number")
  stopifnot2(is.finite(stopband_attenuation) && stopband_attenuation > 0,
             msg = "`stopband_attenuation` must be a finite positive number")

  # suggested bandwidth if bandwidths are NA

  if(is.na(low_pass_trans_freq)) {
    low_pass_trans_freq <- min(max(0.25 * low_pass_freq, 2), low_pass_freq)
  }
  if(is.na(high_pass_trans_freq)) {
    high_pass_trans_freq <- min(max(0.25 * high_pass_freq, 2), sample_rate / 2 - high_pass_freq)
  }

  # determine the filter type
  if( is.na(high_pass_freq) ) {
    ftype <- "low"
    w_pass <- low_pass_freq / nyquist
    w_trans <- low_pass_trans_freq / nyquist

    w_stop <- clamp(w_pass + abs(w_trans), min = 0, max = 1)
  } else if ( is.na(low_pass_freq) ) {
    ftype <- "high"
    w_pass <- high_pass_freq / nyquist
    w_trans <- high_pass_trans_freq / nyquist

    w_stop <- clamp(w_pass - abs(w_trans), min = 0, max = 1)
  } else if ( isTRUE(high_pass_freq <= low_pass_freq) ) {
    ftype <- "pass"
    w_pass <- c(high_pass_freq, low_pass_freq) / nyquist
    w_trans <- c(high_pass_trans_freq, low_pass_trans_freq) / nyquist

    w_stop <- clamp( w_pass + c(-1, 1) * abs(w_trans) , min = 0, max = 1)
  } else {
    ftype <- "stop"
    w_pass <- c(low_pass_freq, high_pass_freq) / nyquist
    w_trans <- c(low_pass_trans_freq, high_pass_trans_freq) / nyquist

    w_stop <- clamp( w_pass + c(1, -1) * abs(w_trans) , min = 0, max = 1)
  }

  w_stop[w_stop < 0] <- 0
  w_stop[w_stop > 1] <- 1

  # get filter specs
  get_specs <- function(fun_ord, fun_make, min_order = 1) {
    spec <- fun_ord(Wp = w_pass, Ws = w_stop, Rp = r_pass, Rs = r_stop)
    validator <- function(n) {
      spec$n <- n
      filter <- gsignal::as.Arma(fun_make(spec))
      reciprocal_condition <- rcond_filter_ar(a = filter$a)
      reciprocal_condition > .Machine$double.eps
    }
    max_order <- guess_max_integer(validator, initial = spec$n, min_v = min_order)
    if(is.na(max_order)) {
      max_order <- 1
    }
    if( spec$n > max_order ) {
      spec$n <- max_order
    }
    spec
  }
  switch(
    method,
    "cheby1" = {
      fun_ord <- gsignal::cheb1ord
      fun_make <- gsignal::cheby1
      spec <- get_specs(fun_ord, fun_make)
      filter <- fun_make(spec)
    },
    "cheby2" = {
      fun_ord <- gsignal::cheb2ord
      fun_make <- gsignal::cheby2
      spec <- get_specs(fun_ord, fun_make, min_order = length(w_stop) * 2)
      filter <- fun_make(spec)
    },
    "ellip" = {
      fun_ord <- gsignal::ellipord
      fun_make <- gsignal::ellip
      spec <- get_specs(fun_ord, fun_make, min_order = 1)
      filter <- fun_make(spec)
    },
    {
      # "butter"
      spec <- gsignal::buttord(Wp = w_pass, Ws = w_stop, Rp = r_pass, Rs = r_stop)
      order <- butter_max_order(w = w_pass, r = r_pass, type = ftype)$n
      if(spec$n < order) {
        order <- spec$n
      }
      w_cutoff <- butter_cutoff(type = ftype, n = order, w = w_pass, r = r_pass)
      spec <- gsignal::FilterSpecs(n = order, Wc = w_cutoff, type = ftype)
      filter <- gsignal::butter(spec)
    }
  )

  checks <- check_filter(
    b = filter$b,
    a = filter$a,
    w = c(w_pass, w_stop),
    r_expected = -rep(c(r_pass, r_stop), each = length(w_pass)),
    fs = sample_rate
  )

  params <- list(
    type = c("iir", ftype),
    method = method,
    order = spec$n,
    passband = w_pass,
    stopband = w_stop,
    cutoff = spec$Wc,
    sample_rate = sample_rate
  )


  filter$parameters <- params
  filter$checks <- checks


  class(filter) <- c(
    "ravetools-design_filter_iir",
    sprintf("ravetools-design_filter_%s", method),
    "ravetools-design_filter",
    "ravetools-printable",
    class(filter)
  )

  filter
}


