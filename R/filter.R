# # Design filters

#' Design a digital filter
#' @description
#' Provides 'FIR' and 'IIR' filter options; default is 'FIR', see also
#' \code{\link{design_filter_fir}}; for 'IIR' filters, see
#' \code{\link{design_filter_iir}}.
#'
#' @param sample_rate data sample rate
#' @param data data to be filtered, can be optional (\code{NULL})
#' @param family filter family, options are \code{"fir"} (default),
#' \code{"butter"}, \code{"cheby1"}, \code{"cheby2"}, and \code{"ellip"}
#' @param high_pass_freq,low_pass_freq high-pass or low-pass frequency,
#' see \code{\link{design_filter_fir}} or \code{\link{design_filter_iir}}
#' @param high_pass_trans_freq,low_pass_trans_freq transition bandwidths,
#' see \code{\link{design_filter_fir}} or \code{\link{design_filter_iir}}
#' @param peak_attenuation desired magnitude (decibel) at low-pass high-pass
#' frequencies, see \code{\link{design_filter_iir}}; ignored in 'FIR' filters.
#' @param trans_attenuation desired magnitude (decibel) at edges of transition
#' frequencies; see \code{\link{design_filter_fir}} or
#' \code{\link{design_filter_iir}}
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
#' oldpar <- par(mfrow = c(2, 1),
#'               mar = c(3.1, 2.1, 3.1, 0.1))
#' # ---- Using FIR ------------------------------------------------
#'
#' data_size <- length(x)
#'
#' # Low-pass filter with no window, with transition bandwidth 0.5Hz
#' # with attenuation 3, i.e. expected -3dB power at 3.5Hz
#' f1 <- design_filter(
#'   data_size = data_size,
#'   sample_rate = sample_rate,
#'   low_pass_freq = 3, low_pass_trans_freq = 0.5,
#'   trans_attenuation = 3
#' )
#'
#' # Band-pass filter 8-12 Hz with custom transition
#' # of 1Hz, resulting in 3 dB decay around 7.75-12.25 Hz
#' f2 <- design_filter(
#'   data_size = data_size,
#'   sample_rate = sample_rate,
#'   low_pass_freq = 12, low_pass_trans_freq = .25,
#'   high_pass_freq = 8, high_pass_trans_freq = .25
#' )
#'
#'
#' f3 <- design_filter(
#'   data_size = data_size,
#'   sample_rate = sample_rate,
#'   low_pass_freq = 80, low_pass_trans_freq = 10,
#'   high_pass_freq = 30, high_pass_trans_freq = 5,
#'   trans_attenuation = 12
#' )
#'
#' # Check power attenuation
#' y1 <- f1(x)
#' y2 <- f2(x)
#' y3 <- f3(x)
#'
#'
#' plot(t, x, type = 'l', xlab = "Time", ylab = "",
#'      main = "Mixture of 0.2, 2, and 60Hz", xlim = c(0,1))
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
#'
#'
#' @export
design_filter <- function(
    sample_rate,
    data = NULL,
    family = c("fir", "butter", "cheby1", "cheby2", "ellip"),
    high_pass_freq = NA, high_pass_trans_freq = NA,
    low_pass_freq = NA, low_pass_trans_freq = NA,
    peak_attenuation = NA, trans_attenuation = 3,
    filter_order = NA,
    ..., data_size = length(data)
) {
  family <- match.arg(family)
  if(data_size <= 0) {
    data_size <- NA
  }
  if(length(data) && is.na(data_size)) {
    data_size <- floor(length(data) - 1) / 3
  }

  if( family == "fir" ) {
    generator <- design_filter_fir(
      sample_rate = sample_rate,
      filter_order = filter_order,
      data_size = data_size,
      high_pass_freq = high_pass_freq,
      high_pass_trans_freq = high_pass_trans_freq,
      low_pass_freq = low_pass_freq,
      low_pass_trans_freq = low_pass_trans_freq,
      trans_attenuation = trans_attenuation
    )
  } else {
    generator <- design_filter_iir(
      family = family,
      sample_rate = sample_rate,
      filter_order = filter_order,
      high_pass_freq = high_pass_freq,
      high_pass_trans_freq = high_pass_trans_freq,
      low_pass_freq = low_pass_freq,
      low_pass_trans_freq = low_pass_trans_freq,
      peak_attenuation = peak_attenuation,
      trans_attenuation = trans_attenuation
    )
  }

  if(length(data)) {
    re <- generator(..., data = data, order = NULL)
  } else {
    re <- generator
  }

  re

}

#' @title Design an 'IIR' filter
#' @param sample_rate sampling frequency
#' @param family filter family name, choices are \code{"butter"},
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
#' @param peak_attenuation expected power attenuation (in decibel) at high-pass
#' low pass frequencies. For 'Butterworth' filters, we recommend set to
#' \code{NA} and let the algorithm estimate automatically;
#' for all other filters, this is the allowable pass-band ripple, and
#' default value is \code{1} dB.
#' @param trans_attenuation power attenuation (in decibel) at
#' transition frequency. For 'Butterworth' filters, default is \code{3} dB
#' (half power); For all other filters, this is the minimum attenuation in the
#' stop band and default is \code{12} dB.
#' @returns A filter generator function.
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
#' iir_generator <- design_filter_iir(
#'   family = "butter",
#'   low_pass_freq = 12,
#'   high_pass_freq = 8,
#'   sample_rate = 500
#' )
#'
#' # magnitude decay is
#' #  8-12 Hz: around -0.2 dB
#' #  7-8 Hz and 12-13 Hz: -3 dB
#' filter <- iir_generator()
#' print(filter$check)
#'
#' my_diagnose(filter)
#'
#' ## explicit bandwidths and attenuation (sharper transition)
#'
#' # Butterworth filter with cut-off frequency
#' # 7.75 ~ 12.25 Hz, at -6 dB
#' iir_generator <- design_filter_iir(
#'   family = "butter",
#'   low_pass_freq = 12, low_pass_trans_freq = 0.25,
#'   high_pass_freq = 8, high_pass_trans_freq = 0.25,
#'   sample_rate = 500, trans_attenuation = 6
#' )
#'
#' # magnitude decay is
#' #  8-12 Hz: around -3 dB
#' #  7-8 Hz and 12-13 Hz: -6 dB
#' filter <- iir_generator()
#' filter$check
#'
#' my_diagnose(filter)
#'
#' # ---- cheby1 --------------------------------
#'
#' # pass-band: 8-12 Hz with ripple 1 dB
#' # stop-band: 7-13 Hz with minimum power attenuation 12 dB (default)
#' filter <- design_filter_iir(
#'   family = "cheby1",
#'   low_pass_freq = 12,
#'   high_pass_freq = 8,
#'   sample_rate = 500
#' )()
#'
#' my_diagnose(filter)
#'
#' # ---- cheby2 --------------------------------
#'
#' # pass-band: 8-12 Hz with ripple 1 dB
#' # stop-band: 7-13 Hz with minimum power attenuation 40 dB
#' filter <- design_filter_iir(
#'   family = "cheby2",
#'   low_pass_freq = 12,
#'   high_pass_freq = 8,
#'   trans_attenuation = 40,
#'   sample_rate = 500
#' )()
#'
#' my_diagnose(filter)
#'
#' # ----- ellip ---------------------------------
#'
#' # pass-band: 8-12 Hz with ripple 1 dB
#' # stop-band: 7.75-12.25 Hz with minimum power attenuation 30 dB
#' filter <- design_filter_iir(
#'   family = "ellip",
#'   low_pass_freq = 12, low_pass_trans_freq = 0.25,
#'   high_pass_freq = 8, high_pass_trans_freq = 0.25,
#'   trans_attenuation = 30,
#'   sample_rate = 500
#' )()
#'
#' my_diagnose(filter)
#'
#'
#'
#' @export
design_filter_iir <- function(
    family = c("butter", "cheby1", "cheby2", "ellip"),
    sample_rate, filter_order = NA,
    high_pass_freq = NA, high_pass_trans_freq = NA,
    low_pass_freq = NA, low_pass_trans_freq = NA,
    peak_attenuation = NA, trans_attenuation = NA
) {
  family <- match.arg(family)

  # DIPSAUS DEBUG START
  # sample_rate <- 500
  # data_size <- 1000
  # list2env(list(sample_rate = sample_rate,
  #                 high_pass_freq = 1, high_pass_trans_freq = 0.45,
  #                 low_pass_freq = 50, low_pass_trans_freq = 10), envir=.GlobalEnv)
  # list2env(list(filter_order = NA, low_pass_freq = sample_rate/2, high_pass_freq = 4,
  #               low_pass_trans_freq = NA, high_pass_trans_freq = NA), envir=.GlobalEnv)
  # trans_attenuation <- 3
  # peak_attenuation <- NA
  # family <- "butter"


  nyquist <- sample_rate / 2

  stopifnot2(
    !is.na(low_pass_freq) || !is.na(high_pass_freq),
    msg = "Please specify low-pass and/or high-pass frequencies"
  )

  # get transition
  if(!is.na(low_pass_freq)) {
    if(is.na(low_pass_trans_freq)) {
      if(family %in% c("cheby2")) {
        low_pass_trans_freq <- max(0.25 * low_pass_freq, 2.0)
      } else {
        low_pass_trans_freq <- min(0.25 * low_pass_freq, 1.0)
      }
    }
    if (low_pass_trans_freq > low_pass_freq) {
      low_pass_trans_freq <- low_pass_freq
    }
  }
  if(!is.na(high_pass_freq)) {
    if(is.na(high_pass_trans_freq)) {
      if(family %in% c("cheby2")) {
        high_pass_trans_freq <- max(0.25 * high_pass_freq, 2.0)
      } else {
        high_pass_trans_freq <- min(0.25 * high_pass_freq, 1.0)
      }
    }
    if (high_pass_trans_freq + high_pass_freq > nyquist) {
      high_pass_trans_freq <- nyquist - high_pass_freq
    }
  }


  if( is.na(high_pass_freq) || high_pass_freq <= 0 ) {
    ftype <- "low"
    w_peak <- low_pass_freq / nyquist
    w_pass <- clamp((low_pass_freq + low_pass_trans_freq) / nyquist, min = 0, max = 1)
  } else if ( is.na(low_pass_freq) || low_pass_freq >= nyquist ) {
    ftype <- "high"
    w_peak <- high_pass_freq / nyquist
    w_pass <- clamp((high_pass_freq - high_pass_trans_freq) / nyquist, min = 0, max = 1)
  } else if ( isTRUE(high_pass_freq < low_pass_freq) ) {
    ftype <- "pass"
    w_peak <- c(high_pass_freq, low_pass_freq) / nyquist
    w_pass <- clamp(c(
      high_pass_freq - high_pass_trans_freq,
      low_pass_freq + low_pass_trans_freq
    ) / nyquist, min = 0, max = 1)
  } else {
    ftype <- "stop"
    w_peak <- c(high_pass_freq, low_pass_freq) / nyquist
    w_pass <- clamp(c(
      high_pass_freq - high_pass_trans_freq,
      low_pass_freq + low_pass_trans_freq
    ) / nyquist, min = 0, max = 1)
    if( w_pass[[1]] < w_pass[[2]] ) {
      if( high_pass_trans_freq + low_pass_trans_freq > 0 ) {
        w_stop_fac <- (high_pass_freq - low_pass_freq) / (high_pass_trans_freq + low_pass_trans_freq) * 0.9
        w_pass <- clamp(c(
          max((high_pass_freq - high_pass_trans_freq * w_stop_fac) / nyquist, 0),
          min((low_pass_freq + low_pass_trans_freq * w_stop_fac) / nyquist, 1)
        ), min = 0, max = 1)
      }
    }
  }

  trans_attenuation <- abs(trans_attenuation)

  filter_max_order <- NULL
  make_filter <- NULL
  switch(
    family,
    "cheby1" = {
      if(is.na(peak_attenuation)) {
        r_peak <- 1
      } else {
        r_peak <- abs(peak_attenuation)
      }
      if(is.na(trans_attenuation)) {
        r_pass <- 12
      } else {
        r_pass <- abs(trans_attenuation)
      }
      filter_ord <- gsignal::cheb1ord
      make_filter_impl <- gsignal::cheby1
    },
    "cheby2" = {
      if(is.na(peak_attenuation)) {
        r_peak <- 1
      } else {
        r_peak <- abs(peak_attenuation)
      }
      if(is.na(trans_attenuation)) {
        r_pass <- 12
      } else {
        r_pass <- abs(trans_attenuation)
      }
      filter_ord <- gsignal::cheb2ord
      make_filter_impl <- gsignal::cheby2
    },
    "ellip" = {
      if(is.na(peak_attenuation)) {
        r_peak <- 1
      } else {
        r_peak <- abs(peak_attenuation)
      }
      if(is.na(trans_attenuation)) {
        r_pass <- 12
      } else {
        r_pass <- abs(trans_attenuation)
      }
      filter_ord <- gsignal::ellipord
      make_filter_impl <- gsignal::ellip
    },
    {
      # "butter"
      if(is.na(trans_attenuation)) {
        r_pass <- 3
      } else {
        r_pass <- abs(trans_attenuation)
      }
      if(is.na(peak_attenuation)) {
        r_peak <- NA
      } else {
        r_peak <- min(abs(peak_attenuation), r_pass / 2)
      }

      filter_ord <- gsignal::buttord
      make_filter_impl <- gsignal::butter

      make_filter <- function(order, ...) {
        Wc <- butter_cutoff(type = ftype, n = order, w = w_pass, r = r_pass)
        gsignal::butter(n = order, w = Wc, type = ftype, ...)
      }

      filter_max_order <- function() {
        butter_max_order(w = w_pass, type = ftype, r = r_pass)$n
      }
    }
  )

  if(!is.function(make_filter)) {
    make_filter <- function(order, ...) {
      spec <- filter_ord(Wp = sort(w_peak), Ws = sort(w_pass), Rp = r_peak, Rs = r_pass)
      spec$n <- order
      make_filter_impl(spec, ...)
    }
  }

  if( !is.function(filter_max_order) ) {
    filter_max_order <- function() {
      spec <- filter_ord(Wp = sort(w_peak), Ws = sort(w_pass), Rp = r_peak, Rs = r_pass)
      validator <- function(n) {
        spec$n <- n
        filter <- gsignal::as.Arma(make_filter_impl(spec))
        check <- check_filter(b = filter$b, a = filter$a)
        check$reciprocal_condition > .Machine$double.eps
      }
      n <- guess_max_integer(validator, min_v = spec$n)
      if(is.na(n)) { n <- 1 }
      n
    }
  }

  # calculate stop-bandwidth relative to transition bandwidth
  stopifnot2(is.na(r_peak) || r_peak > 0, msg = "`peak_attenuation` must be non-zero")
  stopifnot2(is.finite(r_pass) && r_pass > 0, msg = "`trans_attenuation` must be non-zero")


  params <- check_filter_params(
    type = ftype,
    w_peak = w_peak, w_pass = w_pass,
    r_peak = min(0, r_peak, na.rm = TRUE), r_pass = r_pass, r_stop = Inf
  )
  params$r_peak <- r_peak

  # get max_available filter given w_pass
  max_order <- filter_max_order()

  if( is.na(filter_order) ) {
    if( !is.na(r_peak) ) {
      # get n from transfer function
      spec <- filter_ord(Wp = sort(w_peak), Ws = sort(w_pass), Rp = r_peak, Rs = r_pass)
      filter_order <- spec$n
    } else {
      filter_order <- max_order
    }
  }

  filter_order <- min(filter_order, max_order)

  printable_message <- c(
    "<RAVE filter generator function>",
    sprintf("  Filter type : %s (%s)", family, params$type),
    sprintf("  Filter order: %.0f (suggested)", filter_order),
    sprintf("    Max order: %.0f (restricted by reciprocal condition number)",
            max_order),
    sprintf("  Input parameters for `%s`: (not actual filter)", family),
    sprintf("    Peak frequency: [%s] x Nyquist (%s, -%.4g dB)",
            paste(sprintf("%.4g", params$w_peak), collapse = ","),
            paste(sprintf("%.1f Hz", params$w_peak * nyquist), collapse = "~"),
            params$r_peak),
    sprintf("    Transition frequency: [%s] x Nyquist (%s, -%.4g dB)",
            paste(sprintf("%.4g", params$w_pass), collapse = ","),
            paste(sprintf("%.1f Hz", params$w_pass * nyquist), collapse = "~"),
            params$r_pass),
    "",
    "Generator Function Definition:",
    '  function(order, order = NULL, output = c("Arma", "Zpg", "Sos")) { ... }',
    "Arguments:",
    "  - `data`   : (optional) data to be filtered",
    "  - `order`  : (optional) Integer(1) filter order if not default",
    "  - `output` : (optional) String(1) filter form; default is ARMA model, only used if `data` is missing",
    "Values: ",
    "  Filtered data if `data` is provided (using `filtfilt`); otherwise will be the filter."
  )


  structure(
    function(data, order = NULL, output = c("Arma", "Zpg", "Sos")) {
      output <- match.arg(output)

      if(!length(order) || !is.finite(order)) {
        order <- filter_order
      } else if(order > max_order) {
        order <- max_order
      }

      filter <- make_filter(order = order, output = output)

      arma <- gsignal::as.Arma(filter)

      if(missing(data)) {
        filter$check <- check_filter(
          b = arma$b,
          a = arma$a,
          w = c(w_peak, w_pass),
          r_expected = rep(c(r_peak, r_pass), each = length(w_peak))
        )
        re <- filter
      } else {
        re <- filtfilt(b = arma$b, a = arma$a, x = data)
      }


      re

    },
    sample_rate = sample_rate,
    filter_type = c(family, params$type),
    filter_order = list(suggested = filter_order, min = 1, max = max_order),
    iir_args = data.frame(
      w = c(w_peak, w_pass),
      r = rep(c(r_peak, r_pass), each = length(w_peak))
    ),
    checks = params,
    printable_message = printable_message,
    class = c(sprintf("ravetools-design_filter_%s", ftype),
              "ravetools-design_filter_iir",
              "ravetools-design_filter",
              "ravetools-printable", "function")
  )
}


