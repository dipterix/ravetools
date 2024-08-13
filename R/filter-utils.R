# core util functions in signal filter

clamp <- function(x, min = NA, max = NA) {
  if(!is.na(min)) {
    x[x < min] <- min
  }
  if(!is.na(max)) {
    x[x > max] <- max
  }
  x
}

# ---- Shape/Window functions --------------------------------------------------
sineint <- function(x) {
  neg <- x < 0
  x[neg] <- -x[neg]

  re <- Im( pracma::expint( 1i * x ) ) + pi/2
  re[neg] <- -re[neg]
  re[x == 0] <- 0
  re
}

polyval <- function (coef, z)
{
  lz <- length(z)
  if (!lz)
    return(numeric(0))
  n <- length(coef)
  if (!n) {
    z[] <- 0
    return(z)
  }
  if (!(mode(coef) == "numeric") && !(mode(coef) == "complex"))
    stop("Argument 'coef' must be a real or complex vector.")
  d_z <- dim(z)
  dim(z) <- lz
  y <- outer(z, (n - 1):0, "^") %*% coef
  dim(y) <- d_z
  return(y)
}

#' @title Filter window functions
#' @name filter-window
#' @param n number of time-points in window
#' @returns A numeric vector of window with length \code{n}
#' @examples
#'
#' hanning(10)
#' hamming(11)
#' blackmanharris(21)
#'
NULL

#' @rdname filter-window
#' @export
hanning <- function (n) {
  if (!(n == round(n) && n > 0)) {
    stop("n has to be an integer > 0")
  }

  if (n == 1) {
    c <- 1
  } else {
    c <- 0.5 - 0.5 * cos(2 * pi * seq(0, n - 1)/(n - 1))
  }
  c
}

#' @rdname filter-window
#' @export
hamming <- function (n) {
  if (!(n == round(n) && n > 0)) {
    stop("n has to be an integer > 0")
  }

  if (n == 1) {
    c <- 1
  } else {
    n <- n - 1
    c <- 0.54 - 0.46 * cos(2 * pi * (0:n)/n)
  }
  c
}

#' @rdname filter-window
#' @export
blackman <- function(n) {
  if (!(n == round(n) && n > 0)) {
    stop("n has to be an integer > 0")
  }

  if (n == 1) {
    c <- 1
  } else {
    x <- seq.int(0, n - 1) * (2 * pi / (n - 1))
    c <- 0.42 - 0.5 * cos(x) + 0.08 * cos(x * 2)
  }
  c
}


#' @rdname filter-window
#' @export
blackmannuttall <- function(n) {
  if (!(n == round(n) && n > 0)) {
    stop("n has to be an integer > 0")
  }

  if (n == 1) {
    c <- 1
  } else {
    x <- seq.int(0, n - 1) * (2 * pi / (n - 1))
    c <- 0.3635819 - 0.4891775 * cos(x) + 0.1365995 * cos(x * 2) - 0.0106411 * cos(x * 3)
  }
  c
}

#' @rdname filter-window
#' @export
blackmanharris <- function(n) {
  if (!(n == round(n) && n > 0)) {
    stop("n has to be an integer > 0")
  }

  if (n == 1) {
    c <- 1
  } else {
    x <- seq.int(0, n - 1) * (2 * pi / (n - 1))
    c <- 0.35875 - 0.48829 * cos(x) + 0.14128 * cos(x * 2) - 0.01168 * cos(x * 3)
  }
  c
}

#' @rdname filter-window
#' @export
flattopwin <- function(n) {
  if (!(n == round(n) && n > 0)) {
    stop("n has to be an integer > 0")
  }

  if (n == 1) {
    c <- 1
  } else {
    x <- seq.int(0, n - 1) * (2 * pi / (n - 1))
    c <- 0.21557895 - 0.41663158 * cos(x) + 0.277263158 * cos(x * 2) - 0.083578947 * cos(x * 3) + 0.006947368 * cos(x * 4)
  }
  c
}

#' @rdname filter-window
#' @export
bohmanwin <- function(n) {
  if (!(n == round(n) && n > 0)) {
    stop("n has to be an integer > 0")
  }

  if (n == 1) {
    w <- 1
  } else {
    N <- n - 1
    k <- seq(-N / 2, N / 2)
    w <- (1 - 2 * abs(k) / N) * cos(2 * pi * abs(k) / N) + (1 / pi) *
      sin(2 * pi * abs(k) / N)
    w[1] <- w[length(w)] <- 0
  }
  w
}

# ---- fft transform -----------------------------------------------------------
fft <- function(x, inverse = FALSE){
  if(inverse){
    fftw_c2r(x)
  } else {
    fftw_r2c(x)
  }
}
ifft <- function(x){
  fftw_c2r(x) / length(x)
}



# ---- Estimation, validation,  ------------------------------------------------

filtmag_db <- function(b, a, f) {
  # Find filter's magnitude response in decibels at given frequency.
  nb <- length(b)
  na <- length(a)

  sapply(f, function(fi) {
    top <- sum(exp(-1i * seq(0, nb - 1) * pi * fi) * b)
    bot <- sum(exp(-1i * seq(0, na - 1) * pi * fi) * a)

    H <- 20 * log10(abs(top / bot))

    H
  })

}

#' Check filter parameters
#' @description
#' Internally (not exported) used to check the validity of filter inputs.
#' The function only checks if the inputs are ordered correctly, instead of
#' checking the feasibility of the filter.
#'
#' @param type filter type
#' @param w_peak frequency correspond to peak amplitude within 0 to 1, where 1
#' stands for 'Nyquist' frequency. This is the frequency range that users want
#' to "keep". For low-pass and high-pass filters \code{w_peak} is a single
#' number. For band-pass or band-stop filters, \code{w_peak} should be a
#' numeric vector of length 2, with first element being the high-pass frequency
#' and second element low-pass frequency. For band-pass filter, \code{w_peak}
#' is increasing; for band-stop filter, \code{w_peak} is decreasing.
#' @param r_peak expected decibel power (loss) of the peak power
#' @param w_pass pass-band frequency or frequency range, see to \code{w_peak}
#' for format
#' @param r_pass expected band-pass ripple
#' @param w_stop stop-band frequency or frequency range, see to \code{w_peak}
#' for format
#' @param r_stop expected stop-band attenuation
#'
#'
#' @noRd
check_filter_params <- function(
    type = c("low", "high", "pass", "stop"),
    w_peak, w_pass, w_stop,
    r_peak = 0.1, r_pass = 3, r_stop = 60
) {
  type <- match.arg(type)
  switch(
    type,
    "low" = {
      w_stop %?<-% 1
      stopifnot2(
        length(w_peak) == 1 && length(w_pass) == 1 && length(w_stop) == 1,
        msg = "Low-pass frequency cutoff length must be 1"
      )
      stopifnot2(
        0 < w_peak && w_peak <= w_pass && w_pass <= w_stop && w_stop <= 1,
        msg = "Low-pass frequency/nyquist must be within [0,1]. 0 < peak power frequency <= passband frequency <= stopband frequency <= 1"
      )
    },
    "high" = {
      w_stop %?<-% 0
      stopifnot2(
        length(w_peak) == 1 && length(w_pass) == 1 && length(w_stop) == 1,
        msg = "High-pass frequency cutoff length must be 1"
      )
      stopifnot2(
        0 <= w_stop && w_stop <= w_pass && w_pass <= w_peak && w_peak < 1,
        msg = "High-pass frequency/nyquist must be within [0,1]. 0 < stopband frequency <= passband frequency <= peak power frequency <= 1"
      )
    },
    "pass" = {
      w_stop %?<-% c(0, 1)
      stopifnot2(
        length(w_peak) == 2 && length(w_pass) == 2 && length(w_stop) == 2,
        msg = "Band-pass frequency cutoff length must be 2"
      )
      stopifnot2(
        0 <= w_stop[[1]] && w_stop[[1]] <= w_pass[[1]] && w_pass[[1]] <= w_peak[[1]] &&
          w_peak[[1]] <= w_peak[[2]] && w_peak[[2]] <= w_pass[[2]] && w_pass[[2]] <= w_stop[[2]] &&
          w_stop[[2]] <= 1,
        msg = "Band-pass frequency/nyquist must be within [0,1] in non-decreasing order: 0, w_stop[1], w_pass[1], w_peak[1], w_peak[2], w_pass[2], w_stop[2], 1"
      )
    },
    "stop" = {
      w_stop %?<-% c(1, 0)
      stopifnot2(
        length(w_peak) == 2 && length(w_pass) == 2 && length(w_stop) == 2,
        msg = "Band-stop frequency cutoff length must be 2"
      )
      stopifnot2(
        1 >= w_peak[[1]] && w_peak[[1]] >= w_pass[[1]] && w_pass[[1]] >= w_stop[[1]] &&
          w_stop[[1]] >= w_stop[[2]] && w_stop[[2]] >= w_pass[[2]] && w_pass[[2]] >= w_peak[[2]] &&
          w_peak[[2]] >= 0,
        msg = "Band-stop frequency/nyquist must be within [0,1] in non-decreasing order: 0, w_peak[2], w_pass[2], w_stop[2], w_stop[1], w_pass[1], w_peak[1], 1"
      )
    }
  )

  r_peak <- abs(r_peak)
  r_pass <- abs(r_pass)
  r_stop <- abs(r_stop)

  stopifnot2(
    0 <= r_peak && r_peak < r_pass && r_pass < r_stop,
    msg = "Absolute ripple (dB) must be in increasing order: 0 <= r_peak < r_pass < r_stop"
  )

  list(
    type = type,
    w_peak = w_peak,
    w_pass = w_pass,
    w_stop = w_stop,
    r_peak = r_peak,
    r_pass = r_pass,
    r_stop = r_stop
  )
}


#' Computer reciprocal condition number of an 'Arma' filter
#' @description
#' Test whether the filter is numerically stable for \code{\link{filtfilt}}.
#'
#' @param a auto-regression coefficient, numerical vector; the first element
#' must not be zero
#' @returns Reciprocal condition number of matrix \code{z1}, used in
#' \code{\link{filtfilt}}. If the number is less than
#' \code{.Machine$double.eps}, then \code{\link{filtfilt}} will fail.
#' @examples
#'
#'
#' # Butterworth filter with low-pass at 0.1 Hz (order = 4)
#' filter <- butter(4, 0.1, "low")
#'
#' # TRUE
#' rcond_filter_ar(filter$a) > .Machine$double.eps
#'
#' diagnose_filter(filter$b, filter$a, 500)
#'
#' # Bad filter (order is too high)
#' filter <- butter(50, 0.1, "low")
#'
#' rcond_filter_ar(filter$a) > .Machine$double.eps
#'
#' # filtfilt needs to inverse a singular matrix
#' diagnose_filter(filter$b, filter$a, 500)
#'
#' @seealso \code{\link{check_filter}}
#' @export
rcond_filter_ar <- function(a) {
  z1_r <- length(a)
  if(z1_r == 1) { return(1) }

  if(a[[1]] == 0) { stop("ARMA filter AR[1] must not be zero") }
  a <- a / a[[1]]

  if(z1_r == 2) { return( 1 + a[-1] ) }
  z1 <- diag(1, z1_r - 1) - cbind( -a[-1], rbind(diag(1, z1_r - 2), 0) )
  rcon <- abs(rcond(z1))
  rcon
}

#' @title Check 'Arma' filter
#' @param b moving average (\code{MA}) polynomial coefficients.
#' @param a auto-regressive (\code{AR}) polynomial coefficients.
#' @param w normalized frequency, ranging from 0 to 1, where 1 is 'Nyquist'
#' @param r_expected attenuation in decibel of each \code{w}
#' @param fs sample rate, used to infer the frequencies and formatting print
#' message, not used in calculation; leave it blank by default
#' @returns A list of power estimation and the reciprocal condition number
#' of the \code{AR} coefficients.
#' @examples
#'
#'
#' # create a butterworth filter with -3dB (half-power) at [1, 5] Hz
#' # and -60dB stop-band attenuation at [0.5, 6] Hz
#'
#' sample_rate <- 20
#' nyquist <- sample_rate / 2
#'
#' specs <- buttord(
#'   Wp = c(1, 5) / nyquist,
#'   Ws = c(0.5, 6) / nyquist,
#'   Rp = 3,
#'   Rs = 60
#' )
#' filter <- butter(specs)
#'
#' # filter quality is poor because the AR-coefficients
#' # creates singular matrix with unstable inverse,
#' # this will cause `filtfilt` to fail
#' check_filter(
#'   b = filter$b, a = filter$a,
#'
#'   # frequencies (normalized) where power is evaluated
#'   w = c(1, 5, 0.5, 6) / nyquist,
#'
#'   # expected power
#'   r_expected = c(3, 3, 60, 60)
#'
#' )
#'
#'
#' @export
check_filter <- function(b, a, w = NULL, r_expected = NULL, fs = NULL) {

  # r_expected <- abs(r_expected)

  if(length(w)) {
    attenuation <- data.frame(
      frequency = w,
      expected_db = r_expected,
      filtmag_db = filtmag_db(b, a, w)
    )
    attenuation$diff <- attenuation$filtmag_db - attenuation$expected_db
  } else {
    attenuation <- data.frame()
  }

  # atten_passed <- abs(attenuation$diff) < tol
  # if(!all(atten_passed)) {
  #   warn <- c(warn, " * Some frequencies might have unexpected attenuation")
  # }

  rcon <- rcond_filter_ar(a)

  structure(
    list(
      b = b,
      a = a,
      fs = fs,
      attenuation = attenuation,
      reciprocal_condition = rcon
    ),
    class = c("ravetools-test_filter", "ravetools-printable")
  )
}

#' @export
`format.ravetools-test_filter` <- function(x, ...) {
  rcon <- x$reciprocal_condition
  rcond_passed <- rcon > .Machine$double.eps
  attenuation <- x$attenuation
  sample_rate <- x$fs
  if(length(sample_rate) == 1 && is.numeric(sample_rate) && !is.na(sample_rate)) {
    nyquist <- sample_rate / 2
    unit <- "Hz"
  } else {
    nyquist <- 1
    unit <- "xPi rad/s"
  }

  if(!isTRUE(rcon > .Machine$double.eps)) {
    warn <- c("", "WARNING: ", " * Unstable autoregressive (AR) polynomial coefficients")
  } else {
    warn <- NULL
  }
  c(
    "<RAVE filter quality test>",
    "Attenuation: ",
    sprintf("  Freq=%.2g %s, mag=%.4g dB (expected=%.4g dB)",
            attenuation$frequency * nyquist, unit, attenuation$filtmag_db, attenuation$expected_db),
    sprintf("Reciprocal condition number: %.2g %s .Machine$double.eps",
            rcon, ifelse(rcond_passed, ">", "<")),
    warn
  )
}
