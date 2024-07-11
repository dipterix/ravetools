
butter_cutoff <- function(type, n, w, r) {
  # DIPSAUS DEBUG START
  # w <- 60 / 250
  # type <- "high"
  # w <- c(1, 60) / 250
  # type <- "pass"
  # r <- 6
  # n <- 6

  Ws <- sort(w)
  Rs <- r
  Rc <- 10 * log10(2)

  # Wpw <- (2/T) * tan(pi * Wp/T)
  Wsw <- tan(pi * Ws / 2)

  qs <- log(10 ^ (Rs / 10) - 1)
  qc <- log(10 ^ (Rc / 10) - 1)

  # n <- 0.5 * (qs - qp)/log(ws/wp)
  # log(ws/wc) =
  log_ws_wc <- (qs - qc) / ( 2 * n )
  ws_wc <- exp( log_ws_wc )

  switch(
    type,
    "low" = {
      Wcw <- Wsw / ws_wc
    },
    "high" = {
      # wc <- Wsw
      # ws <- Wcw
      Wcw <- Wsw * ws_wc
    },
    "pass" = {
      w02 <- Wsw[[1]] * Wsw[[2]]
      # wp <- Wpw[2] - Wpw[1]
      ws <- Wsw[2] - Wsw[1]

      # ws <- ws/wp
      # wp <- 1
      wc <- ws / ws_wc

      # solve:
      # Wcw[2] - Wcw[1] = wc
      # Wcw[2] * Wcw[1] = w02
      #
      # i.e.
      # Wcw[1]^2 + wc * Wcw[1] - w02 = 0
      Wcw1 <- (-wc + sqrt(wc^2 + 4 * w02)) / 2
      Wcw2 <- Wcw1 + wc
      # 0.005631939 0.441716719
      Wcw <- c(Wcw1, Wcw2)
    },
    "stop" = {
      w02 <- Wsw[1] * Wsw[2]

      # wp <- w02/(Wpw[2] - Wpw[1])
      # ws <- w02/(Wsw[2] - Wsw[1])
      # ws <- ws/wp
      # wp <- 1

      ws <- w02/(Wsw[2] - Wsw[1])
      wc <- ws / ws_wc
      wc <- w02 / wc

      # solve:
      # Wcw[2] - Wcw[1] = wc
      # Wcw[2] * Wcw[1] = w02
      #
      # i.e.
      # Wcw[1]^2 + wc * Wcw[1] - w02 = 0
      Wcw1 <- (-wc + sqrt(wc^2 + 4 * w02)) / 2
      Wcw2 <- Wcw1 + wc
      Wcw <- c(Wcw1, Wcw2)
    }
  )

  Wc <- atan( Wcw ) * 2 / pi

  # filter <- butter(n, Wc, type)
  # hw <- freqz2(filter$b, filter$a, fs = 500)
  # print(hw, cutoffs = -Rs)
  # filtmag_db(filter$b, filter$a, Wc) + Rc
  # filtmag_db(filter$b, filter$a, Ws) + Rs
  Wc

}

butter_check <- function(type, n, w, r) {
  r <- abs(r)
  Wc <- butter_cutoff(type, n, w, r)
  filter <- gsignal::butter(n = n, w = Wc, type = type, output = "Arma")

  half_pow <- 10*log10(0.5)
  check_filter(b = filter$b, a = filter$a, w = c(w, Wc), r_expected = rep(c(-r, half_pow), each = length(w)))
}

#' @title 'Butterworth' filter with maximum order
#' @description
#' Large filter order might not be optimal, but at lease this function
#' provides a feasible upper bound for the order such that the
#' filter has a stable \code{AR} component.
#'
#' @param type filter type
#' @param w scaled frequency ranging from 0 to 1, where 1 is 'Nyquist' frequency
#' @param r decibel attenuation at frequency \code{w}, default is around
#' \code{3 dB} (half power)
#' @param tol tolerance of reciprocal condition number, default is
#' \code{.Machine$double.eps}.
#'
#' @returns 'Butterworth' filter in 'Arma' form.
#'
#' @examples
#'
#' # Find highest order (sharpest transition) of a band-pass filter
#' sample_rate <- 500
#' nyquist <- sample_rate / 2
#'
#' type <- "pass"
#' w <- c(1, 50) / nyquist
#' Rs <- 6     # power attenuation at w
#'
#' # max order filter
#' filter <- butter_max_order(w, "pass", Rs)
#'
#' # -6 dB cutoff should be around 1 ~ 50 Hz
#' diagnose_filter(filter$b, filter$a, fs = sample_rate)
#'
#' @export
butter_max_order <- function(w, type = c("low", "high", "pass", "stop"),
                             r = 10*log10(2), tol = .Machine$double.eps) {

  w <- sort(w)
  r <- abs(r)

  stopifnot2(
    all(is.finite(w) & w >= 0 & w <= 1),
    msg = "`w` must be scaled to [0, 1], where 1 is Nyquist frequency"
  )

  if(type %in% c("pass", "stop")) {
    stopifnot(length(w) == 2)
  } else {
    stopifnot(length(w) == 1)
  }

  validator <- function(n) {
    check <- butter_check(type, n, w, r)
    check$reciprocal_condition > tol
  }

  n <- guess_max_integer(validator, min_v = 1)

  if(is.na(n)) { n <- 1 }

  Wc <- butter_cutoff(type, n, w, r)

  filter <- butter(n, Wc, type)

  filter$n <- n
  filter$Wc <- Wc
  filter$type <- type
  filter$check <- butter_check(type, n, w, r)
  filter

}

