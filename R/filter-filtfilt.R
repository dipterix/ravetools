

filter_initialize <- function(b, a, x) {
  # DIPSAUS DEBUG START
  # filter <- ravetools::butter(1, c(50 / 250), "low")
  # b <- filter$b
  # a <- filter$a
  # x <- sample_signal(500)

  if(a[[1]] != 0) {
    a <- a / a[[1]]
    b <- b / a[[1]]
  }

  # make sure a, b share the same order
  na <- length(a)
  nb <- length(b)

  if( na > nb ) {
    b <- c(b, rep(0, na - nb))
    n <- na
  } else {
    a <- c(a, rep(0, nb - na))
    n <- nb
  }

  # length of edge transients
  nf <- max(1, 3 * (n - 1))

  if(length(x) <= nf) {
    stop("Cannot apply the filter: the input signal is too short")
  }

  # n is nfilt in matlab filtfilt -> getCoeffsAndInitialConditions
  # nf is nfact in matlab filtfilt -> getCoeffsAndInitialConditions

  # compute the initial condition if n > 1
  if( n > 2 ) {
    z1 <- diag(1, n - 1) - cbind( -a[-1], rbind(diag(1, n - 2), 0) )
    z2 <- b[-1] - b[1] * a[-1]
    z <- qr.solve(z1, z2, tol = 1e-30)

  } else if ( n == 2 ) {
    z1 <- 1 + a[-1]
    z2 <- b[-1] - b[1] * a[-1]
    z <- z2 / z1
  } else {
    z <- numeric(0)
  }
  list(
    a = a,
    b = b,
    z = z,
    nfilt = n,
    nfact = nf
  )
}

filtfilt_naive <- function(b, a, y, z, nfact) {

  # pad edges with reverse signals
  ny <- length(y)
  ytmp <- c(
    2.0 * y[[1]] - y[seq(nfact + 1, 2, by = -1)],
    y,
    2.0 * y[[ny]] - y[seq(ny - 1, ny - nfact, by = -1)]
  )

  ytmp <- filter_signal(b = b, a = a, x = ytmp, z = z * ytmp[[1]])[[1]]
  ytmp <- rev(ytmp)
  ytmp <- filter_signal(b = b, a = a, x = ytmp, z = z * ytmp[[1]])[[1]]

  re <- ytmp[seq(length(ytmp) - nfact, nfact + 1, by = -1)]
  return(re)
}

filtfilt_naive2 <- function(b, a, y, z, nfact) {

  xt <- 2.0 * y[[1]] - y[seq(nfact + 1, 2, by = -1)]
  z0 <- filter_signal(b = b, a = a, x = xt, z = z * xt[[1]])[[2]]
  yc2 <- filter_signal(b = b, a = a, x = y, z = z0)
  z0 <- yc2[[2]]
  yc2 <- yc2[[1]]

  ny <- length(y)
  xt <- 2.0 * y[[ny]] - y[seq(ny - 1, ny - nfact, by = -1)]
  yc3 <- filter_signal(b = b, a = a, x = xt, z = z0)[[1]]

  z0 <- filter_signal(b = b, a = a, x = rev(yc3), z = z * yc3[[length(yc3)]])[[2]]
  yc5 <- filter_signal(b = b, a = a, x = rev(yc2), z = z0)[[1]]

  return(rev( yc5 ))
}

#' @title Forward and reverse filter a one-dimensional signal
#' @description The result has been tested against 'Matlab' \code{filtfilt}
#' function. Currently this function only supports one filter at a time.
#' @param b one-dimensional real numerical vector, the moving-average
#' coefficients of an \code{ARMA} filter
#' @param a the auto-regressive (recursive) coefficients of an \code{ARMA} filter
#' @param x numerical vector input (real value)
#' @returns The filtered signal, normally the same length as the input signal
#' \code{x}.
#' @examples
#'
#' t <- seq(0, 1, by = 0.01)
#' x <- sin(2 * pi * t * 2.3)
#' bf <- gsignal::butter(2, c(0.15, 0.3))
#'
#' res <- filtfilt(bf$b, bf$a, x)
#'
#' ## Matlab (2022a) equivalent:
#' # t = [0:0.01:1];
#' # x = sin(2 * pi * t * 2.3);
#' # [b,a] = butter(2,[.15,.3]);
#' # res = filtfilt(b, a, x)
#'
#' @export
filtfilt <- function(b, a, x) {

  nx <- length(x)

  # DIPSAUS DEBUG START
  # filter <- ravetools::butter(5, c(1 / 250, 50 / 250), "pass")
  # b <- filter$b
  # a <- filter$a
  # x <- sample_signal(500)
  init <- filter_initialize(b, a, x)

  if(nx < 10000) {
    re <- filtfilt_naive(b = init$b, a = init$a, y = x, z = init$z, nfact = init$nfact)
  } else {
    re <- filtfilt_naive2(b = init$b, a = init$a, y = x, z = init$z, nfact = init$nfact)
  }
  re
}


#' @title Filter one-dimensional signal
#' @description The function is written from the scratch. The result has been
#' compared against the 'Matlab' \code{filter} function with one-dimensional
#' real inputs. Other situations such as matrix \code{b} or multi-dimensional
#' \code{x} are not implemented. For double filters (forward-backward),
#' see \code{\link{filtfilt}}.
#' @param b one-dimensional real numerical vector, the moving-average
#' coefficients of an \code{ARMA} filter
#' @param a the auto-regressive (recursive) coefficients of an \code{ARMA} filter
#' @param x numerical vector input (real value)
#' @param z initial condition, must have length of \code{n-1}, where \code{n}
#' is the maximum of lengths of \code{a} and \code{b}; default is all zeros
#' @returns A list of two vectors: the first vector is the filtered signal;
#' the second vector is the final state of \code{z}
#'
#' @examples
#'
#'
#' t <- seq(0, 1, by = 0.01)
#' x <- sin(2 * pi * t * 2.3)
#' bf <- gsignal::butter(2, c(0.15, 0.3))
#'
#' res <- filter_signal(bf$b, bf$a, x)
#' y <- res[[1]]
#' z <- res[[2]]
#'
#' ## Matlab (2022a) equivalent:
#' # t = [0:0.01:1];
#' # x = sin(2 * pi * t * 2.3);
#' # [b,a] = butter(2,[.15,.3]);
#' # [y,z] = filter(b, a, x)
#'
#'
#' @export
filter_signal <- function(b, a, x, z) {

  na <- length(a)
  nb <- length(b)

  if( na > nb ) {
    b <- c(b, rep(0, na - nb))
    n <- na
  } else {
    a <- c(a, rep(0, nb - na))
    n <- nb
  }

  if(missing(z)) {
    z <- rep(0.0, n - 1)
  } else {
    if(length(z) < n-1) {
      stop(sprintf("`filter`: initial condition `z` must have length >= %d", n-1))
    }
  }

  if(!is.double(a)) {
    a <- as.double(a)
  }
  if(!is.double(b)) {
    b <- as.double(b)
  }
  if(!is.double(z)) {
    z <- as.double(z)
  }
  if(!is.double(x)) {
    x <- as.double(x)
  }
  return(cpp_filter(b, a, x, z))
}

