
#' Decimate with 'FIR' or 'IIR' filter
#' @param x signal to be decimated
#' @param q integer factor to down-sample by
#' @param n filter order used in the down-sampling; default is \code{30}
#' if \code{ftype='fir'}, or \code{8} if \code{ftype='iir'}
#' @param ftype filter type, choices are \code{'fir'} (default) and
#' \code{'iir'}
#' @returns Decimated signal
#' @details This function is migrated from \code{gsignal} package,
#' but with padding and indexing fixed. The results agree with 'Matlab'.
#'
#' @examples
#'
#' x <- 1:100
#' y <- decimate(x, 2, ftype = "fir")
#' y
#'
#' # compare with signal package
#' z <- gsignal::decimate(x, 2, ftype = "fir")
#'
#' # Compare decimated results
#' plot(x, type = 'l')
#' points(seq(1,100, 2), y, col = "green")
#' points(seq(1,100, 2), z, col = "red")
#'
#'
#' @export
decimate <- function (
  x, q, n = if (ftype == "iir") 8 else 30, ftype = "fir") {
  if (q != round(q))
    stop("decimate only works with integer q.")
  l_x <- length(x)


  if(ftype == "fir"){
    npad <- ceiling(n / 2)
    lpad <- 2*x[1] - x[(npad+1):2]
    rpad <- 2*x[l_x] - x[l_x - (1:npad)]
    inp <- c(lpad, x, rpad)

    b <- fir1(n, 1/q)$b
    y <- fftfilt(b, inp)
    y <- y[ceiling(npad + n/2) + (1:l_x)]
    y <- y[seq(1, length(x), by = q)]
  } else {
    # y <- gsignal::decimate(x, q, n, "iir")

    rip <- 0.05
    w <- 0.8 / q
    ba <- gsignal::cheby1(n, rip, w)

    while(
      n > 1 &&
      (
        all(ba$b == 0) ||
        abs(filtmag_db(ba$b, ba$a, w) + rip) > 1e-6
      )
    ) {
      n <- n - 1
      ba <- gsignal::cheby1(n, rip, w)
    }

    y <- filtfilt(ba$b, ba$a, x)

    ny <- ceiling(l_x / q)

    nbeg <- q - (q * ny - l_x)

    y <- y[seq(nbeg, l_x, by = q)]
  }

  y

}

