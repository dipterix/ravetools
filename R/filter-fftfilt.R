# ---- fftfilt -----------------------------------------------------------------

fftfilt <- function (b, x, n = NULL, legacy = FALSE) {
  if (length(n) > 1) {
    stop("fftfilt: n has to be a scalar")
  }
  # DIPSAUS DEBUG START
  # b <- rnorm(5)
  # x <- rnorm(500)
  # n <- NULL

  N <- n
  l_b <- length(b)

  if(!is.matrix(x)) {
    l_x <- length(x)
    x <- matrix(x, ncol = 1L)
    is_vector <- TRUE
  } else {
    l_x <- nrow(x)
    is_vector <- FALSE
  }
  if(!length(n)) {
    N <- 2 ^ (ceiling(log(l_x + l_b - 1) / log(2)))
    B <- fft(postpad(b, N))
    if( legacy ) {
      y <- apply(x, 2L, function(xi) {
        # y <- ifft(fft(postpad(x, N)) * B)
        ifft(fft(postpad(xi, N)) * B)
      })
    } else {
      y <- stats::mvfft(stats::mvfft(postpad(x, N)) * B, inverse = TRUE) / N
    }
  } else {
    N <- 2 ^ (ceiling(log(max(n, l_b)) / log(2)))
    L <- N - l_b + 1
    B <- fft(postpad(b, N))
    R <- ceiling(l_x/L)

    y <- array(0.0i, c(l_x, ncol(x)))
    for (r in seq_len(R)) {
      lo <- (r - 1) * L + 1
      hi <- min(r * L, l_x)
      tmp <- x[lo:hi, , drop = FALSE]
      if( legacy ) {
        tmp <- apply(tmp, 2L, function(tmpi) {
          ifft(fft(postpad(tmpi, N)) * B)
        })
      } else {
        tmp <- stats::mvfft(stats::mvfft(postpad(tmp, N)) * B, inverse = TRUE) / N
      }
      hi <- min(lo + N - 1, l_x)
      y[lo:hi,] <- y[lo:hi,] + tmp[1:(hi - lo + 1),]
    }
  }
  y <- y[1:l_x, , drop = is_vector]
  if (is.numeric(b) && is.numeric(x)) {
    y <- Re(y)
  }
  if (!any(as.logical(b - round(b)))) {
    idx <- !any(as.logical(x - round(x)))
    y[idx] <- round(y[idx])
  }
  return(y)

  # b <- rnorm(5)
  # x <- rnorm(500)
  # n <- NULL
  # range(fftfilt(b, x) - fftfilt(b, x, legacy = TRUE))
  # range(fftfilt(b, array(x, c(250,2))) - fftfilt(b, array(x, c(250,2)), legacy = TRUE))
  # fftfilt(b, array(x, c(250,2))) - fftfilt(b, array(x, c(250,2)), legacy = TRUE)

}
