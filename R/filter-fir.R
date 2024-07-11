


scale_filter <- function(b, fband, freq, len) {

  if( fband ) {
    b <- b / sum(b)
  } else {
    if( freq[[4]] == 1 ) {
      # length(w) == 1
      f0 <- 1
    } else {
      f0 <- mean(freq[c(3,4)])
    }
    # b = b / abs( exp(-1i*2*pi*(0:L-1)*(f0/2))*(b.') );
    b <- b / Mod(sum(
      exp(-1i*2*pi*seq(0, len - 1)*(f0/2)) * b
    ))
  }

  b

}

fir_validate <- function(n, max_freq, magnitude, odd_allowed = FALSE) {

  # order
  if(missing(n) || length(n) != 1 || !is.numeric(n) || n - round(n) != 0 || n <= 0 ) {
    stop("Filter order `n` must be a positive integer")
  }

  n_is_odd <- n %% 2 == 1

  if( magnitude[[length(magnitude)]] != 0 && max_freq == 1 && n_is_odd && !odd_allowed ) {
    warning("IR filters must have a gain of zero at the Nyquist frequency. Increasing the order by 1")
    n <- n + 1
  }

  n
}


#' Window-based \code{FIR} filter design
#' @description
#' Generate a \code{fir1} filter that is checked against \code{Matlab}
#' \code{fir1} function.
#'
#' @param n filter order
#' @param w band edges, non-decreasing vector in the range 0 to 1, where 1 is
#' the \code{Nyquist} frequency. A scalar for high-pass or low-pass filters,
#' a vector pair for band-pass or band-stop, or a vector for an
#' alternating pass/stop filter.
#' @param type type of the filter, one of \code{"low"} for a low-pass filter,
#' \code{"high"} for a high-pass filter, \code{"stop"} for a stop-band
#' (band-reject) filter, \code{"pass"} for a pass-band filter, \code{"DC-0"}
#' for a band-pass as the first band of a multi-band filter, or \code{"DC-1"}
#' for a band-stop as the first band of a multi-band filter; default \code{"low"}
#' @param window smoothing window function or a numerical vector. The filter is
#' the same shape as the smoothing window. When \code{window} is a function,
#' \code{window(n+1)} will be called, otherwise the length of the window
#' vector needs to have length of \code{n+1}; default: \code{hamming}
#' @param scale whether to scale the filter; default is true
#' @param hilbert whether to use 'Hilbert' transformer; default is false
#' @returns The \code{FIR} filter coefficients with class \code{'Arma'}.
#' The moving average coefficient is a vector of length \code{n+1}.
#' @export
fir1 <- function(
    n, w, type = c("low", "high", "stop", "pass", "DC-0", "DC-1"),
    window = hamming, scale = TRUE, hilbert = FALSE
) {
  type <- match.arg(type)

  nw <- length(w)
  if(!nw || any(w < 0 | w > 1)) {
    stop("`fir1`: w must be a real vector from range 0.0 to 1.0")
  }
  if(nw > 1 && any(diff(w) < 0)) {
    stop("`fir1`: w must be non-decreasing")
  }

  if(w[[length(w)]] >= 1) {
    w <- w[-length(w)]
    if(length(w) == 1) {
      type <- "high"
    } else {
      type <- "DC-0"
    }
  } else if(missing(type)) {
    if(nw == 1) {
      type <- "low"
    } else if (nw == 2) {
      type <- 'pass'
    } else {
      type <- "DC-0"
    }
  } else {
    type <- match.arg(type)
  }

  nw <- length(w)

  nbands <- nw + 1

  # make sure default 3 band filter is bandpass
  if( type == "pass" ) {
    if( nbands > 2 ) {
      type <- "DC-0"
    }
  }

  # frequency vector
  freq <- c(0, rep(w, each = 2), 1)

  # magnitude vector
  use_first_band <- !type %in% c("DC-0", "high")
  magnitude <- rep( (use_first_band + seq(0, nbands - 1)) %% 2, each = 2)

  n <- fir_validate(n, freq[[length(freq)]], magnitude, odd_allowed = hilbert)

  # filter length is order+1
  filter_len <- n + 1

  if(is.function(window)) {
    wind <- window(filter_len)
    if( length(wind) != filter_len ) {
      stop("`fir1`: window length is invalid given requested length: (n+1)")
    }
  } else {
    if(length(window) != filter_len) {
      stop("`fir1`: window must be either a function or a vector of length (n+1)")
    }
    wind <- window
  }

  if(hilbert) {
    hh <- firls(n, freq, magnitude, ftype = "hilbert")$b
  } else {
    hh <- firls(n, freq, magnitude)$b
  }

  b <- wind * hh

  if(scale) {
    b <- scale_filter(b, use_first_band, freq, filter_len)
  }

  return(structure(list(b = b, a = 1), class = "Arma"))
}

sinc <- function(x) {
  v <- pi * x
  re <- sin(v) / v
  if(length(v)) {
    re[ v == 0 ] <- 1
  }
  re
}

#' Least-squares linear-phase \code{FIR} filter design
#' @description
#' Produce a linear phase filter from the weighted mean squared such that error
#' in the specified bands is minimized.
#' @param N filter order, must be even (if odd, then will be increased by one)
#' @param freq vector of frequency points in the range from 0 to 1, where 1
#' corresponds to the \code{Nyquist} frequency.
#' @param A vector of the same length as \code{freq} containing the desired
#' amplitude at each of the points specified in \code{freq}.
#' @param W weighting function that contains one value for each band that
#' weights the mean squared error in that band. \code{W} must be half the
#' length of \code{freq}.
#' @param ftype transformer type; default is \code{""}; alternatively,
#' \code{'h'} or \code{'hilbert'} for 'Hilbert' transformer.
#' @returns The \code{FIR} filter coefficients with class \code{'Arma'}.
#' The moving average coefficient is a vector of length \code{n+1}.
#' @export
firls <- function(N, freq, A, W = NULL, ftype = "") {

  if (!is.numeric(N) || length(N) != 1 || !is.finite(N) || N <= 0) {
    stop("Invalid input for N.")
  }
  N <- as.integer(N)

  if (!is.numeric(freq) || !is.numeric(A) || length(freq) != length(A) || any(freq < 0) || any(freq > 1)) {
    stop("Invalid input for freq or A.")
  }
  freq <- as.double(freq)
  A <- as.double(A)

  if (!is.null(W)) {
    if (length(W) == 0) {
      W <- rep(1, floor(length(freq) / 2))
    } else {
      W <- as.double(W)
    }
  } else {
    W <- rep(1, floor(length(freq) / 2))
  }

  if (identical(ftype, "h") || identical(ftype, "hilbert")) {
    filtype <- 1
    differ <- 0
  } else if (identical(ftype, "d") || identical(ftype, "differentiator")) {
    stop("Filter type 'differentiator' has not been implemented yet")
    filtype <- 1
    differ <- 1
  } else {
    filtype <- 0
    differ <- 0
  }

  if (length(W) == 0) {
    weight <- rep(1, floor(length(freq) / 2))
  } else {
    weight <- W
  }

  freq_length <- length(freq)
  wt <- abs(sqrt(weight))

  N <- N + 1

  Fn <- freq / 2

  dF <- diff(Fn)

  if (any(dF < 0)) {
    stop("Invalid frequency vector.")
  }

  lendF <- freq_length - 1
  fullband <- FALSE

  if (lendF > 1) {
    fullband <- all(dF == 0)
  }

  tempW <- wt - wt[1]
  constant_weights <- sum(tempW) == 0

  L <- (N - 1) / 2
  Nodd <- (N %% 2 == 1)

  b0 <- 0

  if (filtype == 0) {
    if (!Nodd) {
      m <- (0:L) + 0.5
    } else {
      m <- 0:L
    }

    k <- m

    need_matrix <- !fullband || !constant_weights

    if (need_matrix) {
      result <- initMatrices(m)
      I1 <- result$I1
      I2 <- result$I2
      G <- result$G
    }

    if (Nodd) {
      k <- k[-1]
      b0 <- 0
    }

    b <- numeric(length(k))

    for (s in seq(1, length(Fn), by = 2)) {
      if(Fn[s + 1] == Fn[s]) {
        next
      }
      m_s <- (A[s + 1] - A[s]) / (Fn[s + 1] - Fn[s])
      b1 <- A[s] - m_s * Fn[s]

      if (Nodd) {
        b0 <- b0 + (b1 * (Fn[s + 1] - Fn[s]) + m_s / 2 * (Fn[s + 1]^2 - Fn[s]^2)) * wt[(s + 1) / 2]^2
      }

      b <- b + (m_s / (4 * pi^2) * (cos(2 * pi * k * Fn[s + 1]) - cos(2 * pi * k * Fn[s])) / (k^2)) * wt[(s + 1) / 2]^2
      b <- b + (Fn[s + 1] * (m_s * Fn[s + 1] + b1) * sinc(2 * k * Fn[s + 1]) -
                  Fn[s] * (m_s * Fn[s] + b1) * sinc(2 * k * Fn[s])) * wt[(s + 1) / 2]^2

      if (need_matrix) {
        G <- G + (0.5 * Fn[s + 1] * (sinc(2 * I1 * Fn[s + 1]) + sinc(2 * I2 * Fn[s + 1])) -
                    0.5 * Fn[s] * (sinc(2 * I1 * Fn[s]) + sinc(2 * I2 * Fn[s]))) * wt[(s + 1) / 2]^2
      }
    }

    if (Nodd) {
      b <- c(b0, b)
    }

    if (need_matrix) {
      a <- qr.solve(G, b, tol = 1e-30)
    } else {
      a <- (wt[1]^2) * 4 * b
      if (Nodd) {
        a[1] <- a[1] / 2
      }
    }

    if (Nodd) {
      h <- c(a[seq(L+1,2, by = -1)] / 2, a[1], a[2:(L + 1)] / 2)
    } else {
      h <- 0.5 * c(rev(a), a)
    }

  } else if (filtype == 1) {
    if (differ) {
      do_weight <- as.double(abs(A[seq(2, length(A), by = 2)]) + abs(A[seq(1, length(A), by = 2)]) > 0)
    } else {
      do_weight <- numeric(length(Fn))
    }

    need_matrix <- !fullband || any(do_weight) || !constant_weights

    if (Nodd) {
      m <- 1:L
    } else {
      m <- 0:L + 0.5
    }

    k <- m
    b <- numeric(length(k))

    if (need_matrix) {
      result <- initMatrices(m)
      I1 <- result$I1
      I2 <- result$I2
      G <- result$G
    } else {
      G <- matrix(0, nrow = 0, ncol = 0)
    }

    for (s in seq(1, length(Fn), by = 2)) {
      if( Fn[s + 1] == Fn[s] ) { next }
      if (do_weight[(s + 1) / 2]) {
        if (Fn[s] == 0) {
          Fn[s] <- 1e-5
        }

        m_s <- (A[s + 1] - A[s]) / (Fn[s + 1] - Fn[s])
        b1 <- A[s] - m_s * Fn[s]

        snint1 <- sineint(2 * pi * k * Fn[s + 1]) - sineint(2 * pi * k * Fn[s])
        csint1 <- Re((-1/2) * (pracma::expint(1i * 2 * pi * k * Fn[s + 1]) +
                                 pracma::expint(-1i * 2 * pi * k * Fn[s + 1]) -
                                 pracma::expint(1i * 2 * pi * k * Fn[s]) -
                                 pracma::expint(-1i * 2 * pi * k * Fn[s])))

        b <- b + (m_s * snint1 +
                    b1 * 2 * pi * k * (-sinc(2 * k * Fn[s + 1]) + sinc(2 * k * Fn[s]) + csint1)) *
          wt[(s + 1) / 2]^2

        snint1 <- sineint(2 * pi * Fn[s + 1] * (-I2))
        snint2 <- sineint(2 * pi * Fn[s + 1] * I1)
        snint3 <- sineint(2 * pi * Fn[s] * (-I2))
        snint4 <- sineint(2 * pi * Fn[s] * I1)

        G <- G - ((-1/2) * (cos(2 * pi * Fn[s + 1] * (-I2)) / Fn[s + 1] -
                              2 * snint1 * pi * I2 -
                              cos(2 * pi * Fn[s + 1] * I1) / Fn[s + 1] +
                              2 * snint2 * pi * I1) -
                    (-1/2) * (cos(2 * pi * Fn[s] * (-I2)) / Fn[s] -
                                2 * snint3 * pi * I2 -
                                cos(2 * pi * Fn[s] * I1) / Fn[s] +
                                2 * snint4 * pi * I1)) * wt[(s + 1) / 2]^2
      } else {
        m_s <- (A[s + 1] - A[s]) / (Fn[s + 1] - Fn[s])
        b1 <- A[s] - m_s * Fn[s]

        b <- b + (m_s / (4 * pi^2) * (sin(2 * pi * k * Fn[s + 1]) - sin(2 * pi * k * Fn[s])) / (k^2)) * wt[(s + 1) / 2]^2
        b <- b + ((m_s * Fn[s] + b1) * cos(2 * pi * k * Fn[s]) -
                    (m_s * Fn[s + 1] + b1) * cos(2 * pi * k * Fn[s + 1])) / (2 * pi * k) * wt[(s + 1) / 2]^2

        if (need_matrix) {
          G <- G + (0.5 * Fn[s + 1] * (sinc(2 * I1 * Fn[s + 1]) - sinc(2 * I2 * Fn[s + 1])) -
                      0.5 * Fn[s] * (sinc(2 * I1 * Fn[s]) - sinc(2 * I2 * Fn[s]))) * wt[(s + 1) / 2]^2
        }
      }
    }

    if (need_matrix) {
      a <- solve(G, b)
    } else {
      a <- -4 * b * wt[1]^2
    }

    if (Nodd) {
      h <- 0.5 * c(rev(a), 0, -a)
    } else {
      h <- 0.5 * c(rev(a), -a)
    }

    if (differ) {
      h <- -h
    }

  } else {
    h <- 0
  }

  return(structure(list(b = h, a = 1), class = "Arma"))
}


initMatrices <- function(m) {
  k <- m
  x <- matrix(k, ncol = length(m), nrow = length(m), byrow = TRUE)
  y <- matrix(m, ncol = length(k), nrow = length(k), byrow = FALSE)

  I1 <- x + y
  I2 <- x - y
  G <- matrix(0, ncol = length(m), nrow = length(m))

  return(list(I1 = I1, I2 = I2, G = G))
}


# Might not be implemented correctly, not used
fir2 <- function (n, f, m, grid_n = 512,
                  ramp_n = grid_n/20,
                  window = hamming(n + 1))
{
  t <- length(f)
  if (t < 2 || f[1] != 0 || f[t] != 1 || any(diff(f) < 0))
    stop("frequency must be nondecreasing starting from 0 and ending at 1")
  if (t != length(m))
    stop("frequency and magnitude vectors must be the same length")
  if (length(grid_n) > 1 || length(ramp_n) > 1)
    stop("grid_n and ramp_n must be integers")
  if (length(window) != n + 1)
    stop("window must be of length n+1")
  if (2 * grid_n < n + 1)
    grid_n <- 2^ceiling(log2(abs(n + 1)))
  if (ramp_n > 0) {
    basef <- f
    basem <- m
    idx <- which(diff(f) == 0)
    f[idx] <- f[idx] - ramp_n/grid_n/2
    f[idx + 1] <- f[idx + 1] + ramp_n/grid_n/2
    f <- c(f, basef[idx])
    f[f < 0] <- 0
    f[f > 1] <- 1
    f <- sort(unique(c(f, basef[idx])))
    m <- approx(basef, basem, f, ties = "ordered")$y
  }
  grid <- approx(f, m, seq(0, 1, length = grid_n + 1), ties = "ordered")$y
  if ((n%%2) == 0) {
    b <- ifft(c(grid, grid[seq(grid_n, 2, by = -1)]))
    mid <- (n + 1)/2
    b <- Re(c(b[(2 * grid_n - floor(mid) + 1):(2 * grid_n)],
              b[1:ceiling(mid)]))
  } else {
    b <- ifft(c(grid, rep(0, grid_n * 2), grid[seq(grid_n, 2, by = -1)]))
    b <- 2 * Re(c(b[seq(length(b) - n + 1, length(b), by = 2)],
                  b[seq(2, n + 2, by = 2)]))
  }
  b <- b * window
  b
}

