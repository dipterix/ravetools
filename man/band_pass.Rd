% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/band-pass.R
\name{band_pass}
\alias{band_pass}
\alias{band_pass1}
\alias{band_pass2}
\title{Band-pass signals}
\usage{
band_pass1(x, sample_rate, lb, ub, domain = 1, ...)

band_pass2(
  x,
  sample_rate,
  lb,
  ub,
  order,
  method = c("fir", "butter"),
  direction = c("both", "forward", "backward"),
  window = "hamming",
  ...
)
}
\arguments{
\item{x}{input signals, numeric vector or matrix. \code{x} must be
row-major if input is a matrix: each row is a channel, and each column is
a time-point.}

\item{sample_rate}{sampling frequency}

\item{lb}{lower frequency bound of the band-passing filter, must be positive}

\item{ub}{upper frequency bound of the band-passing filter, must be greater
than the lower bound and smaller than the half of sampling frequency}

\item{domain}{1 if \code{x} is in time-domain, or 0 if \code{x} is in
frequency domain}

\item{...}{ignored}

\item{order}{the order of the filter, must be positive integer and be
less than one-third of the sample rate}

\item{method}{filter type, choices are \code{'fir'} and \code{'butter'}}

\item{direction}{filter direction, choices are \code{'forward'},
\code{'backward'}, and \code{'both'} directions}

\item{window}{window type, can be a character, a function, or a vector.
For character, \code{window} is a function name in the
\code{signal} package, for example, \code{'hanning'}; for a function,
\code{window} takes one integer argument and returns a numeric vector
with length of that input; for vectors, \code{window} is a numeric vector
o length \code{order+1}.}
}
\value{
Filtered signals, vector if \code{x} is a vector, or matrix of
the same dimension as \code{x}
}
\description{
Band-pass signals
}
\examples{


t <- seq(0, 1, by = 0.0005)
x <- sin(t * 0.4 * pi) + sin(t * 4 * pi) + 2 * sin(t * 120 * pi)

oldpar <- par(mfrow = c(2, 2), mar = c(3.1, 2.1, 3.1, 0.1))
# ---- Using band_pass1 ------------------------------------------------

y1 <- band_pass1(x, 2000, 0.1, 1)
y2 <- band_pass1(x, 2000, 1, 5)
y3 <- band_pass1(x, 2000, 10, 80)

plot(t, x, type = 'l', xlab = "Time", ylab = "",
     main = "Mixture of 0.2, 2, and 60Hz")
lines(t, y1, col = 'red')
lines(t, y2, col = 'blue')
lines(t, y3, col = 'green')
legend(
  "topleft", c("Input", "Pass: 0.1-1Hz", "Pass 1-5Hz", "Pass 10-80Hz"),
  col = c(par("fg"), "red", "blue", "green"), lty = 1,
  cex = 0.6
)

# plot pwelch
pwelch(x, fs = 2000, window = 4000, noverlap = 2000, plot = 1)
pwelch(y1, fs = 2000, window = 4000, noverlap = 2000,
       plot = 2, col = "red")
pwelch(y2, fs = 2000, window = 4000, noverlap = 2000,
       plot = 2, col = "blue")
pwelch(y3, fs = 2000, window = 4000, noverlap = 2000,
       plot = 2, col = "green")


# ---- Using band_pass2 with FIR filters --------------------------------

order <- floor(2000 / 3)
z1 <- band_pass2(x, 2000, 0.1, 1, method = "fir", order = order)
z2 <- band_pass2(x, 2000, 1, 5, method = "fir", order = order)
z3 <- band_pass2(x, 2000, 10, 80, method = "fir", order = order)

plot(t, x, type = 'l', xlab = "Time", ylab = "",
     main = "Mixture of 0.2, 2, and 60Hz")
lines(t, z1, col = 'red')
lines(t, z2, col = 'blue')
lines(t, z3, col = 'green')
legend(
  "topleft", c("Input", "Pass: 0.1-1Hz", "Pass 1-5Hz", "Pass 10-80Hz"),
  col = c(par("fg"), "red", "blue", "green"), lty = 1,
  cex = 0.6
)

# plot pwelch
pwelch(x, fs = 2000, window = 4000, noverlap = 2000, plot = 1)
pwelch(z1, fs = 2000, window = 4000, noverlap = 2000,
       plot = 2, col = "red")
pwelch(z2, fs = 2000, window = 4000, noverlap = 2000,
       plot = 2, col = "blue")
pwelch(z3, fs = 2000, window = 4000, noverlap = 2000,
       plot = 2, col = "green")

# ---- Clean this demo --------------------------------------------------
par(oldpar)

}
