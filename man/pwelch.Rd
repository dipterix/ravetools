% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pwelch.R
\name{pwelch}
\alias{pwelch}
\alias{print.ravetools-pwelch}
\alias{plot.ravetools-pwelch}
\alias{mv_pwelch}
\title{Calculate 'Welch Periodogram'}
\usage{
pwelch(
  x,
  fs,
  window = 64,
  noverlap = window/2,
  nfft = "auto",
  window_family = hamming,
  col = "black",
  xlim = NULL,
  ylim = NULL,
  main = "Welch periodogram",
  plot = 0,
  log = c("xy", "", "x", "y"),
  ...
)

\method{print}{`ravetools-pwelch`}(x, ...)

\method{plot}{`ravetools-pwelch`}(
  x,
  log = c("xy", "x", "y", ""),
  se = FALSE,
  xticks,
  type = "l",
  add = FALSE,
  col = graphics::par("fg"),
  col.se = "orange",
  alpha.se = 0.5,
  lty = 1,
  lwd = 1,
  cex = 1,
  las = 1,
  main = "Welch periodogram",
  xlab,
  ylab,
  xlim = NULL,
  ylim = NULL,
  xaxs = "i",
  yaxs = "i",
  xline = 1.2 * cex,
  yline = 2 * cex,
  mar = c(2.6, 3.8, 2.1, 0.6) * (0.5 + cex/2),
  mgp = cex * c(2, 0.5, 0),
  tck = -0.02 * cex,
  grid = TRUE,
  ...
)

mv_pwelch(
  x,
  margin,
  fs,
  window = 64,
  noverlap = window/2,
  nfft = "auto",
  window_family = hamming
)
}
\arguments{
\item{x}{numerical vector or a row-major vector, signals.
If \code{x} is a matrix, then each row is a channel. For \code{plot}
function, \code{x} is the instance returned by \code{pwelch} function.}

\item{fs}{sample rate, average number of time points per second}

\item{window}{window length in time points, default size is \code{64}}

\item{noverlap}{overlap between two adjacent windows, measured in time
points; default is half of the \code{window}}

\item{nfft}{number of points in window function; default is automatically
determined from input data and window, scaled up to the nearest power of 2}

\item{window_family}{function generator for generating filter windows,
default is \code{\link{hamming}}. This can be any window function listed in
the filter window family, or any window generator function from package
\code{gsignal}. Default is \code{\link{hamming}}. For 'iEEG' users, both
\code{hamming} and \code{\link{blackmanharris}} are offered by 'EEG-lab';
while \code{blackmanharris} offers better attenuation than Hamming windows,
it also has lower spectral resolution. \code{hamming} has a 42.5 dB side-lobe
attenuation. This may mask spectral content below this value (relative
to the peak spectral content). Choosing different windows enables
you to make trade-off between resolution (e.g., using a rectangular
window) and side-lobe attenuation (e.g., using a \code{\link{hanning}}
window)}

\item{col, xlim, ylim, main, type, cex, las, xlab, ylab, lty, lwd, xaxs, yaxs, mar, mgp, tck}{parameters passed to \code{\link[graphics]{plot.default}}}

\item{plot}{integer, whether to plot the result or not; choices are \code{0}, no plot; \code{1} plot on a new canvas; \code{2} add to existing canvas}

\item{log}{indicates which axis should be \code{log10}-transformed, used by the plot function. For \code{'x'} axis, it's \code{log10}-transform; for \code{'y'} axis, it's \code{10log10}-transform (decibel unit). Choices are \code{"xy"}, \code{"x"}, \code{"y"}, and \code{""}.}

\item{...}{will be passed to \code{plot.pwelch} or ignored}

\item{se}{logical or a positive number indicating whether to plot standard
error of mean; default is false. If provided with a number, then a multiple
of standard error will be drawn. This option is only available when power
is in log-scale (decibel unit)}

\item{xticks}{ticks to show on frequency axis}

\item{add}{logical, whether the plot should be added to existing canvas}

\item{col.se, alpha.se}{controls the color and opacity of the standard error}

\item{xline, yline}{controls how close the axis labels to the corresponding axes}

\item{grid}{whether to draw rectangular grid lines to the plot; only
respected when \code{add=FALSE}; default is true}

\item{margin}{the margin in which \code{pwelch} should be applied to}
}
\value{
A list with class \code{'ravetools-pwelch'} that contains the
following items:
\describe{
\item{\code{freq}}{frequencies used to calculate the 'periodogram'}
\item{\code{spec}}{resulting spectral power for each frequency}
\item{\code{window}}{window function (in numerical vector) used}
\item{\code{noverlap}}{number of overlapping time-points between two adjacent windows}
\item{\code{nfft}}{number of basis functions}
\item{\code{fs}}{sample rate}
\item{\code{x_len}}{input signal length}
\item{\code{method}}{a character string \code{'Welch'}}
}
}
\description{
\code{pwelch} is for single signal trace only; \code{mv_pwelch}
is for multiple traces. Currently \code{mv_pwelch} is experimental and
should not be called directly.
}
\examples{

x <- rnorm(1000)
pwel <- pwelch(x, 100)
pwel

plot(pwel, log = "xy")

}
