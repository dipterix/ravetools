#' Calculate 'Welch Periodogram'
#' @description \code{pwelch} is for single signal trace only; \code{mv_pwelch}
#' is for multiple traces. Currently \code{mv_pwelch} is experimental and
#' should not be called directly.
#' @param x numerical vector, analog voltage signal
#' @param fs sample rate, average number of time points per second
#' @param window window length in time points, default size is \code{64}
#' @param nfft number of basis functions to apply
#' @param noverlap overlap between two adjacent windows, measured in time points; default is \code{8}
#' @param log indicates which axis should be \code{log10}-transformed, used by the plot function. For \code{'x'} axis, it's \code{log10}-transform; for \code{'y'} axis, it's \code{10log10}-transform (decibel unit). Choices are \code{"xy"}, \code{"x"}, \code{"y"}, and \code{""}.
#' @param plot integer, whether to plot the result or not; choices are \code{0}, no plot; \code{1} plot on a new canvas; \code{2} add to existing canvas
#' @param add logical, whether the plot should be added to existing canvas
#' @param ... will be passed to \code{plot.pwelch} or ignored
#' @param x \code{'pwelch'} instance returned by \code{pwelch} function
#' @param col,xlim,ylim,main,type,cex,cex.main,cex.sub,cex.lab,cex.axis,las,xlab,ylab parameters passed to \code{\link[graphics]{plot.default}}
#' @param margin the margin in which \code{pwelch} should be applied to
#' @return A list with class \code{'ravetools-pwelch'} that contains the
#' following items:
#' \describe{
#' \item{\code{freq}}{frequencies used to calculate the 'periodogram'}
#' \item{\code{spec}}{resulting spectral power for each frequency}
#' \item{\code{window}}{window function (in numerical vector) used}
#' \item{\code{noverlap}}{number of overlapping time-points between two adjacent windows}
#' \item{\code{nfft}}{number of basis functions}
#' \item{\code{fs}}{sample rate}
#' \item{\code{x_len}}{input signal length}
#' \item{\code{method}}{a character string \code{'Welch'}}
#' }
#' @examples
#'
#' x <- rnorm(1000)
#' pwel <- pwelch(x, 100)
#' pwel
#'
#' plot(pwel, log = "xy")
#'
#' @export
pwelch <- function (
  x, fs, window = 64, noverlap = 8, nfft = 256,
  col = 'black', xlim = NULL, ylim = NULL, main = 'Welch periodogram',
  plot = 0, log = c("xy", "", "x", "y"), ...) {


  # list2env(list(window = 64, noverlap = 8, nfft = 256,
  #               col = 'black', xlim = NULL, ylim = NULL, main = 'Welch periodogram',
  #               plot = TRUE, log = 'xy', spec_func = stats::spectrum, cex = 1), .GlobalEnv)

  x <- as.vector(x)
  x_len <- length(x)

  nfft <- max(min(nfft, length(x)), window)

  window <- hanning(window)

  # window_norm = norm(window, '2')
  window_len <- length(window)

  # normalization <- mean(window^2)

  step <- max(floor(window_len - noverlap + 0.99), 1)

  ## Average the slices
  offset <- seq(1, x_len-window_len+1, by = step)

  N <- length(offset)

  re <- sapply(seq_len(N), function(i){
    a <- detrend_naive(x[offset[i] - 1 + seq_len(window_len)])
    postpad(a$Y * window, nfft)
    # # HermConj = 0 : without the "Hermitian" redundancy
    # a = fftwtools::fftw_r2c(postpad(a$Y * window, nfft), HermConj = 0)
    # Mod(a)^2
  })

  re <- Mod(mvfftw_r2c(re, inplace = TRUE))^2

  NN <- floor((nfft + 1)/2)
  spec <- rowMeans(re) / (window_len / 2)^2
  spec <- spec[seq_len(NN)]
  freq <- seq(1, fs / 2, length.out = NN)

  res <- structure(list(
    freq = freq,
    spec = spec,
    window = window,
    noverlap = noverlap,
    nfft = nfft,
    fs = fs,
    x_len = length(x),
    method = "Welch"
  ), class = c("ravetools-pwelch", "pwelch"))

  if( plot ) {
    if(!is.null(log)){
      log <- match.arg(log)
    }
    plot(res, col = col, xlim = xlim, ylim = ylim, main = main,
         add = plot >= 2, log = log, ...)
    return(invisible(res))
  }
  return(res)
}

#' @rdname pwelch
#' @export
print.pwelch <- function(x, ...){
  cat(paste0(
    "Welch Periodogram:\n",
    sprintf("  time points: %d\n", x$x_len),
    sprintf("  sample rate: %.2f\n", x$fs),
    sprintf("  window size: %d\n", length(x$window)),
    sprintf("  window overlaps: %d\n", x$noverlap),
    sprintf("  filter count: %d\n", x$nfft)
  ))
  invisible(x)
}

#' @rdname pwelch
#' @export
plot.pwelch <- function(x, log = c("xy", "x", "y", ""), type = 'l', add = FALSE, col = 1, cex = 1, cex.main = cex, cex.sub = cex, cex.lab = cex * 0.8, cex.axis = cex * 0.7, las = 1, main = 'Welch periodogram', xlab, ylab, xlim = NULL, ylim = NULL, ...) {
  if(!is.null(log)){
    log <- match.arg(log)
  } else {
    log <- ''
  }
  freq <- x$freq
  spec <- x$spec

  if(!length(xlim)){
    xlim <- range(freq)
  }

  switch (
    log,
    "xy" = {
      xlab %?<-% 'Log10(Frequency)'
      ylab %?<-% 'Power (dB)'
      freq <- log10(freq)
      spec <- log10(spec) * 10
      xlabel <- pretty(10^xlim)
      xat <- xlabel
      xat[xat <= 0 ] <- 1
      xat <- log10(xat)
      xlim <- range(xat)
    },
    "x" = {
      xlab %?<-% 'Log10(Frequency)'
      ylab %?<-% 'Power'
      freq <- log10(freq)
      xlabel <- pretty(10^xlim)
      xat <- xlabel
      xat[xat <= 0 ] <- 1
      xat <- log10(xat)
      xlim <- range(xat)
    },
    "y" = {
      xlab %?<-% 'Frequency'
      ylab %?<-% 'Power (dB)'
      spec <- log10(spec) * 10
      xlabel <- pretty(xlim)
      xat <- xlabel
    },
    {
      xlab %?<-% 'Frequency'
      ylab %?<-% 'Power'
      xlabel <- pretty(xlim)
      xat <- xlabel
    }
  )
  if(!length(ylim)){
    ylim <- range(spec)
  }

  if(add){
    graphics::points(freq, spec, type = type, col = col, ...)
  } else {

    graphics::plot(
      freq, spec, type = type, col = col, xlab = xlab, ylab = ylab,
      xlim = xlim, ylim = ylim, main = main, las = las,
      cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main,
      cex.sub = cex.sub, axes = FALSE, ...)
    graphics::axis(1, at = xat, labels = xlabel)
    graphics::axis(2, at = pretty(ylim), las = 1)
  }
  invisible()
}



#' @rdname pwelch
#' @export
mv_pwelch <- function(x, margin, fs, nfft){
  xlen <- length(x) / dim(x)[[margin]]
  window_len <- xlen
  window <- hanning(xlen)
  if(missing(nfft)){
    nfft <- 2^ceiling(log2(xlen))
  }
  re <- apply(x, margin, function(s){
    a <- detrend_naive(s)
    postpad(a$Y * window, nfft)
  })
  re <- Mod(mvfftw_r2c(re, inplace = TRUE))^2

  NN <- floor((nfft + 1)/2)
  spec <- rowMeans(re) / (window_len / 2)^2
  spec <- spec[seq_len(NN)]
  freq <- seq(1, fs / 2, length.out = NN)

  res <- structure(list(
    freq = freq,
    spec = spec,
    window = window,
    noverlap = NA,
    nfft = nfft,
    fs = fs,
    x_len = xlen,
    method = "Welch"
  ), class = c("ravetools-pwelch", "pwelch"))
}
