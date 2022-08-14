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
#' @param col,xlim,ylim,main,type,cex,cex.main,cex.sub,cex.lab,cex.axis,las,xlab,ylab,lty,lwd,xaxs,yaxs,mar parameters passed to \code{\link[graphics]{plot.default}}
#' @param se logical or a positive number indicating whether to plot standard
#' error of mean; default is false. If provided with a number, then a multiple
#' of standard error will be drawn. This option is only available when power
#' is in log-scale (decibel unit)
#' @param col.se,alpha.se controls the color and opacity of the standard error
#' @param xticks ticks to show on frequency axis
#' @param xline,yline controls how close the axis labels to the corresponding axes
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
  offset <- seq(1, max(x_len-window_len+1, 1), by = step)

  N <- length(offset)

  re <- sapply(seq_len(N), function(i){
    slice <- x[offset[i] - 1 + seq_len(window_len)]
    slice[is.na(slice)] <- 0
    a <- detrend_naive(slice)
    postpad(a$Y * window, nfft)
    # # HermConj = 0 : without the "Hermitian" redundancy
    # a = fftwtools::fftw_r2c(postpad(a$Y * window, nfft), HermConj = 0)
    # Mod(a)^2
  })

  re <- Mod(mvfftw_r2c(re))^2

  NN <- floor((nfft + 1)/2)

  re <- re[seq_len(NN), , drop = FALSE] / (window_len / 2)^2

  spec <- rowMeans(re)

  # decibel unit so we can calculate sterr of mean
  re_db <- 20 * log10(re)
  spec_db <- rowMeans(re_db)
  if(N > 1) {
    spec_db_se <- apply(re_db, 1, stats::sd) / sqrt(N - 1)
  } else {
    spec_db_se <- rep(N, nrow(re))
  }

  freq <- seq(1, fs / 2, length.out = NN)

  res <- structure(list(
    freq = freq,
    spec = spec,
    spec_db = spec_db,
    spec_db_se = spec_db_se,
    df = N - 1,
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
plot.pwelch <- function(
    x, log = c("xy", "x", "y", ""), se = FALSE, xticks, type = 'l', add = FALSE,
    col = 1, col.se = "orange", alpha.se = 0.5, lty = 1, lwd = 1,
    cex = 1, cex.main = cex, cex.sub = cex, cex.lab = cex * 0.8,
    cex.axis = cex * 0.7, las = 1, main = 'Welch periodogram', xlab, ylab,
    xlim = NULL, ylim = NULL, xaxs ="i", yaxs = "i", xline = 1.5, yline = 2.5,
    mar = c(2.6, 3.5, 1.1, 0.6),
    ...) {
  if(!is.null(log) && !identical(log, "")){
    log <- match.arg(log)
  } else {
    log <- ''
  }
  freq <- x$freq

  if(!length(xlim)){
    xlim <- range(freq)
  } else {
    if(x$fs > 500 && xlim[[1]] == 0 && xlim[[2]] < 2.5) {
      warning("`plot.pwelch`: `xlim` should be the frequency range in native values not log-range. Please check your plot function")
      xlim <- 10^(xlim)
    }
  }

  if(!missing(xticks)) {
    xticks <- xticks[xticks >= min(xlim) & xticks <= max(xlim)]
    xlabel <- c(xticks, xlim)
  } else {
    xlabel <- pretty(xlim)
  }

  if( x$df < 1 || log %in% c("x", "") ) {
    se <- FALSE
    spec_lb <- 0
  }
  if( se ) {
    spec_lb <- x$spec_db - x$spec_db_se * as.numeric(se)
    spec_ub <- x$spec_db + x$spec_db_se * as.numeric(se)
  } else {
    spec_lb <- 0
    spec_ub <- 0
  }

  switch (
    log,
    "xy" = {
      xlab %?<-% 'Log10(Frequency)'
      ylab %?<-% 'Power (dB)'

      xat <- xlabel
      xat[xat <= 0 ] <- min(freq[freq > 0])
      xat <- log10(xat)
      freq <- log10(freq)
      spec <- x$spec_db

      xlim <- range(xat)
      if(xlim[1] < min(freq)) {
        xlim[1] <- min(freq)
      }
    },
    "x" = {
      xlab %?<-% 'Log10(Frequency)'
      ylab %?<-% 'Power'
      xat <- xlabel
      xat[xat <= 0 ] <- min(freq[freq > 0])
      xat <- log10(xat)
      freq <- log10(freq)
      spec <- x$spec

      xlim <- range(xat)
      if(xlim[1] < min(freq)) {
        xlim[1] <- min(freq)
      }
      se <- FALSE
    },
    "y" = {
      xlab %?<-% 'Frequency'
      ylab %?<-% 'Power (dB)'
      spec <- x$spec_db
      xat <- xlabel
    },
    {
      xlab %?<-% 'Frequency'
      ylab %?<-% 'Power'
      spec <- x$spec
      xat <- xlabel
      se <- FALSE
    }
  )
  if(!length(ylim)){
    ylim <- range(pretty(spec))
  }

  if(add){
    graphics::points(freq, spec, type = type, col = col, lty = lty, lwd = lwd, ...)
  } else {

    old_par <- graphics::par(c("mgp", "mar"))
    graphics::par(mgp = c(2, 0.5, 0), mar = mar)
    on.exit({
      do.call(graphics::par, old_par)
    })

    graphics::plot(
      range(freq), range(spec), type = "n", xlab = "", ylab = "",
      xlim = xlim, ylim = ylim, main = main, las = las,
      cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main,
      cex.sub = cex.sub, axes = FALSE,
      xaxs = xaxs, yaxs = yaxs, ...)
    graphics::axis(1, at = xat, labels = xlabel, tck = -0.03, cex.axis = cex.axis)
    graphics::axis(2, at = pretty(ylim), las = 1, tck = -0.03, cex.axis = cex.axis)

    graphics::mtext(side = 2, text = ylab, line = yline, cex = cex.lab)
    graphics::mtext(side = 1, text = xlab, line = xline, cex = cex.lab)

    if(!isFALSE(se)) {
      graphics::polygon(c(freq, rev(freq)), c(spec_lb, rev(spec_ub)),
              density = NA, border = NA,
              col = grDevices::adjustcolor(col.se, alpha.f = alpha.se))
    }
    graphics::lines(x = freq, y = spec, col = col, lty = lty, lwd = lwd, type = type)

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
  re <- Mod(mvfftw_r2c(re))^2

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
