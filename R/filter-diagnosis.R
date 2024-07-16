# wrap phase angles
unwrap <- function (x, tol = pi) {
  if (is.vector(x)) {
    x <- as.matrix(x, ncol = 1)
    vec <- TRUE
  } else {
    vec <- FALSE
  }
  nr <- nrow(x)
  nc <- ncol(x)
  if (!is.numeric(x)) {
    stop("x must be a numeric matrix or vector")
  }
  tol <- abs(tol)
  y <- x
  if (nr > 1) {
    rng <- 2 * pi
    for (col in seq_len(nc)) {
      valid <- is.finite(x[, col])
      d <- diff(x[valid, col])
      p <- round(abs(d)/rng) * rng * (((d > tol) > 0) - ((d < -tol) > 0))
      r <- cumsum(p)
      y[valid, col] <- x[valid, col] - c(0, r)
    }
  }
  if (vec) {
    y <- as.vector(y)
  }
  y
}

#' Frequency response of digital filter
#' @description
#' Compute the z-plane frequency response of an \code{ARMA} model.
#' @param b the moving-average coefficients of an \code{ARMA} model
#' @param a the auto-regressive coefficients of an \code{ARMA} filter; default
#' is \code{1}
#' @param fs sampling frequency in \code{Hz}
#' @param n number of points at which to evaluate the frequency response;
#' default is \code{512}
#' @param whole whether to evaluate beyond \code{Nyquist} frequency; default
#' is false
#' @param ... ignored
#' @returns A list of frequencies and corresponding responses in complex vector
#' @export
freqz2 <- function(b, a = 1, fs = 2 * pi, n = 512, whole = FALSE, ...) {
  freqz_internal <- asNamespace("gsignal")[["freqz.default"]]
  hw <- freqz_internal(unclass(b), a = unclass(a), fs = fs, n = n, whole = whole, ...)

  hw$f <- hw$w
  hw$srate <- fs
  class(hw) <- c("ravetools-freqz2", "ravetools-printable")
  hw
}

#' @export
`summary.ravetools-freqz2` <- function(object, ..., cutoffs = c(-3, -6, -12)) {
  nm <- deparse(substitute(object))
  h <- object$h
  w <- object$w
  rw <- c(0, object$srate / 2)
  mag <- 20 * log10(abs(h))
  phase <- unwrap(Arg(h))
  rp <- range(phase, na.rm = TRUE)

  isfinite <- is.finite(mag)
  if(any(isfinite)) {
    idx <- which.max(mag[isfinite])
    idx <- which(isfinite)[[idx[[1]]]]
    peak <- list(which = idx, freq = object$w[[idx]], mag = mag[[idx]])
  } else {
    peak <- list(which = NA, freq = NA, mag = NA)
  }

  cutoffs <- as.numeric(cutoffs)
  cutoffs <- cutoffs[is.finite(cutoffs)]

  cutoffs <- lapply(cutoffs, function(cutoff) {
    # cutoff <- -3
    idx <- which(diff(as.integer(isfinite & mag > cutoff)) != 0)
    if(!length(idx)) { return(NULL) }

    m1 <- mag[idx]
    m2 <- mag[idx + 1]

    a <- ifelse(m1 == m2, 1, (cutoff - m2) / (m1 - m2))
    freq <- a * w[idx] + (1-a) * w[idx + 1]


    list(
      cutoff = cutoff,
      frequency = freq
    )
  })

  cutoffs <- cutoffs[!vapply(cutoffs, is.null, FALSE)]

  structure(list(name = nm, unit = object$u,
                 frequency_range = rw, phase_range = rp,
                 peak = peak, cutoffs = cutoffs),
            class = c("summary.ravetools-freqz2", "ravetools-printable", "list"))
}

#' @export
`format.ravetools-freqz2` <- function(x, ..., name = deparse(substitute(x))) {
  force(name)
  res <- summary(x, ...)
  res$name <- name
  format(res)
}

#' @export
`format.summary.ravetools-freqz2` <- function(x, ...) {

  s <- c(
    sprintf("<RAVE filter freqz summary> %s", paste(x$name, collapse = "")),
    sprintf("Frequency ranges: %.5g to %.5g %s", x$frequency_range[[1]], x$frequency_range[[2]], x$unit)
  )
  if(!is.na(x$peak[[1]])) {
    s <- c(
      s,
      sprintf("Peak magnitude: %.2f dB at frequency %.5g %s",
              x$peak$mag, x$peak$freq, x$unit)
    )
  }

  s <- c(s, unlist(lapply(x$cutoffs, function(item) {
    # item <- x$cutoffs[[1]]
    if(length(item$frequency) > 4) {
      re <- c(sprintf("%.4g %s", item$frequency[1:4], x$unit), "...")
    } else {
      re <- sprintf("%.4g %s", item$frequency, x$unit)
    }

    sprintf("  %4.4g dB cutoff at %s", item$cutoff, paste(re, collapse = ", "))
  })))

  rpd <- x$phase_range * 360/(2 * pi)
  s <- c(
    s,
    sprintf("Phase ranging from %.4g to %.4g rad (%.3f to %.3f degrees)",
            x$phase_range[[1]], x$phase_range[[2]], rpd[[1]], rpd[[2]])
  )
  s
}

#' @export
`plot.ravetools-freqz2` <- function(
    x, ..., cutoffs = c(-3, -6, -12), vlines = NULL,
    ylim = NULL, xlim = NULL, color_palette = NULL,
    add = FALSE, which = c(1,2,3), draw_legend = TRUE
) {

  # DIPSAUS DEBUG START
  # list2env(list(fs = 500, n = 6000, whole = FALSE), .GlobalEnv)
  # b <- ravetools::fir1(160, c(1 / 250, 70 / 250), "pass")$b
  # a <- 1
  # xlim = NULL
  # vlines <- NULL
  # ylim <- NULL
  # cutoffs = c(-3, -6, -12)
  # color_palette <- NULL
  # which = c(1,2,3)
  # add = FALSE
  # x <- freqz2(b = b, a = a, fs = fs, n = n, whole = whole)
  # draw_legend <- TRUE


  which <- which[which %in% c(1, 2, 3)]
  if(!length(which)) {
    which <- c(1, 2, 3)
  }

  mag <- 20 * log10(abs(x$h))
  argh <- Arg(x$h)

  isfinite <- is.finite(mag)
  if(any(isfinite)) {
    idx <- which.max(mag[isfinite])
    idx <- which(isfinite)[[idx[[1]]]]
    maxmag <- mag[[idx]]
  } else {
    idx <- NULL
    maxmag <- 1
  }

  phase <- unwrap(argh)

  w <- x$w
  if(!length(cutoffs)) { cutoffs <- -3 }
  smry <- `summary.ravetools-freqz2`(x, cutoffs = cutoffs)
  smry_str <- format(smry)

  if(length(smry$cutoffs)) {
    if(!length(color_palette)) {
      color_palette <- grDevices::adjustcolor(seq_along(smry$cutoffs) + 1)
    }
    if(length(color_palette) < length(smry$cutoffs)) {
      color_palette <- c(grDevices::adjustcolor(color_palette), grDevices::adjustcolor(seq_along(smry$cutoffs) + 1))
    }
  }


  legend_text <- NULL
  legend_col <- NULL
  legend_lty <- NULL
  if(length(vlines)) {
    legend_text <- "Filter frequency"
    legend_col <- graphics::par("fg")
    legend_lty <- 3
  }

  if(length(smry$cutoffs)) {
    legend_text <- c(legend_text, sprintf("Cutoff %.2f %s", cutoffs, x$u))
    legend_col <- c(legend_col, color_palette[seq_along(cutoffs)])
    legend_lty <- c(legend_lty, rep(2, length(cutoffs)))
    has_legend <- TRUE
  } else {
    has_legend <- FALSE
    draw_legend <- FALSE
  }

  if(!add) {
    if( !isFALSE(draw_legend) ) {
      graphics::layout(
        mat = matrix(c(seq_along(which) + 1, 1), ncol = 1),
        heights = c(rep(1, length(which)), graphics::lcm(3))
      )
    } else {
      graphics::layout(
        mat = matrix(c(seq_along(which)), ncol = 1),
        heights = c(rep(1, length(which)), graphics::lcm(3))
      )
    }
    oldpar <- graphics::par(c("mfrow", "mar"))
    # mfrow = c(length(which), 1), mar = c(4, 4, 1.5, 1))
    on.exit({
      do.call(graphics::par, oldpar)
    })
  }

  # Legend and summary text
  if( !isFALSE(draw_legend) ) {
    if(!identical(draw_legend, "add")) {
      graphics::par(mar = c(0.1,0.1,0.1,0.1))
      plot(c(0,1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
    }

    graphics::legend(x = 0, y = 0.5, xjust = 0, yjust = 0.5,
                     smry_str[-1],
                     bty = "n")

    if( has_legend ) {
      graphics::legend(x = 0.65, y = 0.5, xjust = 0, yjust = 0.5,
                       legend_text,
                       col = legend_col,
                       lty = legend_lty,
                       bty = "n")
    }

  }

  graphics::par(mar = c(2.1, 3.6, 1.5, 1))

  if(!length(xlim)) {
    xlim <- range(pretty(w))
  }
  if(!length(ylim)) {
    ylim <- range(pretty(mag))
  }

  if(1 %in% which) {
    ylim2 <- c(min(cutoffs) - 1, maxmag)
    graphics::plot(range(pretty(w)), ylim2, ...,
                   xlim = xlim, type = "n", xlab = "", ylab = "", main = "", axes = FALSE)
    graphics::lines(w, mag)
    graphics::axis(side = 1, at = pretty(xlim))
    graphics::axis(side = 2, at = pretty(ylim2), las = 1)
    lapply(seq_along(smry$cutoffs), function(ii) {
      # item <- smry$cutoffs[[1]]
      item <- smry$cutoffs[[ii]]
      if(length(item$frequency) > 4) {
        item$frequency <- item$frequency[1:4]
      }
      graphics::abline(h = item$cutoff, v = item$frequency, lty = 2, col = color_palette[[ii]])
    })
    if(length(vlines)) {
      graphics::abline(v = vlines, lty = 3)
    }
    graphics::title("Pass band (dB)")
  }

  if(2 %in% which) {
    graphics::plot(range(pretty(w)), range(mag, na.rm = TRUE), ylim = ylim, ...,
                   type = "n", xlab = "", ylab = "", main = "", axes = FALSE)
    graphics::lines(w, mag)
    graphics::axis(side = 1, at = pretty(w))
    graphics::axis(side = 2, at = pretty(ylim), las = 1)

    lapply(seq_along(smry$cutoffs), function(ii) {
      # item <- smry$cutoffs[[1]]
      item <- smry$cutoffs[[ii]]
      if(length(item$frequency) > 4) {
        item$frequency <- item$frequency[1:4]
      }
      graphics::abline(h = item$cutoff, v = item$frequency, lty = 2, col = color_palette[[ii]])
    })
    if(length(vlines)) {
      graphics::abline(v = vlines, lty = 3)
    }
    graphics::title("Stop band (dB)")
  }

  if(3 %in% which) {
    d <- phase * 360/(2 * pi)
    xlim <- range(xlim, na.rm = TRUE)
    if(length(xlim)) {
      sel <- w >= xlim[[1]] & w < xlim[[2]]
      ylim <- range(d[sel], na.rm = TRUE)
    }
    xlim <- range(pretty(xlim))
    ylim <- range(pretty(ylim))

    graphics::plot(xlim, ylim, ...,
                   type = "n", xlab = "", ylab = "", main = "", axes = FALSE)
    graphics::axis(side = 1, at = pretty(xlim))
    graphics::axis(side = 2, at = pretty(ylim), las = 1)
    graphics::lines(w, d)
    lapply(seq_along(smry$cutoffs), function(ii) {
      # item <- smry$cutoffs[[1]]
      item <- smry$cutoffs[[ii]]
      if(length(item$frequency) > 4) {
        item$frequency <- item$frequency[1:4]
      }
      graphics::abline(v = item$frequency, lty = 2, col = color_palette[[ii]])
    })
    if(length(vlines)) {
      graphics::abline(v = vlines, lty = 3)
    }
    graphics::title("Phase (degrees)")
    graphics::mtext("Frequency", side = 1, outer = FALSE, line = 1.6, cex = 0.8)
  }
  invisible(smry)
}

sample_signal <- function(n) {
  a0 <- c(0, 5, 30, 1, 0, -30, 110, -100, 0, 14, 25, 2, 0, 0, 0)
  d0 <- c(0, 2, 6, 10, 13, 14, 16, 18, 19, 28, 30, 34, 36, 39, 44)
  a <- a0 / max(a0)
  d <- round(d0 * n / d0[15])
  d[15] <- n

  x <- list()
  for(i in seq_len(14)) {
    m <- seq(d[i], d[i + 1] - 1)
    slope <- (a[i + 1] - a[i]) / (d[i + 1] - d[i])
    x[[length(x) + 1]] <- a[i] + slope * (m - d[i])
  }
  unlist(x)
}

#' Diagnose digital filter
#' @description
#' Generate frequency response plot with sample-data simulation
#' @param b the moving-average coefficients of an \code{ARMA} model
#' @param a the auto-regressive coefficients of an \code{ARMA} filter; default
#' is \code{1}
#' @param fs sampling frequency in \code{Hz}
#' @param n number of points at which to evaluate the frequency response;
#' default is \code{512}
#' @param whole whether to evaluate beyond \code{Nyquist} frequency; default
#' is false
#' @param sample sample signal of length \code{n} for simulation
#' @param vlines additional vertical lines (frequencies) to plot
#' @param xlim frequency limit of frequency response plot; default is
#' \code{"auto"}, can be \code{"full"} or a numeric of length 2
#' @param cutoffs cutoff decibel powers to draw on the frequency plot, also used
#' to calculate the frequency limit when \code{xlim} is \code{"auto"}
#' @returns Nothing
#'
#' @examples
#'
#'
#'
#' library(ravetools)
#'
#' # sample rate
#' srate <- 500
#'
#' # signal length
#' npts <- 1000
#'
#' # band-pass
#' bpass <- c(1, 50)
#'
#' # Nyquist
#' fn <- srate / 2
#' w <- bpass / fn
#'
#' # ---- FIR filter ------------------------------------------------
#' order <- 160
#'
#' # FIR1 is MA filter, a = 1
#' filter <- fir1(order, w, "pass")
#'
#' diagnose_filter(
#'   b = filter$b, a = filter$a, n = npts,
#'   fs = srate, vlines = bpass
#' )
#'
#' # ---- Butter filter --------------------------------------------
#' filter <- butter(3, w, "pass")
#'
#' diagnose_filter(
#'   b = filter$b, a = filter$a, n = npts,
#'   fs = srate, vlines = bpass
#' )
#'
#'
#'
#' @export
diagnose_filter <- function(b, a, fs, n = 512, whole = FALSE,
                            sample = stats::rnorm(n, mean = sample_signal(n), sd = 0.2),
                            vlines = NULL, xlim = "auto",
                            cutoffs = c(-3, -6, -12)) {

  # DIPSAUS DEBUG START
  # list2env(list(fs = 500, n = 6000, whole = FALSE), .GlobalEnv)
  # b <- ravetools::fir1(160, c(1 / 250, 70 / 250), "pass")$b
  # a <- 1
  # sample = sample_signal(n)
  # xlim = "auto"
  # vlines <- c(0, 1, 70)
  # cutoffs = c(-3, -6, -12)

  hw <- freqz2(b = b, a = a, fs = fs, n = n, whole = whole)
  force(sample)
  time <- seq_along(sample) / fs
  # res1 <- filter_signal(b = b, a = a, x = sample)[[1]]
  res2 <- filtfilt(unclass(b), a = a, x = sample)

  h <- hw$h
  w <- hw$w
  mag <- 20 * log10(abs(h))

  cutoffs <- cutoffs[is.finite(cutoffs)]
  if(length(xlim) != 2) {
    if(identical(xlim, "auto")) {
      xlim <- NULL
      if(length(cutoffs)) {
        rg <- range(c(w[mag >= min(cutoffs)], vlines), na.rm = TRUE)
        if(all(is.finite(rg))) {
          rg[[1]] <- floor(rg[[1]])
          rg[[2]] <- ceiling(rg[[2]])
          xlim <- rg
        }
      }
    } else {
      xlim <- NULL
    }
  }

  oldpar <- graphics::par(c("mfrow", "mar", "mgp", "tck", "xaxs", "yaxs", "cex.lab"))
  on.exit({
    do.call(graphics::par, oldpar)
  })

  graphics::layout(
    mat = matrix(
      byrow = FALSE,
      ncol = 2,
      c(
        2,3,4,1,
        5,6,7,1
      )
    ), widths = c(3, 2), heights = c(1, 1, 1, graphics::lcm(3))
  )

  cex <- 1

  # mar = c(2.6, 3.8, 2.1, 0.6) * (0.5 + cex / 2), mgp = cex * c(2, 0.5, 0), tck = -0.02 * cex, xaxs =
  # "i", yaxs = "i",
  graphics::par(
    mar = c(0.1,3.6,0.1,1),
    mgp = cex * c(2, 0.5, 0),
    tck = -0.02 * cex,
    xaxs = "i", yaxs = "i",
    cex.lab = 0.8
  )
  plot(c(0,1.1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "", main = "")

  graphics::legend(x = 0.9, y = 0.5, xjust = 0, yjust = 0.5, lty = 1, bty = "n",
         legend = c("Sample", "Filtered"), col = c("gray", "orange"))
  plot(hw, xlim = xlim, add = TRUE, vlines = vlines, cutoffs = cutoffs, draw_legend = "add")

  graphics::par(
    mar = c(2.1, 3.6, 1.5, 1),
    mgp = cex * c(2, 0.5, 0),
    tck = -0.02 * cex,
    xaxs = "i",
    yaxs = "i",
    cex.lab = 0.8
  )
  ylim <- range(pretty(c(sample, sample)))
  plot(range(pretty(time)), ylim, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
  graphics::title("Simulation")

  graphics::lines(time, sample, col = "gray")
  graphics::lines(time, res2, col = "orange")

  graphics::axis(side = 1, at = pretty(time))
  graphics::axis(side = 2, at = pretty(ylim), las = 1)

  graphics::mtext(text = "Time", side = 1, line = 1.6, cex = 0.8)

  # lines(time, sample, lwd = 2)

  pn <- pwelch(sample, fs = fs, window = fs * 2, noverlap = fs)

  # p0 <- pwelch(sample, fs = fs, window = fs * 2, noverlap = fs)

  pr <- pwelch(res2, fs = fs, window = fs * 2, noverlap = fs)

  specs <- pn$spec
  specs <- specs[is.finite(specs) & specs > 0]
  if(length(specs) <= 1) {
    plim <- c(-100, 0)
  } else {
    plim <- 10 * log10( range(specs) )
    plim[[1]] <- plim[[1]] + min(-3, cutoffs, na.rm = TRUE)
  }

  for(plog in c("y", "xy")) {
    plot(pn, add = FALSE, col = "gray", log = plog, mar = c(2.1, 3.6, 1.5, 1), yline = 1.6, grid = FALSE, ylim = plim, ylab = "")
    # plot(p0, add = TRUE, col = "black", log = plog)
    plot(pr, add = TRUE, col = "orange", log = plog)

    if(length(vlines)) {
      if( grepl("x", plog) ) {
        log_vlines <- log10(vlines[vlines > 0])
        if(0 %in% vlines) {
          log_vlines <- c(log_vlines, 0)
        }
        graphics::abline(v = log_vlines, lty = 3)
      } else {
        graphics::abline(v = vlines, lty = 3)
      }
    }
  }

  invisible()

}

