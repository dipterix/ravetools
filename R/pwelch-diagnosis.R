diagnose_signal <- function(
  s1, s2 = NULL, sc = NULL, srate, name = '', try_compress = TRUE,
  max_freq = 300, window = ceiling(srate * 2), noverlap = window / 2, std = 3,
  cex = 1.5, lwd = 0.5, flim = NULL, nclass = 100,
  main = 'Channel Inspection', col = c('black', 'red'),
  which = NULL, start_time = 0, boundary = NULL, mar = c(5.2, 5.1, 4.1, 2.1),
  ...){

  # is sc not specified, and srate is too high, compress s1
  if(try_compress && (is.null(sc) || (srate > 200 && length(s1) / srate > 300))){
    sratec = 100
    sc <- s1[round(seq(1, length(s1), by = srate/sratec))]
  }else{
    sc %?<-% s1
    sratec = srate / length(s1) * length(sc)
  }
  max_freq = min(max_freq, floor(srate/ 2))
  xlim = c(0, max_freq)

  # Calculate boundary to draw
  if(is.null(boundary)){
    boundary = std* stats::sd(s1)
  }
  ylim = max(abs(s1), boundary)

  # Grid layout
  if(length(which) == 0){
    # grid::grid.newpage()
    lay <- rbind(c(1,1,1), c(2,3,4))
    graphics::par(mar = mar)
    graphics::layout(mat = lay)
    # mai = graphics::par('mai');
    # on.exit({graphics::par(mai = mai)}, add = T)
    # graphics::par(mai = c(1.1, 0.8 ,0.4, 0.25))
  }

  # First plot: plot sc directly with col[1]
  if(length(which) == 0 || 1 %in% which){
    graphics::plot(start_time + (seq_along(sc) / sratec), sc, xlab = 'Time (seconds)', ylab = 'Voltage',
                   main = main, lwd = lwd,
                   type = 'l', ylim = c(-ylim-1, ylim+1), yaxt="n", col = col[1],
                   cex.axis = cex * 0.7, cex.lab = cex *0.8, cex.main = cex, cex.sub = cex, ...)
    graphics::abline(h = c(-1,1) * boundary, col = 'red')
    ticks<-c(-ylim, -boundary,0,boundary, ylim)
    graphics::axis(2,at=ticks,labels=round(ticks), las = 1,
                   cex.axis = cex*0.7, cex.lab = cex *0.8, cex.main = cex, cex.sub = cex)
  }

  # plot 2, 3 too slow, need to be faster - pwelch periodogram
  if(length(which) == 0 || 2 %in% which){
    if(!is.null(s2)){
      pwelch(s2, fs = srate, window = window,
             noverlap = noverlap, plot = 1, col = col[2], cex = cex, ylim = flim,
             log = 'y', xlim = xlim)
      pwelch(s1, fs = srate, window = window, noverlap = noverlap, cex = cex, ylim = flim,
             plot = 2, col = col[1], log = 'y', xlim = xlim)
      graphics::legend('topright', sprintf('%s %s', c('Before', 'After'), name), col = rev(col), lty = 1, cex = cex * 0.7)
    }else{
      pwelch(s1, fs = srate, window = window,
             noverlap = noverlap, plot = 1, col = col[1], cex = cex, ylim = flim,
             log = 'y', xlim = xlim)
    }
  }


  if(length(which) == 0 || 3 %in% which){
    log_xlim = log10(sapply(xlim, max, 1))
    if(!is.null(s2)){
      pwelch(s2, fs = srate, window = window,
             noverlap = noverlap, plot = 1, col = col[2], cex = cex, ylim = flim,
             log = 'xy', xlim = log_xlim)
      pwelch(s1, fs = srate, window = window, noverlap = noverlap, cex = cex, ylim = flim,
             plot = 2, col = col[1], log = 'xy', xlim = log_xlim)
      graphics::legend('topright', paste0(c('Before ', 'After '), name), col = rev(col), lty = 1, cex = cex * 0.8)
    }else{
      pwelch(s1, fs = srate, window = window,
             noverlap = noverlap, plot = 1, col = col[1], cex = cex, ylim = flim,
             log = 'xy', xlim = log_xlim)
    }
  }


  if(length(which) == 0 || 4 %in% which){
    # Plot 4:
    graphics::hist(s1, nclass = nclass,
                   xlab = 'Signal Voltage Histogram', main = paste0('Histogram ', name),
                   cex.axis = cex * 0.7, cex.lab = cex*0.8, cex.main = cex, cex.sub = cex)
  }

  return(list(
    ylim = ylim,
    boundary = boundary
  ))
}