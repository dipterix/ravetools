#' Export diagnostic plot (\code{\link{pwelch}}) for 'RAVE' subjects
#' @param subject character, or \code{\link[raveio]{RAVESubject}}
#' @param electrodes electrodes to generate plots
#' @param blocks the session block number
#' @param h5_names the 'HDF5' data names used
#' @param save_dir where to save the 'PDF' files
#' @param width,height,useDingbats,onefile,... passed to \code{\link[grDevices]{pdf}}
#' @param winlens,freq_lims see \code{\link{pwelch}}
#' @param nclass number of bins in the histogram
#' @param cex font size
#' @param col,fore_col,back_col the color to use for different types of signals;
#' by default signals after the 'Notch' filters uses \code{fore_col} color
#' (black), and original trace uses \code{back_col} (grey) color
#' @return Nothing, however, one or multiple 'PDF' files will be generated at
#' \code{save_dir} directory.
#' @examples
#'
#' # Please download RAVE demo subject first
#' if(interactive()){
#'   export_path <- tempfile()
#'   export_diagnose_voltage("demo/DemoSubject", 13:16,
#'                           "008", save_dir = export_path)
#'   # open the following path
#'   normalizePath(export_path)
#' }
#'
#' @export
export_diagnose_voltage <- function(
  subject, electrodes, blocks,
  h5_names = c("/notch/{block}", "/raw/{block}"),
  save_dir = './export', width = 12, height = 7,
  useDingbats = FALSE, onefile = TRUE, winlens, freq_lims,
  nclass = 20, cex = 1,
  col = c(fore_col, back_col), ...,
  fore_col = 'black', back_col = 'grey80'
){
  stopifnot(length(h5_names) <= 2)
  subject <- raveio::as_rave_subject(subject, strict = FALSE)

  blocks <- subject$preprocess_settings$blocks
  if(missing(electrodes)){
    electrodes <- subject$preprocess_settings$electrodes
  }
  sel <- subject$preprocess_settings$electrodes %in% electrodes
  srates <- subject$preprocess_settings$sample_rates[sel]
  winlens %?<-% (2 * srates)
  freq_lims %?<-% (srates / 2)

  raveio::dir_create2(save_dir)

  col <- rep(col, length(h5_names))[seq_along(h5_names)]
  winlens <- rep(winlens, ceiling(length(electrodes) / length(winlens)))[seq_along(electrodes)]
  freq_lims <- rep(freq_lims, ceiling(length(electrodes) / length(freq_lims)))[seq_along(electrodes)]

  # dipsaus::make_forked_clusters(raveio::raveio_getopt("max_worker"))

  dipsaus::lapply_async2(seq_along(electrodes), function(i){
    e <- electrodes[[i]]
    file = file.path(save_dir, sprintf('Welch_Periodogram_%d.pdf', e))
    grDevices::pdf(file = file, width = width, height = height, useDingbats = useDingbats, ..., onefile = onefile)
    tryCatch({
      srate <- srates[[i]]
      winlen <- winlens[[i]]
      freq_lim <- freq_lims[[i]]
      h5_file <- file.path(subject$preprocess_path, "voltage", sprintf("Electrode_%d.h5", e))

      for(block in blocks){

        s <- list()

        for(ii in seq_along(h5_names)){
          try({
            h5_name <- raveio::glue(h5_names[[ii]])
            s[[ii]] <- raveio::load_h5(h5_file, h5_name, ram = TRUE)
          }, silent = TRUE)
        }

        sel <- which(!vapply(s, is.null, FUN.VALUE = FALSE))

        if(setequal(sel, c(1L, 2L))){
          name = 'Notch'
          main = sprintf('Notch Filtered Signal - Block: %s, Electrode: %d', block, e)

          diagnose_signal(
            s1 = s[[1]], s2 = s[[2]], col = col,
            name = name,
            max_freq = freq_lim,
            srate = srate,
            window = winlen,
            noverlap = winlen / 2,
            nclass = nclass,
            cex = cex,
            std = 3,
            lwd = 0.3,
            main = main
          )
        } else if(sel == 1) {
          name = 'Notch'
          main = sprintf('Notch Filtered Signal - Block: %s, Electrode: %d', block, e)

          diagnose_signal(
            s1 = s[[1]],s2 = NULL,  col = col[[1]],
            name = name,
            max_freq = freq_lim,
            srate = srate,
            window = winlen,
            noverlap = winlen / 2,
            nclass = nclass,
            cex = cex,
            std = 3,
            lwd = 0.3,
            main = main
          )

        } else {
          name = 'Raw'
          main = sprintf('Raw Voltage Signal - Block: %s, Electrode: %d', block, e)
          diagnose_signal(
            s1 = s[[2]],s2 = NULL,  col = col[[2]],
            name = name,
            max_freq = freq_lim,
            srate = srate,
            window = winlen,
            noverlap = winlen / 2,
            nclass = nclass,
            cex = cex,
            std = 3,
            lwd = 0.3,
            main = main
          )
        }
      }

    }, error = function(err){
      raveio::catgl('Cannot generate pwelch forthe electrode {e}. Reason: {err$message}', level = 'ERROR')
    })
    grDevices::dev.off()
  }, callback = function(i){
    sprintf("Generating Welch Periodogram|Electrode - %d", electrodes[i])
  }, plan = FALSE)

  invisible()
}