#' High-performance 'Morlet'-wavellet on 'RAVE' subject
#' @description This function is intended to be called within builtin
#' pipelines. Please check \url{https://rave.wiki} to see details.
#' @param subject character of subject ID or \code{\link[raveio]{RAVESubject}}
#' instance
#' @param electrode length of one, which electrode should be applied
#' @param frequencies,kernel_cycles,demean see \code{\link{morlet_wavelet}}
#' @param downsample_rate down-sample rate after wavelet
#' @param signal_source which 'HDF5' name should be used. Current default is
#' \code{"/notch/{block}"} for 'Notch' filtered signals
#' @return Current time-stamp. The wavelet power and phase will be saved to
#' 'RAVE' structure automatically
#' @export
wavelet_rave_subject <- function(
  subject, electrode, frequencies, kernel_cycles,
  downsample_rate = 100, demean = TRUE,
  signal_source = "/notch/{block}"){
  stopifnot2(length(electrode) == 1, msg = "`wavelet_rave_subject` only wavelet one electrode each time")
  subject <- raveio::as_rave_subject(subject, strict = TRUE)
  preprocess_instance <- subject$preprocess_settings
  blocks = preprocess_instance$blocks
  sel <- preprocess_instance$electrodes %in% electrode
  if(!any(sel)){
    stop("Cannot find electrode ", electrode, " in the subject.")
  }
  sample_rate <- preprocess_instance$sample_rates[sel]
  preprocess_path <- preprocess_instance$subject$preprocess_path
  data_path <- preprocess_instance$subject$data_path
  cache_path <- preprocess_instance$subject$cache_path
  cache_nparts <- max(subject$electrodes)
  raveio::dir_create2(file.path(cache_path, "wavelet"))

  e <- electrode

  h5_path <- file.path(preprocess_path, 'voltage', sprintf('electrode_%d.h5', e))
  h5_power <- file.path(data_path, "power", sprintf("%d.h5", e))
  h5_phase <- file.path(data_path, "phase", sprintf("%d.h5", e))
  # h5_voltage <- file.path(data_path, "voltage", sprintf("%d.h5", e))
  raveio::dir_create2(dirname(h5_path))
  raveio::dir_create2(dirname(h5_power))
  raveio::dir_create2(dirname(h5_phase))
  # raveio::dir_create2(dirname(h5_voltage))

  # load all data
  signals <- structure(lapply(blocks, function(block){
    h5_name <- raveio::glue(signal_source)
    raveio::load_h5(h5_path, h5_name, ram = TRUE)
  }), names = blocks)

  srate_volt <- sample_rate
  srate_power <- downsample_rate

  nfreq <- length(frequencies)

  chunk_size <- floor(1024 / nfreq)
  if( chunk_size < 1 ){ chunk_size <- 1 }
  chunk_size <- c(nfreq, chunk_size)

  lapply(blocks, function(block){
    result <- ravetools::morlet_wavelet(
      data = signals[[block]],
      freqs = frequencies,
      srate = srate_volt,
      wave_num = kernel_cycles,
      demean = demean
    )
    len <- nrow(result)
    npoints <- ceiling(len * srate_power / srate_volt)
    result <- t(result[round(seq(1, len, length.out = npoints)), ])

    # save HDF5 file - power
    raveio::save_h5(
      x = Mod(result) ^ 2,
      file = h5_power,
      name = sprintf("/raw/power/%s", block),
      chunk = chunk_size,
      replace = TRUE,
      level = 7
    )

    # save HDF5 file - phase
    raveio::save_h5(
      x = Arg(result),
      file = h5_phase,
      name = sprintf("/raw/phase/%s", block),
      chunk = chunk_size,
      replace = TRUE,
      level = 7
    )

    # save cache
    cache_tensor_path <- file.path(cache_path, "wavelet", block)
    arr <- NULL
    if(dir.exists(cache_tensor_path)){
      tryCatch({
        arr <- filearray::filearray_load(cache_tensor_path, "readwrite")
        dim <- arr$dimension()
        nparts <- dim[[length(dim)]]
        if( cache_nparts > nparts ){
          # need to expand the array
          meta <- file.path(arr$.filebase, "meta")
          fid <- file(meta, "r+b")
          on.exit({
            try({close(fid)}, silent = TRUE)
          }, after = FALSE, add = TRUE)
          seek(fid, where = 36, origin = "start", rw = "write")
          writeBin(as.double(length(result) * cache_nparts), fid, size = 8L, endian = "little")

          seek(fid, where = 64, origin = "start", rw = "write")
          writeBin(as.double(cache_nparts), fid, size = 8L, endian = "little")
          seek(fid, where = 0, origin = "start", rw = "write")
          close(fid)
          arr <- filearray::filearray_load(cache_tensor_path, "readwrite")
        }
      }, error = function(e){
        unlink(cache_tensor_path, recursive = TRUE)
      })
    }
    if(!dir.exists(cache_tensor_path)){
      arr <- filearray::filearray_create(
        filebase = cache_tensor_path,
        dimension = c(npoints, nfreq, cache_nparts),
        type = "complex",
        partition_size = 1L
      )
    }
    arr[,,e] <- t(result)
    NULL
  })

  # preprocess_instance$data[[electrode]]$has_wavelet <- TRUE
  # preprocess_instance$save()

  Sys.time()
}

