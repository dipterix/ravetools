# Double-precision Morlet wavelet implementation.
#
# Public entry point is `morlet_wavelet()` (defined in R/wavelet.R), which
# dispatches here when `precision = "double"`.

wavelet_kernels2_double <- function(freqs, srate, wave_num, data_length, signature = NULL) {
  freqs <- as.double(freqs)
  srate <- as.double(srate)
  wave_num <- as.double(wave_num)
  data_length <- as.integer(data_length)
  kernel_info <- wavelet_kernels(freqs = freqs, srate = srate, wave_num = wave_num)
  digest <- digest::digest(list(freqs, srate, wave_num, data_length, signature = signature))
  root_dir <- file.path(tempdir2(check = TRUE), "ravetools")
  if (!dir.exists(root_dir)) {
    dir.create(root_dir, showWarnings = FALSE, recursive = TRUE)
  }
  real_path <- file.path(root_dir, sprintf("wavelet-real-%s", digest))
  imag_path <- file.path(root_dir, sprintf("wavelet-imag-%s", digest))


  arr_dim <- c(data_length, length(kernel_info$kernels))
  tryCatch({
    return(list(
      real = filearray::filearray_checkload(
        filebase = real_path, mode = "readonly", symlink_ok = FALSE,
        freqs = freqs, srate = srate, wave_num = wave_num, data_length = data_length,
        arr_dim = arr_dim, rave_data_type = "rave-wavelet-kernels-double-real"
      ),
      imag = filearray::filearray_checkload(
        filebase = imag_path, mode = "readonly", symlink_ok = FALSE,
        freqs = freqs, srate = srate, wave_num = wave_num, data_length = data_length,
        arr_dim = arr_dim, rave_data_type = "rave-wavelet-kernels-double-imag"
      )
    ))
  }, error = function(e) {
    if (getOption("ravetools.debug", FALSE)) {
      print(real_path)
      print(imag_path)
      warning(e)
    }

    if (dir.exists(real_path)) { unlink(real_path, recursive = TRUE) }
    if (dir.exists(imag_path)) { unlink(imag_path, recursive = TRUE) }

  })

  arr_real <- filearray::filearray_create(
    filebase = real_path,
    dimension = arr_dim,
    type = "double",
    partition_size = 1
  )
  arr_imag <- filearray::filearray_create(
    filebase = imag_path,
    dimension = arr_dim,
    type = "double",
    partition_size = 1
  )
  arr_real$.mode <- "readwrite"
  arr_imag$.mode <- "readwrite"

  tmp <- complex(data_length)
  lapply(seq_along(kernel_info$kernels), function(ii) {
    tmp_wavelet <- kernel_info$kernels[[ii]]
    w_l <- length(tmp_wavelet)
    n_pre  <- ceiling(data_length / 2) - floor(w_l / 2)
    n_post <- data_length - n_pre - w_l
    x <- c(rep(0, n_pre), tmp_wavelet, rep(0, n_post))
    fftw_c2c(data = x, inverse = 0L, ret = tmp)
    conjugate(tmp)
    arr_real[, ii] <- Re(tmp)
    arr_imag[, ii] <- Im(tmp)
    NULL
  })

  arr_real$.header$freqs <- freqs
  arr_real$.header$srate <- srate
  arr_real$.header$wave_num <- wave_num
  arr_real$.header$data_length <- data_length
  arr_real$.header$arr_dim <- arr_dim
  arr_real$.header$freqs <- freqs
  arr_real$.header$rave_data_type <- "rave-wavelet-kernels-double-real"
  arr_real$.save_header()

  arr_imag$.header$freqs <- freqs
  arr_imag$.header$srate <- srate
  arr_imag$.header$wave_num <- wave_num
  arr_imag$.header$data_length <- data_length
  arr_imag$.header$arr_dim <- arr_dim
  arr_imag$.header$freqs <- freqs
  arr_imag$.header$rave_data_type <- "rave-wavelet-kernels-double-imag"
  arr_imag$.save_header()


  arr_real$.mode <- "readonly"
  arr_imag$.mode <- "readonly"

  return(list(
    real = arr_real,
    imag = arr_imag
  ))
}


morlet_wavelet_double <- function(data, freqs, srate, wave_num,
                                  trend = c("constant", "linear", "none"),
                                  signature = NULL, segment_length = NULL, ...) {

  trend <- match.arg(trend)
  freqs <- as.double(freqs)
  srate <- as.double(srate)
  wave_num <- as.double(wave_num)
  more_args <- list(...)
  data_digest <- digest::digest(data)

  f_l <- length(freqs)
  d_l <- length(data)

  segment_length_norm <- wavelet_normalize_segment_length(segment_length, d_l)
  segmented <- !is.na(segment_length_norm)

  # detrend on the full signal first; preserves legacy semantics
  if (trend != "none") {
    data <- as.vector(detrend(data, trend = trend, ...))
  }

  if (segmented) {
    kern_info <- wavelet_kernels(freqs = freqs, srate = srate, wave_num = wave_num)
    max_kernel_len <- max(vapply(kern_info$kernels, length, integer(1L)))
    plan <- wavelet_segment_plan(d_l, segment_length_norm, max_kernel_len)
    kern_data_length <- segment_length_norm
  } else {
    plan <- NULL
    kern_data_length <- d_l
  }

  # Kernel cache differentiates by `data_length`; output cache deliberately
  # omits `segment_length` since both paths produce the same coefficients.
  fft_waves <- wavelet_kernels2_double(freqs, srate, wave_num, kern_data_length,
                                       signature = signature)

  real_path <- tempfile2(pattern = "wavelet-double-real-")
  imag_path <- tempfile2(pattern = "wavelet-double-imag-")
  output <- tryCatch({
    list(
      real = filearray::filearray_checkload(
        filebase = real_path, symlink_ok = FALSE,
        freqs = freqs, srate = srate, wave_num = wave_num,
        data_digest = data_digest, trend = trend,
        more_args = more_args,
        rave_data_type = "rave-wavelet-coefficients-real",
        precision = "double"
      ),
      imag = filearray::filearray_checkload(
        filebase = imag_path, symlink_ok = FALSE,
        freqs = freqs, srate = srate, wave_num = wave_num,
        data_digest = data_digest, trend = trend,
        more_args = more_args,
        rave_data_type = "rave-wavelet-coefficients-imag",
        precision = "double"
      )
    )
  }, error = function(e) {
    if (getOption("ravetools.debug", FALSE)) {
      print(real_path)
      print(imag_path)
      warning(e)
    }
    if (dir.exists(real_path)) { unlink(real_path, recursive = TRUE) }
    if (dir.exists(imag_path)) { unlink(imag_path, recursive = TRUE) }
    NULL
  })
  if (length(output) == 2) { return(output) }

  output_real <- filearray::filearray_create(
    filebase = real_path, dimension = c(d_l, f_l),
    type = "double", partition_size = 1
  )
  output_imag <- filearray::filearray_create(
    filebase = imag_path, dimension = c(d_l, f_l),
    type = "double", partition_size = 1
  )
  output_real$.mode <- "readwrite"
  output_imag$.mode <- "readwrite"

  wave_len <- nrow(fft_waves$real)
  ind <- seq_len(ceiling(wave_len / 2))
  norm_factor <- wave_len * sqrt(srate / 2)

  if (segmented) {
    seg_mat <- wavelet_build_segment_matrix(data, plan)
    fft_seg <- mvfftw_r2c(seg_mat, HermConj = 1L)
    n_seg <- plan$n_seg
    local_starts <- plan$local_starts
    local_ends <- plan$local_ends

    new_ind <- seq_len(plan$segment_length)
    new_ind <- c(new_ind[-ind], new_ind[ind])

    # Pre-allocated buffer avoids a full segment_length * n_seg complex alloc
    # per frequency. mvfftw_c2c writes in-place via `ret`; norm_factor division
    # is deferred to the small per-segment subsets.
    buffer <- array(0 + 0i, dim = c(plan$segment_length, n_seg))
    ii <- 1L
    fmap_fun <- function(input) {
      kern <- input[[1]] + 1i * input[[2]]
      mvfftw_c2c(kern * fft_seg, inverse = 1L, ret = buffer)
      col_vec <- unlist(lapply(seq_len(n_seg), function(k) {
        buffer[new_ind[local_starts[k]:local_ends[k]], k] / norm_factor
      }))
      output_real[, ii] <- Re(col_vec)
      output_imag[, ii] <- Im(col_vec)
      ii <<- ii + 1L
      NULL
    }
  } else {
    fft_data <- fftw_r2c(data)
    tmp <- complex(length(fft_data))
    ii <- 1L
    fmap_fun <- function(input) {
      kern <- input[[1]] + 1i * input[[2]]
      fftw_c2c(data = kern * fft_data, inverse = 1L, ret = tmp)
      shifted <- c(tmp[-ind], tmp[ind]) / norm_factor
      output_real[, ii] <- Re(shifted)
      output_imag[, ii] <- Im(shifted)
      ii <<- ii + 1L
      NULL
    }
  }

  filearray::fmap2(
    x = fft_waves,
    fun = fmap_fun,
    .buffer_count = f_l,
    .simplify = FALSE
  )

  output_real$.header$freqs <- freqs
  output_real$.header$srate <- srate
  output_real$.header$wave_num <- wave_num
  output_real$.header$data_digest <- data_digest
  output_real$.header$trend <- trend
  output_real$.header$more_args <- more_args
  output_real$.header$rave_data_type <- "rave-wavelet-coefficients-real"
  output_real$.header$precision <- "double"
  output_real$.save_header()

  output_imag$.header$freqs <- freqs
  output_imag$.header$srate <- srate
  output_imag$.header$wave_num <- wave_num
  output_imag$.header$data_digest <- data_digest
  output_imag$.header$trend <- trend
  output_imag$.header$more_args <- more_args
  output_imag$.header$rave_data_type <- "rave-wavelet-coefficients-imag"
  output_imag$.header$precision <- "double"
  output_imag$.save_header()

  output_real$.mode <- "readonly"
  output_imag$.mode <- "readonly"
  list(
    real = output_real,
    imag = output_imag
  )
}
