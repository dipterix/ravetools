# Float-precision Morlet wavelet implementation.
#
# Public entry point is `morlet_wavelet()` (defined in R/wavelet.R), which
# dispatches here when `precision = "float"`.

wavelet_kernels2_float <- function(freqs, srate, wave_num,
                             data_length, signature = NULL) {
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
  path <- file.path(root_dir, sprintf("wavelet-%s", digest))


  arr_dim <- c(data_length, length(kernel_info$kernels))
  tryCatch({
    return(filearray::filearray_checkload(
      filebase = path, mode = "readonly", symlink_ok = FALSE,
      freqs = freqs, srate = srate, wave_num = wave_num, data_length = data_length,
      arr_dim = arr_dim, rave_data_type = "rave-wavelet-kernels-float"
    ))
  }, error = function(e) {

    if (getOption("ravetools.debug", FALSE)) {
      print(path)
      warning(e)
    }

    if (dir.exists(path)) {
      unlink(path, recursive = TRUE)
    }

  })

  arr <- filearray::filearray_create(
    filebase = path,
    dimension = arr_dim,
    type = "complex",
    partition_size = 1
  )
  arr$.mode <- "readwrite"

  tmp <- complex(data_length)
  lapply(seq_along(kernel_info$kernels), function(ii) {
    tmp_wavelet <- kernel_info$kernels[[ii]]
    w_l <- length(tmp_wavelet)
    n_pre  <- ceiling(data_length / 2) - floor(w_l / 2)
    n_post <- data_length - n_pre - w_l
    x <- c(rep(0, n_pre), tmp_wavelet, rep(0, n_post))
    fftw_c2c(data = x, inverse = 0L, ret = tmp)
    conjugate(tmp)
    arr[, ii] <- tmp
    NULL
  })


  arr$.header$freqs <- freqs
  arr$.header$srate <- srate
  arr$.header$wave_num <- wave_num
  arr$.header$data_length <- data_length
  arr$.header$arr_dim <- arr_dim
  arr$.header$freqs <- freqs
  arr$.header$rave_data_type <- "rave-wavelet-kernels-float"
  arr$.save_header()

  arr$.mode <- "readonly"

  return(arr)
}

morlet_wavelet_float <- function(data, freqs, srate, wave_num,
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

  # detrend on the full signal preserves legacy semantics regardless of path
  if (trend != "none") {
    data <- as.vector(detrend(data, trend = trend, ...))
  }

  # In segmented mode we need the segmentation plan up-front (it both
  # validates `segment_length` against the longest wavelet kernel and supplies
  # row indices used by the per-frequency fmap closure below).
  if (segmented) {
    kern_info <- wavelet_kernels(freqs = freqs, srate = srate, wave_num = wave_num)
    max_kernel_len <- max(vapply(kern_info$kernels, length, integer(1L)))
    plan <- wavelet_segment_plan(d_l, segment_length_norm, max_kernel_len)
    kern_data_length <- segment_length_norm
  } else {
    plan <- NULL
    kern_data_length <- d_l
  }

  # Kernel cache differentiates by `data_length`, so single-shot and segmented
  # builds get distinct kernel files. Output cache deliberately omits
  # `segment_length` since both paths yield (up to fp) identical coefficients.
  fft_waves <- wavelet_kernels2_float(freqs, srate, wave_num, kern_data_length,
                                      signature = signature)

  out_path <- tempfile2()
  output <- tryCatch({
    filearray::filearray_checkload(
      filebase = out_path, symlink_ok = FALSE,
      freqs = freqs, srate = srate, wave_num = wave_num,
      data_digest = data_digest, trend = trend,
      more_args = more_args,
      rave_data_type = "rave-wavelet-coefficients",
      precision = "float"
    )
  }, error = function(e) {
    if (getOption("ravetools.debug", FALSE)) {
      print(out_path)
      warning(e)
    }
    if (dir.exists(out_path)) { unlink(out_path, recursive = TRUE) }
    NULL
  })
  if (inherits(output, "FileArray")) { return(output) }

  output <- filearray::filearray_create(
    filebase = out_path, dimension = c(d_l, f_l),
    type = "complex", partition_size = 1
  )
  output$.mode <- "readwrite"

  wave_len <- nrow(fft_waves)
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

    # Per frequency: batch-convolve the kernel against all segment columns,
    # then gather the fft-shifted interior of each segment into a length-d_l
    # vector that fmap writes to the corresponding output column. We index
    # `out_cplx` through the precomputed shift-permutation `new_ind` instead
    # of materializing a `shifted` matrix (avoids a large per-frequency alloc).
    
    # Creating buffer to avoid repeated allocations in the fmap loop
    buffer <- array(0 + 0i, dim = c(plan$segment_length, n_seg))
    fmap_fun <- function(input) {
      # out_cplx <- mvfftw_c2c(kern * fft_seg, inverse = 1L) / norm_factor
      # Hence out_cplx = buffer * norm_factor
      mvfftw_c2c(input[[1]] * fft_seg, inverse = 1L, ret = buffer)

      # The math is: shifted = rbind(out_cplx[-ind, , drop = FALSE], out_cplx[ind, , drop = FALSE])
      # or shifted = rbind(buffer[-ind, , drop = FALSE], buffer[ind, , drop = FALSE]) / norm_factor
      unlist(lapply(seq_len(n_seg), function(k) {
        # equivalent to: shifted[local_starts[k]:local_ends[k], k]
        # Division by norm_factor deferred until after subsetting to avoid repeated large allocs.
        buffer[new_ind[local_starts[k]:local_ends[k]], k] / norm_factor
      }))
    }
  } else {
    fft_data <- fftw_r2c(data)
    tmp <- complex(length(fft_data))
    fmap_fun <- function(input) {
      fftw_c2c(data = input[[1]] * fft_data, inverse = 1L, ret = tmp)
      c(tmp[-ind], tmp[ind]) / norm_factor
    }
  }

  filearray::fmap(x = fft_waves, fun = fmap_fun,
                  .y = output, .buffer_count = ncol(output))

  output$.header$freqs <- freqs
  output$.header$srate <- srate
  output$.header$wave_num <- wave_num
  output$.header$data_digest <- data_digest
  output$.header$trend <- trend
  output$.header$more_args <- more_args
  output$.header$rave_data_type <- "rave-wavelet-coefficients"
  output$.header$precision <- "float"
  output$.save_header()

  output$.mode <- "readonly"
  output
}
