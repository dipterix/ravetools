hertz_to_erb_scale <- function(f) {
  21.4 * log10(1 + 4.37e-3 * f)
}

erb_scale_to_hertz <- function(erb_scale) {
  (10 ^ (erb_scale / 21.4) - 1) / 4.37e-3
}

hertz_to_erb <- function(f) {
  24.7 * (4.37 * f / 1000 + 1)
}

nextpow2 <- function(x) {
  ceiling(log2(x))
}

gammatone_fast <- function(x, sample_rate, center_frequencies, n_bands,
                           use_hilbert = FALSE) {

  # x <- as.vector(ieegio::io_read_mat("~/rave_data/raw_dir/DemoSubject/008/DemoSubjectDatafile008_ch13.mat")$analogTraces)[1:(2000 * 5)]
  # sample_rate <- 2000
  # center_frequencies <- c(20, 1000)
  # n_bands <- 128

  filter_order <- 4
  center_frequencies <- unname(as.double(center_frequencies))
  if(missing(n_bands)) {
    n_bands <- length(center_frequencies)
  } else {
    if(length(center_frequencies) != 2) {
      stop("`gammatone_fast`: `n_bands` is provided but length of `center_frequencies` is not two. Please either leave `n_bands` empty or set `center_frequencies` to be the frequency ranges.")
    }
    center_frequencies <- sort(center_frequencies)
    centr_erbs <- hertz_to_erb_scale(center_frequencies)
    centr_erbs <- seq(centr_erbs[[1]], centr_erbs[[2]], length.out = n_bands)
    center_frequencies <- erb_scale_to_hertz(centr_erbs)
  }

  gammatone_length <- 2 ^ nextpow2(0.128 * sample_rate)

  # Bandwidth factor for gammatone is 1.019
  erb <- hertz_to_erb(center_frequencies)
  b <- 1.019 * erb

  phase <- 0

  two_pi_per_timepoint <- 2 * pi / sample_rate

  # closed form (?) integral of impulse
  # gain <- ((1.019 * b * two_pi_per_timepoint) ^ filter_order) / 6

  time_seq <- (seq_len(gammatone_length) - 1L) / sample_rate

  # for each frequency center calculate IR
  tcs <- NULL
  gammatone_filters <- sapply(center_frequencies, function(f) {
    # f <- center_frequencies[[1]]
    erb <- hertz_to_erb(f)
    b <- 1.019 * erb
    gain <- ((1.019 * b * two_pi_per_timepoint) ^ filter_order) / 6
    # make sure no phase shift
    tc <- (filter_order - 1) / (2 * pi * b)
    tcs <<- c(tcs, tc)
    phase = -2 * pi * f * tc

    gain * sample_rate^3 * time_seq^(filter_order - 1) * exp(-2 * pi * b * time_seq) * cos(2 * pi * f * time_seq + phase)
  })

  # Use fftfilt() but more efficient for large vectors
  l_filter <- gammatone_length
  x_is_mat <- is.matrix(x)
  if( !x_is_mat ) {
    x <- matrix(x, ncol = 1)
  }
  l_x <- nrow(x)
  c_x <- ncol(x)
  N <- 2 ^ nextpow2(l_x + l_filter - 1)


  root_dir <- file.path(tempdir2(check = TRUE), "ravetools")
  digest <- digest::digest(list(
    center_frequencies = center_frequencies,
    filter_order = filter_order,
    N = N,
    sample_rate = sample_rate
  ))
  x_digest <- digest::digest(list(digest::digest(x), digest))

  # if the filters < 256MB
  if( N * n_bands <= 33554432 ) {
    gammatone_filters_fft <- stats::mvfft(postpad(gammatone_filters, N))
  } else {
    if(!dir.exists(root_dir)){
      dir.create(root_dir, showWarnings = FALSE, recursive = TRUE)
    }
    # use filearray to store the signals on disk
    filters_path <- file.path(root_dir, sprintf("gammatone-%s", digest))
    gammatone_filters_fft <- filearray::filearray_load_or_create(
      filebase = filters_path,
      dimension = c(N, n_bands),
      type = "complex",
      symlink_ok = FALSE,
      partition_size = 1L,
      initialize = TRUE,
      ready = TRUE,
      on_missing = function(arr) {
        lapply(seq_len(n_bands), function(ii) {
          arr[, ii] <- fft(postpad(gammatone_filters[, ii], N))
          NULL
        })
        arr$set_header("ready", TRUE)
        arr
      }
    )
  }

  # delay
  delay <- round(tcs * sample_rate)

  if( N * n_bands * c_x <= 33554432 ) {
    # using around 2-3GB RAM so it's fine
    result <- apply(x, 2L, function(x_slice) {
      # print(ii <<- ii + 1)
      # x_slice <- x[,1]
      result <- stats::mvfft(gammatone_filters_fft[drop = FALSE] * fft(postpad(x_slice, N)), inverse = TRUE) / N
      result <- Re(result[seq_len(l_x), , drop = FALSE])
      if(use_hilbert) {
        # apply hilbert transform
        result <- abs(hilbert(result))
      }
      # offset delay
      result <- sapply(seq_len(n_bands), function(jj) {
        res <- result[, jj]
        c(res[-seq_len(delay[[jj]])], res[seq_len(delay[[jj]])])
      })
    })
    dim(result) <- c(l_x, n_bands, c_x)
    result <- aperm(result, c(1, 3, 2))

    # samp <- result[, 1, ]
    # n_plot <- 10000
    # matplot(x = seq_len(n_plot) / sample_rate, y = samp[seq_len(n_plot), ], lty = 1, type = "l")
    # diagnose_channel(s1 = as.vector(x), s2 = rowSums(samp), srate = sample_rate, which = 1, col = c("black", "red"))
    # plot_signals(
    #   t(cbind(x, rowSums(samp))), sample_rate = sample_rate, space = 0, duration = 1
    # )

  } else {
    if(!dir.exists(root_dir)){
      dir.create(root_dir, showWarnings = FALSE, recursive = TRUE)
    }
    result_path <- file.path(root_dir, sprintf("gammatone-results-%s", x_digest))
    result <- filearray::filearray_load_or_create(
      filebase = result_path,
      dimension = c(l_x, c_x, n_bands),
      type = "double",
      symlink_ok = FALSE,
      partition_size = 1L,
      initialize = FALSE,
      ready = FALSE,
      verbose = TRUE,
      on_missing = function(arr) {
        fft_x <- stats::mvfft(postpad(x, N))
        lapply(seq_len(n_bands), function(jj) {
          result <- stats::mvfft(fft_x * gammatone_filters_fft[, jj], inverse = TRUE) / N
          result <- Re(result[seq_len(l_x), , drop = FALSE])
          if(use_hilbert) {
            # apply hilbert transform
            result <- abs(hilbert(result))
          }
          # result is nrow(x) x ncol(x)
          # offset delay
          result <- rbind(
            result[-seq_len(delay[[jj]]), , drop = FALSE],
            result[seq_len(delay[[jj]]), , drop = FALSE]
          )

          arr[, , jj] <- result
          NULL
        })
        arr$set_header("center_frequencies", center_frequencies, save = FALSE)
        arr$set_header("filter_order", as.integer(filter_order), save = FALSE)
        arr$set_header("sample_rate", sample_rate, save = FALSE)
        arr$set_header("ready", TRUE, save = TRUE)
        arr
      }
    )
  }

  result


  # n_plot <- 1:min(10000, length(x))
  # samp <- result[n_plot, 1, ]
  # matplot(x = n_plot / sample_rate, y = samp, lty = 1, type = "l")
  # sp <- ravetools::plot_signals(
  #   as.vector(scale(as.vector(x)[n_plot], scale = FALSE)), sample_rate = sample_rate, space = 1, duration = 1
  # )
  # ravetools::plot_signals(
  #   rowMeans(samp) * 10, sample_rate = sample_rate, space = sp$space, space_mode = "absolute", duration = 1, new_plot = F,
  #   col = "red"
  # )
  # abline(h = sp$space , col = "gray", lty=2)
  # for(ii in 1:128) {
  #   ravetools::plot_signals(
  #     samp[, ii], sample_rate = sample_rate, space = sp$space, space_mode = "absolute", duration = 1, new_plot = F,
  #     col = ii+1
  #   )
  # }

}
