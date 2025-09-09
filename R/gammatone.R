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

#' @title Apply gamma-tone filters to obtain auditory envelopes
#' @param x a numeric vector or matrix; if \code{x} is a matrix, it should
#' be column-major (each column is a sound track)
#' @param sample_rate sampling frequency
#' @param center_frequencies center frequencies at which the envelopes will
#' be derived; can be either a length of two defining the lower and
#' upper bound, and using \code{n_bands} to interpolate automatically, or
#' a length of multiple, with the frequencies specified explicitly
#' @param n_bands number of the center frequencies, can be missing if
#' \code{center_frequencies} is explicit and no interpolation is needed;
#' if specified, then the frequencies will be interpolated using
#' equivalent rectangular bandwidth rate (\code{'ERB'})
#' @param use_hilbert whether to apply 'Hilbert' transform; default is true,
#' which calculates the magnitude; set to false when only the filter is needed
#' @param downsample whether to down-sample the envelopes after the filters;
#' default is \code{NA} (no down-sample).
#' @param downsample_before_hilbert whether the down-sample happens before
#' or after the 'Hilbert' transform so speed up the computation if the signal
#' is too long; only used when \code{downsample} is greater than 1; default
#' is \code{FALSE}. Use with caution, especially when the voice center
#' frequency is close to the 'Nyquist' frequency. However, if used properly,
#' there will be significant performance boost on large signals with
#' high sampling rates
#' @returns A file-array object of filtered and potentially down-sampled
#' data; see 'Examples' on how to use this function.
#' @examples
#'
#'
#' fs <- 4000
#' time <- seq_len(8000) / fs
#' x <- sin(160 * pi * time) +
#'   sin(1000 * pi * time) * dnorm(time, mean = 1, sd = 0.1) +
#'   0.5 * rnorm(length(time))
#'
#' # envelope
#' result <- gammatone_fast(
#'   x,
#'   sample_rate = fs,
#'   center_frequencies = c(20, 1000),
#'   n_bands = 128,
#'
#'   # default downsample happens after hilbert
#'   downsample = 40
#' )
#'
#'
#' oldpar <- par(mfrow = c(2, 1))
#'
#' plot(
#'   time,
#'   x,
#'   type = "l",
#'   xlab = "Time",
#'   ylab = "",
#'   main = "Original mixed 80Hz and 500Hz"
#' )
#'
#' # only one channel
#' envelope <- subset(result, Channel ~ Channel == 1, drop = TRUE)
#' dnames <-  dimnames(envelope)
#' image(
#'   x = as.numeric(dnames$Time),
#'   y = as.numeric(dnames$Frequency),
#'   z = envelope,
#'   xlab = "Time",
#'   ylab = "Frequency",
#'   main = "Envelope from 20Hz to 1000Hz"
#' )
#'
#'
#'
#' par(oldpar) # reset graphics state
#'
#'
#' @export
gammatone_fast <- function(x, sample_rate, center_frequencies, n_bands,
                           use_hilbert = TRUE, downsample = NA,
                           downsample_before_hilbert = FALSE) {

  # x <- as.vector(ieegio::io_read_mat("~/rave_data/raw_dir/DemoSubject/008/DemoSubjectDatafile008_ch13.mat")$analogTraces)[]
  # sample_rate <- 2000
  # center_frequencies <- c(20, 1000)
  # n_bands <- 128
  # downsample <- 20
  # use_hilbert <- TRUE

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
    phase <- -2 * pi * f * tc

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


  downsample <- as.double(downsample)
  final_sample_rate <- sample_rate

  n_timepoints <- l_x

  if(!is.na(downsample)) {
    if(downsample <= 1) {
      downsample <- NA_real_
    } else {
      n_timepoints <- ceiling(n_timepoints / downsample)
      final_sample_rate <- sample_rate / downsample
    }
  } else {
    downsample <- NA_real_
  }


  root_dir <- file.path(tempdir2(check = TRUE), "ravetools")
  digest <- digest::digest(list(
    center_frequencies = center_frequencies,
    filter_order = filter_order,
    N = N,
    sample_rate = sample_rate
  ))
  x_digest <- digest::digest(list(
    x_digest = digest::digest(x),
    filter_digest = digest,
    use_hilbert = use_hilbert,
    downsample = downsample
  ))

  # if the filters < 256MB
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

  # delay
  delay <- round(tcs * sample_rate)

  # result <- apply(x, 2L, function(x_slice) {
  #   # print(ii <<- ii + 1)
  #   # x_slice <- x[,1]
  #   result <- stats::mvfft(gammatone_filters_fft[drop = FALSE] * fft(postpad(x_slice, N)), inverse = TRUE) / N
  #   result <- Re(result[seq_len(l_x), , drop = FALSE])
  #
  #   # offset delay
  #   result <- sapply(seq_len(n_bands), function(jj) {
  #     res <- result[, jj]
  #     res <- c(res[-seq_len(delay[[jj]])], res[seq_len(delay[[jj]])])
  #
  #     if(isTRUE(downsample > 1)) {
  #       # downsample the filtered results before hilbert or returning to speed up
  #       res <- decimate(res, q = downsample)
  #     }
  #     res
  #   })
  #
  #   if(use_hilbert) {
  #     # apply hilbert transform
  #     result <- abs(hilbert(result))
  #   }
  #
  #   result
  #
  # })
  # dim(result) <- c(l_x, n_bands, c_x)
  # result <- aperm(result, c(1, 3, 2))


  # if( N * n_bands * c_x <= 33554432 ) {
  #   # using around 2-3GB RAM so it's fine
  #   result <- apply(x, 2L, function(x_slice) {
  #     # print(ii <<- ii + 1)
  #     # x_slice <- x[,1]
  #     result <- stats::mvfft(gammatone_filters_fft[drop = FALSE] * fft(postpad(x_slice, N)), inverse = TRUE) / N
  #     result <- Re(result[seq_len(l_x), , drop = FALSE])
  #     if(use_hilbert) {
  #       # apply hilbert transform
  #       result <- abs(hilbert(result))
  #     }
  #     # offset delay
  #     result <- sapply(seq_len(n_bands), function(jj) {
  #       res <- result[, jj]
  #       c(res[-seq_len(delay[[jj]])], res[seq_len(delay[[jj]])])
  #     })
  #   })
  #   dim(result) <- c(l_x, n_bands, c_x)
  #   result <- aperm(result, c(1, 3, 2))
  #
  #   # samp <- result[, 1, ]
  #   # n_plot <- 10000
  #   # matplot(x = seq_len(n_plot) / sample_rate, y = samp[seq_len(n_plot), ], lty = 1, type = "l")
  #   # diagnose_channel(s1 = as.vector(x), s2 = rowSums(samp), srate = sample_rate, which = 1, col = c("black", "red"))
  #   # plot_signals(
  #   #   t(cbind(x, rowSums(samp))), sample_rate = sample_rate, space = 0, duration = 1
  #   # )
  #
  # } else {

  result_path <- file.path(root_dir, sprintf("gammatone-results-%s", x_digest))

  result <- filearray::filearray_load_or_create(
    filebase = result_path,
    dimension = c(n_timepoints, c_x, n_bands),
    type = "double",
    symlink_ok = FALSE,
    partition_size = 1L,
    initialize = FALSE,
    ready = FALSE,
    verbose = FALSE,
    on_missing = function(arr) {

      design_table <- expand.grid(
        x = seq_len(c_x),
        f = seq_len(n_bands)
      )

      fft_x <- filearray::as_filearray(stats::mvfft(postpad(x, N)))
      on.exit({
        fft_x$.mode <- "readwrite"
        fft_x$delete()
      })

      # Avoid using gsignal functions
      ravetools <- asNamespace("ravetools")

      lapply(seq_len(nrow(design_table)), function(row_ii) {
        # row_ii <- 1
        design_row <- design_table[row_ii, ]

        ii <- design_row$x
        jj <- design_row$f

        result <- stats::fft(fft_x[, ii] * gammatone_filters_fft[, jj], inverse = TRUE) / N
        result <- Re(result[seq_len(l_x)])

        # offset delay
        result <- c( result[-seq_len(delay[[jj]])], result[seq_len(delay[[jj]])] )

        # downsample
        # TODO: make sure n_timepoints is accurate and not with possible 1-off?
        if( downsample_before_hilbert && isTRUE(downsample > 1) ){
          result <- ravetools$decimate(x = result, q = downsample)
          if(length(result) > n_timepoints) {
            result <- result[seq_len(n_timepoints)]
          } else if (length(result) < n_timepoints) {
            result <- postpad(result, n_timepoints)
          }
        }

        if(use_hilbert) {
          # apply hilbert transform
          result <- abs(hilbert(result))
        }

        # downsample
        # TODO: make sure n_timepoints is accurate and not with possible 1-off?
        if( !downsample_before_hilbert && isTRUE(downsample > 1) ){
          result <- ravetools$decimate(x = result, q = downsample)
          if(length(result) > n_timepoints) {
            result <- result[seq_len(n_timepoints)]
          } else if (length(result) < n_timepoints) {
            result <- postpad(result, n_timepoints)
          }
        }

        arr[, ii, jj] <- result

        return()
      })

      arr$set_header("center_frequencies", center_frequencies, save = FALSE)
      arr$set_header("filter_order", as.integer(filter_order), save = FALSE)
      arr$set_header("orig_sample_rate", sample_rate, save = FALSE)
      arr$set_header("use_hilbert", use_hilbert, save = FALSE)
      arr$set_header("downsample", downsample, save = FALSE)

      time <- seq_len(n_timepoints) / sample_rate
      if(!is.na(downsample)) {
        time <- time * downsample
      }

      dimnames(arr) <- list(
        Time = time,
        Channel = seq_len(c_x),
        Frequency = center_frequencies
      )
      arr$set_header("ready", TRUE, save = TRUE)
      arr
    }
  )
  # }

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
