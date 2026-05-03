require(testthat)

# Helpers --------------------------------------------------------------------

make_signal <- function(n, srate = 200, seed = 1) {
  set.seed(seed)
  t <- seq_len(n) / srate
  sin(2 * pi * 10 * t) +
    0.5 * cos(2 * pi * 5 * t) +
    0.25 * sin(2 * pi * 25 * t) +
    0.1 * rnorm(n)
}

# Read coefficient power from either precision return type.
as_matrix2 <- function(x, ncol_) {
  if (is.null(dim(x))) matrix(x, ncol = ncol_) else x
}
coef_power <- function(re) {
  if (inherits(re, "FileArray")) {
    z <- re[]
    z <- as_matrix2(z, dim(re)[2])
    Mod(z) ^ 2
  } else {
    r <- as_matrix2(re$real[], dim(re$real)[2])
    i <- as_matrix2(re$imag[], dim(re$imag)[2])
    r ^ 2 + i ^ 2
  }
}

coef_complex <- function(re) {
  if (inherits(re, "FileArray")) {
    z <- re[]
    as_matrix2(z, dim(re)[2])
  } else {
    r <- as_matrix2(re$real[], dim(re$real)[2])
    i <- as_matrix2(re$imag[], dim(re$imag)[2])
    r + 1i * i
  }
}

interior_compare <- function(a, b, max_kernel_len) {
  pad <- as.integer(ceiling(max_kernel_len / 2))
  rows <- (pad + 1):(nrow(a) - pad)
  list(a = a[rows, , drop = FALSE], b = b[rows, , drop = FALSE])
}

with_freshtemp <- function(expr) {
  # Force unique kernel cache + output paths so we don't pick up a stale cache.
  old <- options(ravetools.debug = FALSE)
  on.exit(options(old), add = TRUE)
  expr
}

# Plan helper ----------------------------------------------------------------

test_that("wavelet_segment_plan covers the full data range", {
  plan <- ravetools:::wavelet_segment_plan(
    data_length = 1000L, segment_length = 200L, max_kernel_len = 50L
  )
  expect_identical(plan$out_starts[1], 1L)
  expect_identical(tail(plan$out_ends, 1L), 1000L)
  # consecutive segments are contiguous
  diffs <- plan$out_starts[-1] - plan$out_ends[-plan$n_seg]
  expect_true(all(diffs == 1L))
  # local rows have correct length
  lens <- plan$local_ends - plan$local_starts + 1L
  expect_identical(lens, plan$out_ends - plan$out_starts + 1L)
})

test_that("wavelet_segment_plan errors when segment_length <= max_kernel_len", {
  expect_error(
    ravetools:::wavelet_segment_plan(1000L, 50L, 50L),
    "segment_length"
  )
  expect_error(
    ravetools:::wavelet_segment_plan(1000L, 49L, 50L),
    "segment_length"
  )
})

# Numerical equivalence: float ----------------------------------------------

test_that("morlet_wavelet float: segmented matches legacy on the interior", {
  with_freshtemp({
    srate <- 200
    n <- 6000
    x <- make_signal(n, srate)
    freqs <- c(4, 8, 16, 32)

    leg <- ravetools::morlet_wavelet(x, freqs = freqs, srate = srate,
                                     wave_num = c(3, 7))
    seg <- ravetools::morlet_wavelet(x, freqs = freqs, srate = srate,
                                     wave_num = c(3, 7),
                                     segment_length = 1024L)

    # discover max kernel length to choose the comparison interior
    kern <- ravetools:::wavelet_kernels(freqs, srate, c(3, 7))
    max_kl <- max(vapply(kern$kernels, length, integer(1L)))

    pa <- coef_power(leg)
    pb <- coef_power(seg)
    cmp <- interior_compare(pa, pb, max_kl)
    rel <- max(abs(cmp$a - cmp$b) / pmax(abs(cmp$a), 1e-12))
    expect_lt(rel, 1e-4)
  })
})

test_that("morlet_wavelet float: multiple segment_length values agree", {
  with_freshtemp({
    srate <- 200
    n <- 5000
    x <- make_signal(n, srate, seed = 7)
    freqs <- c(5, 12, 30)
    kern <- ravetools:::wavelet_kernels(freqs, srate, c(3, 8))
    max_kl <- max(vapply(kern$kernels, length, integer(1L)))

    leg <- ravetools::morlet_wavelet(x, freqs, srate, c(3, 8))
    pa <- coef_power(leg)

    for (L in c(512L, 1024L, 2048L)) {
      seg <- ravetools::morlet_wavelet(x, freqs, srate, c(3, 8),
                                       segment_length = L)
      pb <- coef_power(seg)
      cmp <- interior_compare(pa, pb, max_kl)
      rel <- max(abs(cmp$a - cmp$b) / pmax(abs(cmp$a), 1e-12))
      expect_lt(rel, 1e-4)
    }
  })
})

# Numerical equivalence: double ---------------------------------------------

test_that("morlet_wavelet double: segmented matches legacy on the interior", {
  with_freshtemp({
    srate <- 200
    n <- 5000
    x <- make_signal(n, srate, seed = 3)
    freqs <- c(4, 10, 20)
    kern <- ravetools:::wavelet_kernels(freqs, srate, c(3, 7))
    max_kl <- max(vapply(kern$kernels, length, integer(1L)))

    leg <- ravetools::morlet_wavelet(x, freqs, srate, c(3, 7),
                                     precision = "double")
    seg <- ravetools::morlet_wavelet(x, freqs, srate, c(3, 7),
                                     precision = "double",
                                     segment_length = 1024L)

    za <- coef_complex(leg)
    zb <- coef_complex(seg)
    cmp <- interior_compare(za, zb, max_kl)
    rel <- max(Mod(cmp$a - cmp$b) / pmax(Mod(cmp$a), 1e-12))
    expect_lt(rel, 1e-9)
  })
})

# Trend modes ---------------------------------------------------------------

test_that("morlet_wavelet: segmented matches legacy across trend modes", {
  with_freshtemp({
    srate <- 200
    n <- 4000
    x <- make_signal(n, srate, seed = 11) + seq_len(n) * 0.001  # add a slope
    freqs <- c(6, 12)
    kern <- ravetools:::wavelet_kernels(freqs, srate, c(3, 7))
    max_kl <- max(vapply(kern$kernels, length, integer(1L)))

    for (tr in c("constant", "linear", "none")) {
      leg <- ravetools::morlet_wavelet(x, freqs, srate, c(3, 7), trend = tr)
      seg <- ravetools::morlet_wavelet(x, freqs, srate, c(3, 7), trend = tr,
                                       segment_length = 1024L)
      pa <- coef_power(leg)
      pb <- coef_power(seg)
      cmp <- interior_compare(pa, pb, max_kl)
      rel <- max(abs(cmp$a - cmp$b) / pmax(abs(cmp$a), 1e-12))
      expect_lt(rel, 1e-4)
    }
  })
})

# Edge cases ----------------------------------------------------------------

test_that("morlet_wavelet: segment_length <= max_kernel_len errors", {
  with_freshtemp({
    srate <- 200
    x <- rnorm(2000)
    # freqs with a low component -> long kernel
    expect_error(
      ravetools::morlet_wavelet(x, freqs = c(2, 10), srate = srate,
                                wave_num = c(3, 7),
                                segment_length = 16L),
      "segment_length"
    )
  })
})

test_that("morlet_wavelet: segment_length >= length(data) silently falls back", {
  with_freshtemp({
    srate <- 200
    n <- 1000
    x <- make_signal(n, srate, seed = 2)
    freqs <- c(8, 16)
    leg <- ravetools::morlet_wavelet(x, freqs, srate, c(3, 7))
    fallback <- ravetools::morlet_wavelet(x, freqs, srate, c(3, 7),
                                          segment_length = n)
    fallback2 <- ravetools::morlet_wavelet(x, freqs, srate, c(3, 7),
                                           segment_length = n + 100L)
    expect_equal(coef_power(leg), coef_power(fallback))
    expect_equal(coef_power(leg), coef_power(fallback2))
  })
})

test_that("morlet_wavelet: segment_length NULL/NA both use legacy path", {
  with_freshtemp({
    srate <- 200
    n <- 1500
    x <- make_signal(n, srate, seed = 5)
    freqs <- c(8, 16)
    leg <- ravetools::morlet_wavelet(x, freqs, srate, c(3, 7))
    null_path <- ravetools::morlet_wavelet(x, freqs, srate, c(3, 7),
                                           segment_length = NULL)
    na_path <- ravetools::morlet_wavelet(x, freqs, srate, c(3, 7),
                                         segment_length = NA)
    expect_equal(coef_power(leg), coef_power(null_path))
    expect_equal(coef_power(leg), coef_power(na_path))
  })
})

test_that("morlet_wavelet: single frequency works", {
  with_freshtemp({
    srate <- 200
    n <- 3000
    x <- make_signal(n, srate, seed = 9)
    freqs <- 10
    kern <- ravetools:::wavelet_kernels(freqs, srate, 5)
    max_kl <- length(kern$kernels[[1]])
    leg <- ravetools::morlet_wavelet(x, freqs, srate, 5)
    seg <- ravetools::morlet_wavelet(x, freqs, srate, 5,
                                     segment_length = 512L)
    pa <- coef_power(leg)
    pb <- coef_power(seg)
    cmp <- interior_compare(pa, pb, max_kl)
    rel <- max(abs(cmp$a - cmp$b) / pmax(abs(cmp$a), 1e-12))
    expect_lt(rel, 1e-4)
  })
})

test_that("morlet_wavelet: segment_length must be a single value", {
  expect_error(
    ravetools::morlet_wavelet(rnorm(1000), 10, 200, 5,
                              segment_length = c(100L, 200L)),
    "single value"
  )
})
