test_that("fftw_r2c_2d", {
  # even margins
  x <- array(rnorm(1000), c(20, 50))

  a <- fftwtools::fftw_r2c_2d(x, HermConj = 1)
  b <- ravetools:::fftw_r2c_2d(x, HermConj = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_r2c_2d(x, HermConj = 0)
  b <- ravetools:::fftw_r2c_2d(x, HermConj = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1, 1000)
  a <- fftwtools::fftw_r2c_2d(x, HermConj = 1)
  b <- ravetools:::fftw_r2c_2d(x, HermConj = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_r2c_2d(x, HermConj = 0)
  b <- ravetools:::fftw_r2c_2d(x, HermConj = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1000, 1)
  a <- t(fftwtools::fftw_r2c_2d(t(x), HermConj = 1))
  b <- ravetools:::fftw_r2c_2d(x, HermConj = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- t(fftwtools::fftw_r2c_2d(t(x), HermConj = 0))[1:501, 1, drop = FALSE]
  b <- ravetools:::fftw_r2c_2d(x, HermConj = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)


  # odd
  x <- array(rnorm(1071), c(21, 51))

  a <- fftwtools::fftw_r2c_2d(x, HermConj = 1)
  b <- ravetools:::fftw_r2c_2d(x, HermConj = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_r2c_2d(x, HermConj = 0)
  b <- ravetools:::fftw_r2c_2d(x, HermConj = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1, 1071)
  a <- fftwtools::fftw_r2c_2d(x, HermConj = 1)
  b <- ravetools:::fftw_r2c_2d(x, HermConj = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_r2c_2d(x, HermConj = 0)
  b <- ravetools:::fftw_r2c_2d(x, HermConj = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1071, 1)
  a <- t(fftwtools::fftw_r2c_2d(t(x), HermConj = 1))
  b <- ravetools:::fftw_r2c_2d(x, HermConj = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- t(fftwtools::fftw_r2c_2d(t(x), HermConj = 0))[1:536, 1, drop = FALSE]
  b <- ravetools:::fftw_r2c_2d(x, HermConj = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  # zero length/margin
  x <- numeric(0)
  dim(x) <- c(0, 100)

  b <- ravetools:::fftw_r2c_2d(x, HermConj = 1)
  testthat::expect_equal(dim(b), dim(x))

  b <- ravetools:::fftw_r2c_2d(x, HermConj = 0)
  testthat::expect_equal(dim(b), dim(x))

  dim(x) <- c(100, 0)
  b <- ravetools:::fftw_r2c_2d(x, HermConj = 1)
  testthat::expect_equal(dim(b), dim(x))

  b <- ravetools:::fftw_r2c_2d(x, HermConj = 0)
  testthat::expect_equal(dim(b), c(51, 0))

})



test_that("fftw_c2c_2d", {
  # even margins
  x <- array(rnorm(1000), c(20, 50))

  a <- fftwtools::fftw_c2c_2d(x, inverse = 1)
  b <- ravetools:::fftw_c2c_2d(x, inverse = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_c2c_2d(x, inverse = 0)
  b <- ravetools:::fftw_c2c_2d(x, inverse = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1, 1000)
  a <- fftwtools::fftw_c2c_2d(x, inverse = 1)
  b <- ravetools:::fftw_c2c_2d(x, inverse = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_c2c_2d(x, inverse = 0)
  b <- ravetools:::fftw_c2c_2d(x, inverse = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1000, 1)
  a <- fftwtools::fftw_c2c_2d(x, inverse = 1)
  b <- ravetools:::fftw_c2c_2d(x, inverse = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_c2c_2d(x, inverse = 0)
  b <- ravetools:::fftw_c2c_2d(x, inverse = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)


  # odd
  x <- array(rnorm(1071), c(21, 51))

  a <- fftwtools::fftw_c2c_2d(x, inverse = 1)
  b <- ravetools:::fftw_c2c_2d(x, inverse = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_c2c_2d(x, inverse = 0)
  b <- ravetools:::fftw_c2c_2d(x, inverse = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1, 1071)
  a <- fftwtools::fftw_c2c_2d(x, inverse = 1)
  b <- ravetools:::fftw_c2c_2d(x, inverse = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_c2c_2d(x, inverse = 0)
  b <- ravetools:::fftw_c2c_2d(x, inverse = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1071, 1)
  a <- fftwtools::fftw_c2c_2d(x, inverse = 1)
  b <- ravetools:::fftw_c2c_2d(x, inverse = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_c2c_2d(x, inverse = 0)
  b <- ravetools:::fftw_c2c_2d(x, inverse = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  # zero length/margin
  x <- numeric(0)
  dim(x) <- c(0, 100)

  b <- ravetools:::fftw_c2c_2d(x, inverse = 1)
  testthat::expect_equal(dim(b), dim(x))

  b <- ravetools:::fftw_c2c_2d(x, inverse = 0)
  testthat::expect_equal(dim(b), dim(x))

  dim(x) <- c(100, 0)
  b <- ravetools:::fftw_c2c_2d(x, inverse = 1)
  testthat::expect_equal(dim(b), dim(x))

  b <- ravetools:::fftw_c2c_2d(x, inverse = 0)
  testthat::expect_equal(dim(b), dim(x))

})

test_that("fftw_c2c_2d input immutability and roundtrip", {
  set.seed(7)
  x <- matrix(complex(real = rnorm(200), imaginary = rnorm(200)),
              nrow = 20, ncol = 10)
  x_orig <- x + 0  # deep copy

  fwd <- ravetools:::fftw_c2c_2d(x, inverse = 0)
  testthat::expect_equal(x, x_orig)

  fwd_meas <- ravetools:::fftw_c2c_2d(x, inverse = 0, fftwplanopt = 1L)
  testthat::expect_equal(x, x_orig)
  testthat::expect_equal(fwd_meas, fwd, tolerance = 1e-10)

  inv <- ravetools:::fftw_c2c_2d(fwd, inverse = 1)
  testthat::expect_equal(inv / length(x), x, tolerance = 1e-10)

  testthat::expect_equal(x, x_orig)
})

test_that("fftw_r2c_2d input immutability", {
  set.seed(8)
  x <- matrix(rnorm(200), nrow = 20, ncol = 10)
  x_orig <- x + 0

  ravetools:::fftw_r2c_2d(x, HermConj = 1)
  testthat::expect_equal(x, x_orig)

  ravetools:::fftw_r2c_2d(x, HermConj = 0)
  testthat::expect_equal(x, x_orig)

  ravetools:::fftw_r2c_2d(x, HermConj = 1, fftwplanopt = 1L)
  testthat::expect_equal(x, x_orig)
})

test_that("fftw_*_2d match stats::fft", {
  set.seed(91)
  # real input
  X <- matrix(rnorm(120), nrow = 10, ncol = 12)
  ref <- stats::fft(X + 0i)
  testthat::expect_equal(
    ravetools:::fftw_r2c_2d(X, HermConj = 1), ref, tolerance = 1e-10
  )
  testthat::expect_equal(
    ravetools:::fftw_c2c_2d(X + 0i, inverse = 0), ref, tolerance = 1e-10
  )

  # complex input
  Z <- matrix(complex(real = rnorm(120), imaginary = rnorm(120)),
              nrow = 10, ncol = 12)
  testthat::expect_equal(
    ravetools:::fftw_c2c_2d(Z, inverse = 0),
    stats::fft(Z),
    tolerance = 1e-10
  )
  testthat::expect_equal(
    ravetools:::fftw_c2c_2d(Z, inverse = 1),
    stats::fft(Z, inverse = TRUE),
    tolerance = 1e-10
  )
})
