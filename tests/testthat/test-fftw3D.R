test_that("fftw_r2c_3d", {
  # even margins
  x <- array(rnorm(10000), c(20, 50, 10))
  y <- as.vector(x)

  a <- fftwtools::fftw_r2c_3d(x, HermConj = 1)
  b <- ravetools:::fftw_r2c_3d(x, HermConj = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_r2c_3d(x, HermConj = 0)
  b <- ravetools:::fftw_r2c_3d(x, HermConj = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1, 1000, 10)
  x1 <- x
  a <- fftwtools::fftw_r2c_3d(x, HermConj = 1)
  b <- ravetools:::fftw_r2c_3d(x, HermConj = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_r2c_3d(x, HermConj = 0)
  b <- ravetools:::fftw_r2c_3d(x, HermConj = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1000, 1, 10)
  a <- aperm(fftwtools::fftw_r2c_3d(x1, HermConj = 1), c(2, 1, 3))
  b <- ravetools:::fftw_r2c_3d(x, HermConj = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- a[1:501, 1, 1:10, drop = FALSE]
  b <- ravetools:::fftw_r2c_3d(x, HermConj = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1000, 10, 1)
  a <- aperm(fftwtools::fftw_r2c_3d(x1, HermConj = 1), c(2, 3, 1))
  b <- ravetools:::fftw_r2c_3d(x, HermConj = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- a[1:501, 1:10, 1, drop = FALSE]
  b <- ravetools:::fftw_r2c_3d(x, HermConj = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  # make sure x is not altered
  testthat::expect_equal(as.vector(x), y)

  # odd
  x <- array(rnorm(11781), c(21, 51, 11))
  y <- as.vector(x)
  a <- fftwtools::fftw_r2c_3d(x, HermConj = 1)
  b <- ravetools:::fftw_r2c_3d(x, HermConj = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_r2c_3d(x, HermConj = 0)
  b <- ravetools:::fftw_r2c_3d(x, HermConj = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1, 1071, 11)
  x1 <- x
  a <- fftwtools::fftw_r2c_3d(x, HermConj = 1)
  b <- ravetools:::fftw_r2c_3d(x, HermConj = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_r2c_3d(x, HermConj = 0)
  b <- ravetools:::fftw_r2c_3d(x, HermConj = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1071, 1, 11)
  a <- aperm(fftwtools::fftw_r2c_3d(x1, HermConj = 1), c(2, 1, 3))
  b <- ravetools:::fftw_r2c_3d(x, HermConj = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- a[1:536, 1, 1:11, drop = FALSE]
  b <- ravetools:::fftw_r2c_3d(x, HermConj = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1071, 11, 1)
  a <- aperm(fftwtools::fftw_r2c_3d(x1, HermConj = 1), c(2, 3, 1))
  b <- ravetools:::fftw_r2c_3d(x, HermConj = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- a[1:536, 1:11, 1, drop = FALSE]
  b <- ravetools:::fftw_r2c_3d(x, HermConj = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  testthat::expect_equal(as.vector(x), y)

  # zero length/margin
  x <- numeric(0)
  dim(x) <- c(0, 100, 10)

  b <- ravetools:::fftw_r2c_3d(x, HermConj = 1)
  testthat::expect_equal(dim(b), dim(x))

  b <- ravetools:::fftw_r2c_3d(x, HermConj = 0)
  testthat::expect_equal(dim(b), dim(x))

  dim(x) <- c(100, 0, 10)
  b <- ravetools:::fftw_r2c_3d(x, HermConj = 1)
  testthat::expect_equal(dim(b), dim(x))

  b <- ravetools:::fftw_r2c_3d(x, HermConj = 0)
  testthat::expect_equal(dim(b), c(51, 0, 10))

})


test_that("fftw_c2c_3d", {
  # even margins
  x <- array(rnorm(10000), c(20, 50, 10))

  a <- fftwtools::fftw_c2c_3d(x, inverse = 1)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_c2c_3d(x, inverse = 0)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1, 1000, 10)
  a <- fftwtools::fftw_c2c_3d(x, inverse = 1)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_c2c_3d(x, inverse = 0)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1000, 1, 10)
  a <- fftwtools::fftw_c2c_3d(x, inverse = 1)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_c2c_3d(x, inverse = 0)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1000, 10, 1)
  a <- fftwtools::fftw_c2c_3d(x, inverse = 1)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_c2c_3d(x, inverse = 0)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)


  # odd
  x <- array(rnorm(11781), c(21, 51, 11))

  a <- fftwtools::fftw_c2c_3d(x, inverse = 1)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_c2c_3d(x, inverse = 0)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1, 1071, 11)
  a <- fftwtools::fftw_c2c_3d(x, inverse = 1)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_c2c_3d(x, inverse = 0)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1071, 1, 11)
  a <- fftwtools::fftw_c2c_3d(x, inverse = 1)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_c2c_3d(x, inverse = 0)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  dim(x) <- c(1071, 11, 1)
  a <- fftwtools::fftw_c2c_3d(x, inverse = 1)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 1)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  a <- fftwtools::fftw_c2c_3d(x, inverse = 0)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 0)
  testthat::expect_equal(b, a, tolerance = 1e-10)

  # zero length/margin
  x <- numeric(0)
  dim(x) <- c(0, 100, 10)

  b <- ravetools:::fftw_c2c_3d(x, inverse = 1)
  testthat::expect_equal(dim(b), dim(x))

  b <- ravetools:::fftw_c2c_3d(x, inverse = 0)
  testthat::expect_equal(dim(b), dim(x))

  dim(x) <- c(100, 0, 10)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 1)
  testthat::expect_equal(dim(b), dim(x))

  b <- ravetools:::fftw_c2c_3d(x, inverse = 0)
  testthat::expect_equal(dim(b), dim(x))

  dim(x) <- c(100, 10, 0)
  b <- ravetools:::fftw_c2c_3d(x, inverse = 1)
  testthat::expect_equal(dim(b), dim(x))

  b <- ravetools:::fftw_c2c_3d(x, inverse = 0)
  testthat::expect_equal(dim(b), dim(x))

})

test_that("fftw_c2c_3d input immutability and roundtrip", {
  set.seed(11)
  x <- array(complex(real = rnorm(60), imaginary = rnorm(60)), c(5, 4, 3))
  x_orig <- x + 0

  fwd <- ravetools:::fftw_c2c_3d(x, inverse = 0)
  testthat::expect_equal(x, x_orig)

  fwd_meas <- ravetools:::fftw_c2c_3d(x, inverse = 0, fftwplanopt = 1L)
  testthat::expect_equal(x, x_orig)
  testthat::expect_equal(fwd_meas, fwd, tolerance = 1e-10)

  inv <- ravetools:::fftw_c2c_3d(fwd, inverse = 1)
  testthat::expect_equal(inv / length(x), x, tolerance = 1e-10)

  testthat::expect_equal(x, x_orig)
})

test_that("fftw_r2c_3d input immutability", {
  set.seed(12)
  x <- array(rnorm(60), c(5, 4, 3))
  x_orig <- x + 0

  ravetools:::fftw_r2c_3d(x, HermConj = 1)
  testthat::expect_equal(x, x_orig)

  ravetools:::fftw_r2c_3d(x, HermConj = 0)
  testthat::expect_equal(x, x_orig)

  ravetools:::fftw_r2c_3d(x, HermConj = 1, fftwplanopt = 1L)
  testthat::expect_equal(x, x_orig)
})

test_that("fftw_*_3d match stats::fft", {
  set.seed(95)
  # real input
  X <- array(rnorm(120), dim = c(5, 4, 6))
  ref <- stats::fft(X + 0i)
  testthat::expect_equal(
    ravetools:::fftw_r2c_3d(X, HermConj = 1), ref, tolerance = 1e-10
  )
  testthat::expect_equal(
    ravetools:::fftw_c2c_3d(X + 0i, inverse = 0), ref, tolerance = 1e-10
  )

  # complex input
  Z <- array(complex(real = rnorm(120), imaginary = rnorm(120)),
             dim = c(5, 4, 6))
  testthat::expect_equal(
    ravetools:::fftw_c2c_3d(Z, inverse = 0),
    stats::fft(Z),
    tolerance = 1e-10
  )
  testthat::expect_equal(
    ravetools:::fftw_c2c_3d(Z, inverse = 1),
    stats::fft(Z, inverse = TRUE),
    tolerance = 1e-10
  )
})
