require(testthat)
test_that("fftw_r2c", {
  set.seed(1)
  x <- rnorm(1000) + 1
  c <- complex(1000)

  expect_error(ravetools:::fftw_r2c(data = x, HermConj = 1, ret = complex(999)))
  expect_error(ravetools:::fftw_r2c(data = x, HermConj = 0, ret = complex(500)))

  expect_equal(
    ravetools:::fftw_r2c(x, HermConj = 1),
    stats::fft(x)
  )
  expect_equal(
    ravetools:::fftw_r2c(x, HermConj = 0, fftwplanopt = 1),
    stats::fft(x)[1:501]
  )

  # make sure not edit in-place
  set.seed(1)
  xx <- rnorm(1000)
  expect_equal(x, xx + 1)

  c <- complex(1000)
  ravetools:::fftw_r2c(data = x, HermConj = 1, ret = c)
  expect_equal(
    c,
    stats::fft(x)
  )

  expect_equal(
    c[1:501],
    fftwtools::fftw_r2c(x, 0)
  )

})

test_that("mvfftw_r2c", {
  set.seed(1)
  x <- rnorm(1000)
  dim(x) <- c(100,10)
  c <- apply(x, 2, stats::fft)[1:51,]
  a <- ravetools:::mvfftw_r2c(x, 0)
  cc <- fftwtools::mvfftw_r2c(x, 1)[1:51,]
  expect_equal(a, c)
  expect_equal(a, cc)

  set.seed(1)
  xx <- rnorm(1000)
  dim(xx) <- c(100,10)
  expect_equal(x, xx)

  e <- complex(length(a))
  ravetools:::mvfftw_r2c(x, ret = e)
  expect_equal(e, as.vector(c))

  set.seed(1)
  x <- rnorm(1000)
  dim(x) <- c(100,10)
  b <- ravetools:::mvfftw_r2c(x, 1L)
  # d <- fftwtools::mvfftw_r2c(x, 1, 1)[1:51,]
  expect_equal(b, c)
  # expect_equal(b, d)

  set.seed(1)
  xx <- rnorm(1000)
  dim(xx) <- c(100,10)
  expect_equal(x, xx)
})
