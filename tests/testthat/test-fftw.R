require(testthat)
test_that("fftw_r2c", {
  x <- rnorm(1000)
  c <- complex(1000)

  expect_error(ravetools:::fftw_r2c(x, 1, complex(999)))
  expect_error(ravetools:::fftw_r2c(x, 0, complex(500)))

  expect_equal(
    ravetools:::fftw_r2c(x, 1),
    fftwtools::fftw_r2c(x, 1)
  )
  expect_equal(
    ravetools:::fftw_r2c(x, 0),
    fftwtools::fftw_r2c(x, 0)
  )

  c <- complex(1000)
  ravetools:::fftw_r2c(x, 1, c)
  expect_equal(
    c,
    fftwtools::fftw_r2c(x, 1)
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
  a <- ravetools:::mvfftw_r2c(x)
  b <- ravetools:::mvfftw_r2c(x, 1)
  # c <- fftwtools::mvfftw_r2c(x, 1)[1:51,]
  # d <- fftwtools::mvfftw_r2c(x, 1, 1)[1:51,]
  c <- apply(x, 2, stats::fft)[1:51,]
  # expect_equal(a, c)
  expect_equal(b, c)

  e <- complex(length(a))
  ravetools:::mvfftw_r2c(x, ret = e)

  expect_equal(e, as.vector(c))
})
