require(testthat)
test_that("fftw_r2c", {
  set.seed(1)
  x <- rnorm(1000) + 1
  c <- complex(1000)

  expect_error(ravetools:::fftw_r2c(
    data = x,
    HermConj = 1,
    ret = complex(999)
  ))
  expect_error(ravetools:::fftw_r2c(
    data = x,
    HermConj = 0,
    ret = complex(500)
  ))

  expect_equal(ravetools:::fftw_r2c(x, HermConj = 1), stats::fft(x))
  expect_equal(ravetools:::fftw_r2c(x, HermConj = 0, fftwplanopt = 1),
               stats::fft(x)[1:501])

  # make sure not edit in-place
  set.seed(1)
  xx <- rnorm(1000)
  expect_equal(x, xx + 1)

  c <- complex(1000)
  ravetools:::fftw_r2c(data = x,
                       HermConj = 1,
                       ret = c)
  expect_equal(c, stats::fft(x))

  expect_equal(c[1:501], fftwtools::fftw_r2c(x, 0))


  expect_equal(x, xx + 1)


})


test_that("mvfftw_r2c", {
  set.seed(1)
  x <- rnorm(1000)
  dim(x) <- c(100, 10)
  c <- apply(x, 2, stats::fft)[1:51, ]
  a <- ravetools:::mvfftw_r2c(x, 0)
  cc <- fftwtools::mvfftw_r2c(x, 1)[1:51, ]
  expect_equal(a, c)
  expect_equal(a, cc)

  set.seed(1)
  xx <- rnorm(1000)
  dim(xx) <- c(100, 10)
  expect_equal(x, xx)

  e <- complex(length(a))
  ravetools:::mvfftw_r2c(x, ret = e)
  expect_equal(e, c)

  set.seed(1)
  x <- rnorm(1000)
  dim(x) <- c(100, 10)
  b <- ravetools:::mvfftw_r2c(x, 1L)
  # d <- fftwtools::mvfftw_r2c(x, 1, 1)[1:51,]
  expect_equal(b, c)
  # expect_equal(b, d)

  set.seed(1)
  xx <- rnorm(1000)
  dim(xx) <- c(100, 10)
  expect_equal(x, xx)
})

test_that("mvfftw_r2c HermConj", {
  set.seed(1)
  x <- matrix(rnorm(1000), nrow = 100, ncol = 10)
  xx <- x + 1

  ref <- stats::mvfft(x)  # 100 x 10 complex

  # --- HermConj=1 returns full n x m matrix matching stats::mvfft ---
  full <- ravetools:::mvfftw_r2c(x, HermConj = 1L)
  expect_true(is.matrix(full))
  expect_equal(dim(full), c(100L, 10L))
  expect_equal(full, ref)

  # --- input immutability ---
  expect_equal(xx - 1, x)

  # --- half-spectrum (default HermConj=0) is consistent prefix ---
  half <- ravetools:::mvfftw_r2c(x)
  expect_equal(half, full[1:51, ])

  # --- pre-allocated ret for HermConj=1 ---
  ret_full <- matrix(0+0i, nrow = 100, ncol = 10)
  ravetools:::mvfftw_r2c(x, HermConj = 1L, ret = ret_full)
  expect_equal(ret_full, ref)

  expect_equal(xx - 1, x)

  # --- fftwplanopt=1 with HermConj=1 (exercises non-ESTIMATE path) ---
  full_meas <- ravetools:::mvfftw_r2c(x, fftwplanopt = 1L, HermConj = 1L)
  expect_equal(full_meas, ref)

  expect_equal(xx - 1, x)

  # --- odd nrows (101): conjugate mirror for odd-length DFT ---
  x_odd <- matrix(rnorm(1010), nrow = 101, ncol = 10)
  ref_odd <- stats::mvfft(x_odd)
  full_odd <- ravetools:::mvfftw_r2c(x_odd, HermConj = 1L)
  expect_equal(dim(full_odd), c(101L, 10L))
  expect_equal(full_odd, ref_odd)
})


test_that("fftw_c2r", {
  set.seed(1)
  x <- rnorm(1000) + 1i * rnorm(1000)
  xx <- x + 1
  c <- double(1000)

  expect_error(ravetools:::fftw_c2r(
    data = x,
    HermConj = 1,
    ret = double(999)
  ))
  expect_error(ravetools:::fftw_c2r(
    data = x,
    HermConj = 0,
    ret = double(1000)
  ))

  expect_equal(x + 1, xx)

  expect_equal(ravetools:::fftw_c2r(x, HermConj = 1),
               fftwtools::fftw_c2r((xx - 1), HermConj = 1, n = 1000))
  expect_equal(x + 1, xx)

  expect_equal(
    ravetools:::fftw_c2r(x, HermConj = 0, fftwplanopt = 1),
    fftwtools::fftw_c2r((xx - 1), HermConj = 0, n = 1998)
  )
  expect_equal(x + 1, xx)

  ravetools:::fftw_c2r(data = x,
                       HermConj = 1,
                       ret = c)
  expect_equal(c, fftwtools::fftw_c2r((xx - 1), HermConj = 1))
  expect_equal(x + 1, xx)

  expect_equal(c, fftwtools::fftw_c2r((xx - 1), HermConj = 1))

  expect_equal(x + 1, xx)


})

test_that("fftw_c2c", {
  set.seed(1)
  x <- rnorm(1000) + 1i * rnorm(1000)
  xx <- x + 1
  c <- complex(1000)

  expect_error(ravetools:::fftw_c2c(data = x, ret = double(1000)))
  expect_error(ravetools:::fftw_c2c(data = x, ret = complex(999)))

  expect_equal(x + 1, xx)

  expect_equal(ravetools:::fftw_c2c(x), fftwtools::fftw_c2c((xx - 1)))
  expect_equal(x + 1, xx)

  expect_equal(ravetools:::fftw_c2c(x, fftwplanopt = 1L, ret = c),
               fftwtools::fftw_c2c((xx - 1)))
  expect_equal(x + 1, xx)

  expect_equal(
    ravetools:::fftw_c2c(
      x,
      inverse = TRUE,
      fftwplanopt = 1L,
      ret = c
    ),
    fftwtools::fftw_c2c((xx - 1), inverse = TRUE)
  )
  expect_equal(x + 1, xx)

  expect_equal(c, fftwtools::fftw_c2c((xx - 1), inverse = TRUE))

  # inplace with ret == data
  expect_val <- fftwtools::fftw_c2c((xx - 1), inverse = TRUE)
  expect_equal(ravetools:::fftw_c2c(x, inverse = TRUE, ret = x), expect_val)

  expect_equal(x, expect_val)


})

test_that("mvfftw_c2c", {
  set.seed(1)
  x <- matrix(rnorm(1000) + 1i * rnorm(1000), nrow = 100, ncol = 10)
  xx <- x + 1

  # --- input validation ---
  expect_error(ravetools:::mvfftw_c2c(data = x, ret = double(1000)))
  expect_error(ravetools:::mvfftw_c2c(data = x, ret = complex(999)))

  # input not modified
  expect_equal(x + 1, xx)

  # --- forward transform matches stats::mvfft ---
  fwd <- ravetools:::mvfftw_c2c(x)
  expect_true(is.matrix(fwd))
  expect_equal(dim(fwd), dim(x))
  expect_equal(fwd, stats::mvfft(xx - 1))

  # input still unmodified
  expect_equal(x + 1, xx)

  # --- inverse transform matches stats::mvfft(inverse = TRUE) ---
  inv <- ravetools:::mvfftw_c2c(x, inverse = 1L)
  expect_true(is.matrix(inv))
  expect_equal(dim(inv), dim(x))
  expect_equal(inv, stats::mvfft(xx - 1, inverse = TRUE))

  expect_equal(x + 1, xx)

  # --- round-trip: inverse(forward(x)) / nrow == x ---
  rt <- ravetools:::mvfftw_c2c(
    ravetools:::mvfftw_c2c(xx - 1),
    inverse = 1L
  ) / nrow(x)
  expect_equal(rt, xx - 1)

  # --- pre-allocated ret (forward) ---
  ret_fwd <- matrix(0 + 0i, nrow = 100, ncol = 10)
  result <- ravetools:::mvfftw_c2c(x, ret = ret_fwd)
  expect_equal(result, stats::mvfft(xx - 1))
  expect_equal(ret_fwd, stats::mvfft(xx - 1))   # ret filled in-place

  expect_equal(x + 1, xx)

  # --- pre-allocated ret (inverse) ---
  ret_inv <- matrix(0 + 0i, nrow = 100, ncol = 10)
  result_inv <- ravetools:::mvfftw_c2c(x, inverse = 1L, ret = ret_inv)
  expect_equal(result_inv, stats::mvfft(xx - 1, inverse = TRUE))
  expect_equal(ret_inv, stats::mvfft(xx - 1, inverse = TRUE))

  expect_equal(x + 1, xx)

  # --- non-ESTIMATE plan (fftwplanopt = 1) ---
  fwd_meas <- ravetools:::mvfftw_c2c(x, fftwplanopt = 1L)
  expect_equal(fwd_meas, stats::mvfft(xx - 1))

  inv_meas <- ravetools:::mvfftw_c2c(x, inverse = 1L, fftwplanopt = 1L)
  expect_equal(inv_meas, stats::mvfft(xx - 1, inverse = TRUE))

  expect_equal(x + 1, xx)
})

test_that("mvfftw_c2r", {
  set.seed(1)
  x <- matrix(rnorm(1000), nrow = 100, ncol = 10)  # real input
  xx <- x + 1

  # --- half-spectrum from r2c (nrows=51 → retrows=100) ---
  half <- ravetools:::mvfftw_r2c(x)                 # 51 × 10 complex
  half2 <- ravetools:::mvfftw_r2c(xx - 1)

  # input validation
  expect_error(ravetools:::mvfftw_c2r(data = half, ret = complex(510)))  # wrong type
  expect_error(ravetools:::mvfftw_c2r(data = half, ret = double(509)))   # wrong length
  expect_error(ravetools:::mvfftw_c2r(data = half, retrows = 99L))       # invalid retrows (nc=51 → only 100 or 101 allowed)

  expect_equal(xx - 1, x)
  expect_equal(half, half2)
  gc()
  expect_equal(half, half2)

  # --- default retrows=0 auto-resolves to even (100) ---
  result <- ravetools:::mvfftw_c2r(half)
  expect_equal(half, half2)

  expect_true(is.matrix(result))
  expect_equal(dim(result), c(100L, 10L))
  # unnormalised inverse: matches stats::fft round-trip (Re only, no /n)
  expect_equal(result, Re(stats::mvfft(stats::mvfft(x), inverse = TRUE)))

  expect_equal(xx - 1, x)

  # --- round-trip: c2r(r2c(x)) / nrow == x ---
  rt <- ravetools:::mvfftw_c2r(ravetools:::mvfftw_r2c(x)) / nrow(x)
  expect_equal(rt, x)

  # --- explicit even retrows ---
  expect_equal(half, half2)
  result_even <- ravetools:::mvfftw_c2r(half, retrows = 100L)
  expect_equal(half, half2)
  expect_equal(result_even, result)

  # --- odd retrows (101 rows of output) ---
  # build half-spectrum from a 101-row real matrix
  x_odd <- matrix(rnorm(1010), nrow = 101, ncol = 10)
  half_odd <- ravetools:::mvfftw_r2c(x_odd)         # 51 × 10 (nc = 51)
  result_odd <- ravetools:::mvfftw_c2r(half_odd, retrows = 101L)
  expect_true(is.matrix(result_odd))
  expect_equal(dim(result_odd), c(101L, 10L))
  rt_odd <- ravetools:::mvfftw_c2r(half_odd, retrows = 101L) / 101
  expect_equal(rt_odd, x_odd)

  # --- pre-allocated ret ---
  ret_mat <- matrix(0.0, nrow = 100, ncol = 10)
  result2 <- ravetools:::mvfftw_c2r(half, ret = ret_mat)
  expect_equal(result2, result)
  expect_equal(ret_mat, result)      # filled in-place

  expect_equal(xx - 1, x)

  # --- non-ESTIMATE plan (fftwplanopt = 1) exercises the fixed memcpy ---
  result_meas <- ravetools:::mvfftw_c2r(half, fftwplanopt = 1L)
  expect_equal(result_meas, result)

  rt_meas <- ravetools:::mvfftw_c2r(ravetools:::mvfftw_r2c(x, fftwplanopt = 1L),
                                    fftwplanopt = 1L) / nrow(x)
  expect_equal(rt_meas, x)

  expect_equal(xx - 1, x)
})
