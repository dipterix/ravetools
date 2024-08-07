test_that("C++ signal filter", {

  initialize_filter <- function(b, a, x) {
    # make sure a, b share the same order
    na <- length(a)
    nb <- length(b)

    if( na > nb ) {
      b <- c(b, rep(0, na - nb))
      n <- na
    } else {
      a <- c(a, rep(0, nb - na))
      n <- nb
    }

    # length of edge transients
    nf <- max(1, 3 * (n - 1))

    # compute the initial condition if n > 1
    if( n > 1 ) {
      z1 <- diag(1, n - 1) - cbind( -a[-1], rbind(diag(1, n - 2), 0) )
      z2 <- b[-1] - b[1] * a[-1]
      z <- solve(z1, z2)
    } else {
      z <- numeric(0)
    }
    list(
      a = a,
      b = b,
      z = z,
      nfilt = n,
      nfact = nf
    )
  }
  myFilter <- function(b, a, x, z) {
    # make sure a, b share the same order
    na <- length(a)
    nb <- length(b)
    if( na > nb ) {
      b <- c(b, rep(0, na - nb))
      n <- na
    } else {
      a <- c(a, rep(0, nb - na))
      n <- nb
    }
    b <- b / a[1]
    a <- a / a[1]
    y <- rep(0, length(x))
    if(missing(z)) {
      z <- rep(0, n - 1)
    }

    for(m in seq_along(y)) {
      xm <- x[m]
      y[m] <- b[1] * xm + z[1]
      ym <- y[m]
      for( i in 2: (n-1)) {
        z[ i-1 ] <- b[i] * xm + z[i] - a[i] * ym
      }
      z[n-1] <- b[n] * xm - a[n] * ym
    }
    list(y, z)
  }

  bf <- gsignal::butter(10, c(0.15, 0.3))
  t <- seq(0, 1, by = 0.005)
  x <- as.double(sin(2*pi*t*2.3)) + rnorm(length(t), mean = t)
  b <- as.double(bf$b)
  a <- as.double(bf$a)
  nx <- length(x)
  init <- initialize_filter(b, a, x)
  nfact <- init$nfact
  z <- as.double(init$z)

  expected <- myFilter(b,a,x,z)
  my_result <- ravetools:::cpp_filter(b,a,x,z)

  testthat::expect_lt(
    max(abs((expected[[1]] - my_result[[1]]) / (expected[[1]] + my_result[[1]] + 1))),
    1e-3
  )
  testthat::expect_lt(
    max(abs((expected[[2]] - my_result[[2]]) / (expected[[2]] + my_result[[2]] + 1))),
    1e-3
  )

})


test_that("multi-signal fftfilt", {
  # sanity-check
  x <- array(sample_signal(5000) + 1:5000, c(1000, 5))
  filter <- design_filter_fir(sample_rate = 2000, high_pass_freq = 0.5, high_pass_trans_freq = 0.4, data_size = 1000)
  b <- filter$b
  expected <- apply(x, 2, function(xi) {
    re <- fftfilt(b, xi, legacy = TRUE)
    expect_vector(re)
    expect_length(re, length(xi))
    re
  })
  expect_equal( fftfilt(b, x, legacy = TRUE), expected )
  expect_equal( fftfilt(b, x, legacy = FALSE), expected )


  # scaled FIR
  filter$a <- 2
  expected <- apply(x, 2, function(xi) {
    filtfilt(filter$b, filter$a, xi)
  })


  actual <- filtfilt(filter$b, filter$a, x)

  testthat::expect_equal(actual, expected)


  # # speed check
  # microbenchmark::microbenchmark(
  #   expected = apply(x, 2, function(xi) {
  #     filtfilt(filter$b, filter$a, xi)
  #   }),
  #   actual = filtfilt(filter$b, filter$a, x),
  #   times = 10
  # )
})
