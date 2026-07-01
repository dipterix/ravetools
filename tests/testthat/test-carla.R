test_that("carla_zmin_boot matches the R reference", {
  # Pure-R reference for the C++ bootstrap reduction: for each column of `ind`,
  # trial-average the resampled channels, re-reference by the channel mean,
  # correlate (over time) every unreferenced channel against every
  # re-referenced one, Fisher-z, and return the mean z of the most globally
  # anti-correlated channel.
  ref_zmin <- function(sub, ind) {
    d <- dim(sub); ii <- d[1]; n_t <- d[2]; n_tr <- d[3]
    M <- matrix(sub, ii * n_t, n_tr)
    vapply(seq_len(ncol(ind)), function(b) {
      w <- tabulate(ind[, b], n_tr) / nrow(ind)
      Useg_m <- matrix(M %*% w, ii, n_t)
      Uref_m <- Useg_m - rep(colMeans(Useg_m), each = ii)
      r <- stats::cor(t(Useg_m), t(Uref_m))
      diag(r) <- NA_real_
      z <- atanh(r)
      kk <- which.min(rowMeans(z, na.rm = TRUE))
      mean(z[kk, ], na.rm = TRUE)
    }, numeric(1))
  }

  set.seed(101)
  # A few shapes, including the smallest valid subset size (ii = 2).
  for (cfg in list(c(2, 10, 4), c(5, 40, 12), c(16, 120, 25), c(30, 200, 40))) {
    ii <- cfg[1]; n_t <- cfg[2]; n_tr <- cfg[3]; nboot <- 20L
    sub <- array(stats::rnorm(ii * n_t * n_tr), dim = c(ii, n_t, n_tr))
    ind <- matrix(sample.int(n_tr, n_tr * nboot, replace = TRUE),
                  nrow = n_tr, ncol = nboot)
    got <- carla_zmin_boot(sub, ind)
    expect_equal(length(got), nboot)
    expect_equal(got, ref_zmin(sub, ind), tolerance = 1e-9)
  }
})

test_that("carla_zmin_boot is deterministic and thread-count invariant", {
  set.seed(202)
  sub <- array(stats::rnorm(20 * 150 * 30), dim = c(20, 150, 30))
  ind <- matrix(sample.int(30, 30 * 40, replace = TRUE), nrow = 30, ncol = 40)

  ravetools_threads(n_threads = 1L)
  one <- carla_zmin_boot(sub, ind)
  ravetools_threads(n_threads = 4L)
  four <- carla_zmin_boot(sub, ind)
  expect_equal(one, four, tolerance = 1e-12)
})

test_that("carla() runs end-to-end and selects the non-responsive channels", {
  set.seed(303)
  nchan <- 16L; n_t <- 200L; ntrial <- 30L
  x <- array(stats::rnorm(nchan * n_t * ntrial), dim = c(nchan, n_t, ntrial))
  # Channels 1:3 carry a phase-locked evoked response; the rest are noise.
  ep <- 6 * exp(-seq(0, 2, length.out = n_t)) * sin(2 * pi * 5 * seq(0, 0.2, length.out = n_t))
  for (ch in 1:3) for (k in seq_len(ntrial)) x[ch, , k] <- x[ch, , k] + ep

  fit <- carla(x, nboot = 50L)
  expect_type(fit, "list")
  expect_true(all(c("channels", "order", "vars", "n_optimum", "zmin_mean") %in% names(fit)))
  # zmin_mean: n_good x nboot, first row structurally NA, rest finite.
  expect_equal(dim(fit$zmin_mean), c(length(fit$order), 50L))
  expect_true(all(is.na(fit$zmin_mean[1, ])))
  expect_true(all(is.finite(fit$zmin_mean[-1, ])))
  # Responsive channels should NOT be picked into the common-average reference.
  expect_false(any(1:3 %in% fit$channels))
})

test_that("carla() determinism under a fixed seed", {
  set.seed(404)
  x <- array(stats::rnorm(12 * 100 * 20), dim = c(12, 100, 20))
  set.seed(1); a <- carla(x, nboot = 30L, sensitive = TRUE)
  set.seed(1); b <- carla(x, nboot = 30L, sensitive = TRUE)
  expect_equal(a$zmin_mean, b$zmin_mean)
  expect_equal(a$channels, b$channels)
  expect_equal(a$n_optimum, b$n_optimum)
})

test_that("carla() handles single-trial and nboot = 1 inputs", {
  set.seed(505)
  x <- array(stats::rnorm(10 * 80 * 15), dim = c(10, 80, 15))

  # Single trial (matrix input) -> variance ranking, single-pass.
  fit1 <- carla(x[, , 1])
  expect_equal(length(fit1$zmin_mean), length(fit1$order))
  expect_true(fit1$n_optimum >= 1L)

  # Multiple trials but nboot = 1 -> single-pass, no bootstrap.
  fitb1 <- carla(x, nboot = 1L)
  expect_equal(ncol(matrix(fitb1$zmin_mean, ncol = 1L)), 1L)
  expect_true(fitb1$n_optimum >= 1L)
})
