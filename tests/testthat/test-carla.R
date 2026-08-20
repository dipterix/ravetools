test_that("carla_zmin_all matches the R reference", {
  # Pure-R reference for the C++ shared-bootstrap evaluation. `x_ord` is
  # time x trials x channels (channels in ranking order). For each column of
  # `ind`: trial-average every channel once (-> Useg, time x channel); then for
  # every subset size ii = 2..n_good re-reference the first ii channels by their
  # channel mean, correlate (over time) each unreferenced channel against every
  # re-referenced one, Fisher-z, and record the mean z of the most globally
  # anti-correlated channel. Row 1 (subset size 1) is structurally NA.
  ref_zmin_all <- function(x_ord, ind) {
    d <- dim(x_ord)
    n_t <- d[1]
    n_tr <- d[2]
    n_good <- d[3]
    nboot <- ncol(ind)
    out <- matrix(NA_real_, n_good, nboot)
    for (b in seq_len(nboot)) {
      w <- tabulate(ind[, b], n_tr) / nrow(ind)
      Useg <- vapply(seq_len(n_good), function(c) {
        matrix(x_ord[, , c], n_t, n_tr) %*% w
      }, numeric(n_t))                          # n_t x n_good (time x channel)
      for (ii in seq.int(2L, n_good)) {
        U    <- Useg[, seq_len(ii), drop = FALSE]  # n_t x ii
        car  <- rowMeans(U)                        # length n_t (candidate CAR)
        Uref <- U - car                            # re-reference
        r <- stats::cor(U, Uref)                   # ii x ii
        diag(r) <- NA_real_
        z <- atanh(r)
        kk <- which.min(rowMeans(z, na.rm = TRUE))
        out[ii, b] <- mean(z[kk, ], na.rm = TRUE)
      }
    }
    out
  }

  set.seed(101)
  # A few shapes, read as (n_t, n_tr, n_good); smallest subset size ii = 2.
  for (cfg in list(c(30, 8, 4), c(60, 20, 6), c(120, 40, 12), c(200, 50, 20))) {
    n_t <- cfg[1]
    n_tr <- cfg[2]
    n_good <- cfg[3]
    nboot <- 20L
    x_ord <- array(stats::rnorm(n_t * n_tr * n_good), dim = c(n_t, n_tr, n_good))
    ind <- matrix(sample.int(n_tr, n_tr * nboot, replace = TRUE),
                  nrow = n_tr, ncol = nboot)
    got <- carla_zmin_all(x_ord, ind)
    expect_equal(dim(got), c(n_good, nboot))
    expect_true(all(is.na(got[1, ])))
    expect_equal(got, ref_zmin_all(x_ord, ind), tolerance = 1e-8)
  }
})

test_that("carla_zmin_all is deterministic and thread-count invariant", {
  set.seed(202)
  # time x trials x channels
  x_ord <- array(stats::rnorm(150 * 30 * 20), dim = c(150, 30, 20))
  ind <- matrix(sample.int(30, 30 * 40, replace = TRUE), nrow = 30, ncol = 40)

  ravetools_threads(n_threads = 1L)
  one <- carla_zmin_all(x_ord, ind)
  ravetools_threads(n_threads = 4L)
  four <- carla_zmin_all(x_ord, ind)
  expect_equal(one, four, tolerance = 1e-12)
})

test_that("carla() runs end-to-end and selects the non-responsive channels", {
  set.seed(303)
  nchan <- 16L
  n_t <- 200L
  ntrial <- 30L
  # time x trials x channels
  x <- array(stats::rnorm(n_t * ntrial * nchan), dim = c(n_t, ntrial, nchan))
  # Channels 1:3 carry a phase-locked evoked response; the rest are noise.
  ep <- 6 * exp(-seq(0, 2, length.out = n_t)) * sin(2 * pi * 5 * seq(0, 0.2, length.out = n_t))
  for (ch in 1:3) for (k in seq_len(ntrial)) x[, k, ch] <- x[, k, ch] + ep

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
  # time x trials x channels
  x <- array(stats::rnorm(100 * 20 * 12), dim = c(100, 20, 12))
  set.seed(1)
  a <- carla(x, nboot = 30L, sensitive = TRUE)
  set.seed(1)
  b <- carla(x, nboot = 30L, sensitive = TRUE)
  expect_equal(a$zmin_mean, b$zmin_mean)
  expect_equal(a$channels, b$channels)
  expect_equal(a$n_optimum, b$n_optimum)
})

test_that("carla() handles single-trial and nboot = 1 inputs", {
  set.seed(505)
  # time x trials x channels
  x <- array(stats::rnorm(80 * 15 * 10), dim = c(80, 15, 10))

  # Single trial (matrix input, time x channels) -> variance ranking, single-pass.
  fit1 <- carla(x[, 1, ])
  expect_equal(length(fit1$zmin_mean), length(fit1$order))
  expect_true(fit1$n_optimum >= 1L)

  # Multiple trials but nboot = 1 -> single-pass, no bootstrap.
  fitb1 <- carla(x, nboot = 1L)
  expect_equal(ncol(matrix(fitb1$zmin_mean, ncol = 1L)), 1L)
  expect_true(fitb1$n_optimum >= 1L)
})

test_that("carla() accepts a FileArray and matches the in-memory result", {
  skip_if_not_installed("filearray")

  set.seed(606)
  nchan <- 14L
  n_t <- 120L
  ntrial <- 25L
  x <- array(stats::rnorm(n_t * ntrial * nchan), dim = c(n_t, ntrial, nchan))
  ep <- 5 * exp(-seq(0, 2, length.out = n_t)) * sin(2 * pi * 6 * seq(0, 0.2, length.out = n_t))
  for (ch in 1:3) for (k in seq_len(ntrial)) x[, k, ch] <- x[, k, ch] + ep

  xf <- filearray::as_filearray(x)
  on.exit(xf$delete(force = TRUE), add = TRUE)

  for (vref in c(FALSE, TRUE)) {
    set.seed(7)
    mem <- carla(x, nboot = 40L, virtual_reference = vref)
    set.seed(7)
    disk <- carla(xf, nboot = 40L, virtual_reference = vref)
    expect_equal(disk$channels,  mem$channels)
    expect_equal(disk$n_optimum, mem$n_optimum)
    expect_equal(disk$zmin_mean, mem$zmin_mean, tolerance = 1e-9)
    expect_equal(disk$vars,      mem$vars,      tolerance = 1e-9)
  }
})
