test_that("Baseline", {
  set.seed(0)
  dm <- c(5, 40, 10, 6)
  x <- array(rnorm(prod(dm))^2, dm)

  # trial-baseline
  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 3, 4),
    baseline_indexpoints = c(20:30),
    method = "percentage"
  )
  exp <- aperm(apply(x, c(1, 3, 4), function(slice) {
    (slice / mean(slice[c(20:30)]) - 1) * 100
  }), c(2, 1, 3, 4))
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 3, 4),
    baseline_indexpoints = c(20:30),
    method = "sqrt_percentage"
  )
  exp <- aperm(apply(x, c(1, 3, 4), function(slice) {
    slice <- sqrt(slice)
    (slice / mean(slice[c(20:30)]) - 1) * 100
  }), c(2, 1, 3, 4))
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 3, 4),
    baseline_indexpoints = c(20:30),
    method = "decibel"
  )
  exp <- aperm(apply(x, c(1, 3, 4), function(slice) {
    slice <- 10 * log10(slice)
    slice - mean(slice[c(20:30)])
  }), c(2, 1, 3, 4))
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 3, 4),
    baseline_indexpoints = c(20:30),
    method = "zscore"
  )
  exp <- aperm(apply(x, c(1, 3, 4), function(slice) {
    bl <- slice[c(20:30)]
    (slice - mean(bl)) / sd(bl)
  }), c(2, 1, 3, 4))
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 3, 4),
    baseline_indexpoints = c(20:30),
    method = "sqrt_zscore"
  )
  exp <- aperm(apply(x, c(1, 3, 4), function(slice) {
    slice <- sqrt(slice)
    bl <- slice[c(20:30)]
    (slice - mean(bl)) / sd(bl)
  }), c(2, 1, 3, 4))
  expect_equal(act, exp)

  # trial-baseline, also with zero values

  x <- array(0, dm)
  exp <- x
  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 3, 4),
    baseline_indexpoints = c(20:30),
    method = "percentage"
  )
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 3, 4),
    baseline_indexpoints = c(20:30),
    method = "sqrt_percentage"
  )
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 3, 4),
    baseline_indexpoints = c(20:30),
    method = "decibel"
  )
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x + 1,
    along_dim = 2,
    unit_dims = c(1, 3, 4),
    baseline_indexpoints = c(20:30),
    method = "zscore"
  )
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x + 1,
    along_dim = 2,
    unit_dims = c(1, 3, 4),
    baseline_indexpoints = c(20:30),
    method = "sqrt_zscore"
  )
  expect_equal(act, exp)


  # global baseline
  set.seed(0)
  dm <- c(5, 40, 10, 6)
  x <- array(rnorm(prod(dm))^2, dm)

  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 4),
    baseline_indexpoints = c(20:30),
    method = "percentage"
  )
  exp <- aperm(array(apply(x, c(1, 4), function(slice) {
    (slice / mean(slice[c(20:30), ]) - 1) * 100
  }), dm[c(2, 3, 1, 4)]), c(3, 1, 2, 4))
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 4),
    baseline_indexpoints = c(20:30),
    method = "sqrt_percentage"
  )
  exp <- aperm(array(apply(x, c(1, 4), function(slice) {
    slice <- sqrt(slice)
    (slice / mean(slice[c(20:30), ]) - 1) * 100
  }), dm[c(2, 3, 1, 4)]), c(3, 1, 2, 4))
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 4),
    baseline_indexpoints = c(20:30),
    method = "decibel"
  )
  exp <- aperm(array(apply(x, c(1, 4), function(slice) {
    slice <- 10 * log10(slice)
    (slice - mean(slice[c(20:30), ]))
  }), dm[c(2, 3, 1, 4)]), c(3, 1, 2, 4))
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 4),
    baseline_indexpoints = c(20:30),
    method = "zscore"
  )
  exp <- aperm(array(apply(x, c(1, 4), function(slice) {
    bl <- slice[c(20:30), ]
    (slice - mean(bl)) / sd(bl)
  }), dm[c(2, 3, 1, 4)]), c(3, 1, 2, 4))
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 4),
    baseline_indexpoints = c(20:30),
    method = "sqrt_zscore"
  )
  exp <- aperm(array(apply(x, c(1, 4), function(slice) {
    slice <- sqrt(slice)
    bl <- slice[c(20:30), ]
    (slice - mean(bl)) / sd(bl)
  }), dm[c(2, 3, 1, 4)]), c(3, 1, 2, 4))
  expect_equal(act, exp)

  # global-baseline, also with zero values

  x <- array(0, dm)
  exp <- x
  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 4),
    baseline_indexpoints = c(20:30),
    method = "percentage"
  )
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 4),
    baseline_indexpoints = c(20:30),
    method = "sqrt_percentage"
  )
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 4),
    baseline_indexpoints = c(20:30),
    method = "decibel"
  )
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x + 1,
    along_dim = 2,
    unit_dims = c(1, 4),
    baseline_indexpoints = c(20:30),
    method = "zscore"
  )
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x + 1,
    along_dim = 2,
    unit_dims = c(1, 4),
    baseline_indexpoints = c(20:30),
    method = "sqrt_zscore"
  )
  expect_equal(act, exp)
})

test_that("Baseline - db_zscore and guard conditions", {

  # ── Normal cases ──────────────────────────────────────────────────────────

  set.seed(0)
  dm <- c(5, 40, 10, 6)
  x <- array(rnorm(prod(dm))^2, dm)

  # trial baseline
  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 3, 4),
    baseline_indexpoints = 20:30,
    method = "db_zscore"
  )
  exp <- aperm(apply(x, c(1, 3, 4), function(slice) {
    bl <- log10(slice[20:30])
    (log10(slice) - mean(bl)) / sd(bl)
  }), c(2, 1, 3, 4))
  expect_equal(act, exp)

  # global baseline
  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 4),
    baseline_indexpoints = 20:30,
    method = "db_zscore"
  )
  exp <- aperm(array(apply(x, c(1, 4), function(slice) {
    bl <- log10(slice[20:30, ])
    (log10(slice) - mean(bl)) / sd(bl)
  }), dm[c(2, 3, 1, 4)]), c(3, 1, 2, 4))
  expect_equal(act, exp)

  # ── db_zscore guard: all-zero baseline (bl_sd = 0) ───────────────────────

  x <- array(0, dm)
  exp <- x

  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 3, 4),
    baseline_indexpoints = 20:30,
    method = "db_zscore"
  )
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 4),
    baseline_indexpoints = 20:30,
    method = "db_zscore"
  )
  expect_equal(act, exp)

  # ── db_zscore guard: constant baseline (bl_sd = 0) ───────────────────────
  # All baseline values identical (> 0) -> sd(log10(bl)) = 0 -> output 0

  x <- array(1, dm)
  exp <- array(0, dm)

  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 3, 4),
    baseline_indexpoints = 20:30,
    method = "db_zscore"
  )
  expect_equal(act, exp)

  act <- ravetools::baseline_array(
    x = x,
    along_dim = 2,
    unit_dims = c(1, 4),
    baseline_indexpoints = 20:30,
    method = "db_zscore"
  )
  expect_equal(act, exp)

  # ── decibel guard fix: geometric mean of baseline = 1 ────────────────────
  # mean(log10(bl)) == 0 but there are valid baseline values.
  # Old guard (bl_mean == 0.0) would incorrectly zero out all output.
  # New guard (bl_valid_len == 0) correctly computes 10*(log10(x) - 0).

  dm_small <- c(1, 10, 1, 1)
  # Baseline points 1-4: c(10, 0.1, 10, 0.1) -> log10 values c(1,-1,1,-1) -> mean = 0
  x_small <- array(c(10, 0.1, 10, 0.1, 10, 100, 1000, 0.1, 0.01, 0.001), dm_small)
  exp_small <- array(10 * log10(as.vector(x_small)), dm_small)

  act_small <- ravetools::baseline_array(
    x = x_small,
    along_dim = 2,
    unit_dims = c(1, 3, 4),
    baseline_indexpoints = 1:4,
    method = "decibel"
  )
  expect_equal(act_small, exp_small)
})
