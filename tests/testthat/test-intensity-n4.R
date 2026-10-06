# Acceptance tests for bias_correction_n4 (plan checks C1-C7).
#
# Phantom: 64^3, three nested spherical classes (50 / 100 / 150) plus Gaussian
# noise (SD 2.5), multiplied by a known smooth bias with range [0.8, 1.3].
# Tests use shrink_factor = 2 (the plan allows it) so they stay fast.

n4_phantom <- function(nd = c(64L, 64L, 64L), seed = 42L, noise_sd = 2.5) {
  set.seed(seed)
  g <- expand.grid(x = 0:(nd[1] - 1), y = 0:(nd[2] - 1), z = 0:(nd[3] - 1))
  c0 <- (nd - 1) / 2
  r <- sqrt((g$x - c0[1])^2 + (g$y - c0[2])^2 + (g$z - c0[3])^2)
  label <- integer(nrow(g))
  label[r <= 26] <- 1L
  label[r <= 18] <- 2L
  label[r <= 9] <- 3L
  clean <- c(0, 50, 100, 150)[label + 1L]
  noisy <- clean + rnorm(length(clean), sd = noise_sd)

  # smooth bias: gentle gradients plus one wide bump, scaled so that it spans
  # [0.8, 1.3] inside the head (label > 0)
  gx <- (g$x - c0[1]) / nd[1]
  gy <- (g$y - c0[2]) / nd[2]
  gz <- (g$z - c0[3]) / nd[3]
  field <- 0.7 * gx + 0.5 * gy - 0.4 * gz +
    0.6 * exp(-((gx - 0.15)^2 + (gy + 0.1)^2 + gz^2) / (2 * 0.4^2))
  inside <- label > 0L
  bias <- 0.8 + 0.5 * (field - min(field[inside])) / (max(field[inside]) - min(field[inside]))

  list(
    dim = nd,
    clean = array(clean, nd),
    noisy = array(noisy, nd),
    biased = array(noisy * bias, nd),
    bias = array(bias, nd),
    label = array(label, nd),
    mask = array(label > 0L, nd)
  )
}

n4_phantom_vox2ras <- function(nd) {
  v <- diag(4)
  v[1:3, 4] <- -(nd - 1) / 2
  v
}

# within-class coefficient of variation of one class
class_cv <- function(x, label, k) {
  v <- x[label == k]
  stats::sd(v) / mean(v)
}

# pooled within-class coefficient of variation: every class is divided by its
# own mean, so the spread that remains is the within-class one
pooled_cv <- function(x, label) {
  z <- numeric(0)
  for (k in sort(unique(label[label > 0L]))) {
    v <- x[label == k]
    z <- c(z, v / mean(v))
  }
  stats::sd(z) / mean(z)
}

with_threads <- function(n, expr) {
  old <- Sys.getenv("RAVETOOLS_NUM_THREADS", unset = NA)
  on.exit({
    if (is.na(old)) {
      Sys.unsetenv("RAVETOOLS_NUM_THREADS")
    } else {
      Sys.setenv(RAVETOOLS_NUM_THREADS = old)
    }
  }, add = TRUE)
  ravetools_threads(n_threads = n)
  force(expr)
}

test_that("C1/C2: N4 recovers a smooth multiplicative bias and flattens classes", {
  ph <- n4_phantom()
  v2r <- n4_phantom_vox2ras(ph$dim)

  bias_est <- bias_correction_n4(
    ph$biased, vox2ras = v2r, mask = ph$mask, intensity_truncation = NULL,
    shrink_factor = 2, return_bias_field = TRUE)
  expect_equal(dim(bias_est), ph$dim)
  expect_equal(attr(bias_est, "vox2ras"), v2r)
  expect_true(all(is.finite(bias_est)))
  expect_true(all(bias_est > 0))

  # C1: log-bias correlation inside the mask
  r <- stats::cor(log(bias_est[ph$mask]), log(ph$bias[ph$mask]))
  expect_gte(r, 0.95)

  corrected <- bias_correction_n4(
    ph$biased, vox2ras = v2r, mask = ph$mask, intensity_truncation = NULL,
    shrink_factor = 2)
  expect_equal(dim(corrected), ph$dim)
  expect_equal(attr(corrected, "vox2ras"), v2r)
  # corrected == input / bias (same bias field as the return_bias_field call)
  expect_equal(as.vector(corrected), as.vector(ph$biased / bias_est), tolerance = 1e-10)

  # C2: the within-class CV (pooled over the three classes, each normalized
  # by its mean) falls by >= 60% versus the biased input, and every class
  # ends within 1.5x of the unbiased phantom's CV. The 60% criterion is
  # evaluated on the pooled CV because the noise floor alone (2.5 / 150 =
  # 1.7% for the brightest class) exceeds 40% of what a bias confined to
  # [0.8, 1.3] can add inside a 9-voxel core, so a per-class reading of the
  # reduction cannot be met even by a perfect correction.
  cv_in <- pooled_cv(ph$biased, ph$label)
  cv_out <- pooled_cv(corrected, ph$label)
  cv_ref <- pooled_cv(ph$noisy, ph$label)
  expect_lte(cv_out, 0.4 * cv_in)
  expect_lte(cv_out, 1.5 * cv_ref)
  for (k in 1:3) {
    expect_lte(class_cv(corrected, ph$label, k), 1.5 * class_cv(ph$noisy, ph$label, k),
               label = sprintf("class %d: CV vs unbiased", k))
    expect_lt(class_cv(corrected, ph$label, k), class_cv(ph$biased, ph$label, k),
              label = sprintf("class %d: CV reduced", k))
  }
})

test_that("C3: unbiased input yields a nearly flat bias field", {
  ph <- n4_phantom()
  v2r <- n4_phantom_vox2ras(ph$dim)
  bias_est <- bias_correction_n4(
    ph$noisy, vox2ras = v2r, mask = ph$mask, intensity_truncation = NULL,
    shrink_factor = 2, return_bias_field = TRUE)
  expect_lte(stats::sd(log(bias_est[ph$mask])), 0.02)
})

test_that("C4: scaling the input by 10 scales the output and keeps the bias", {
  ph <- n4_phantom()
  v2r <- n4_phantom_vox2ras(ph$dim)
  args <- list(vox2ras = v2r, mask = ph$mask, intensity_truncation = NULL,
               shrink_factor = 2)
  b1 <- do.call(bias_correction_n4,
                c(list(ph$biased, return_bias_field = TRUE), args))
  b2 <- do.call(bias_correction_n4,
                c(list(ph$biased * 10, return_bias_field = TRUE), args))
  expect_lte(max(abs(b2 / b1 - 1)), 1e-4)

  c1 <- do.call(bias_correction_n4, c(list(ph$biased), args))
  c2 <- do.call(bias_correction_n4, c(list(ph$biased * 10), args))
  rel <- abs(c2 - 10 * c1) / (abs(10 * c1) + 1e-8)
  expect_lte(max(rel[ph$mask]), 1e-4)

  # with the default truncation and whole-image mask too, for a non-negative
  # (magnitude-like) image: the quantiles scale with the data and every voxel
  # stays positive after truncation. (With negative voxels left after
  # truncation, these drive the fit with their raw value, as in ITK, so the
  # exact scale invariance holds only for non-negative images.)
  mag <- abs(ph$biased)
  b3 <- bias_correction_n4(mag, vox2ras = v2r, shrink_factor = 2,
                                      return_bias_field = TRUE)
  b4 <- bias_correction_n4(mag * 10, vox2ras = v2r, shrink_factor = 2,
                                      return_bias_field = TRUE)
  expect_lte(max(abs(b4 / b3 - 1)), 1e-4)
})

# R reference of the TruncateIntensity rule of ANTsPy 0.6.3 (probed against
# the installed binary): quantiles from a `bins`-bin histogram of every finite
# voxel spanning [min, max], the minimum being raised to 1e-6 when it is
# exactly 0 (voxels below that bound are left out), linearly interpolated
# within the bin.
ref_truncation_counted <- function(x) {
  v <- x[is.finite(x)]
  lowb <- if (min(v) == 0) 1e-6 else min(v)
  v[v >= lowb]
}
ref_truncation_bounds <- function(x, lo, hi, bins) {
  v <- ref_truncation_counted(x)
  mn <- if (min(x[is.finite(x)]) == 0) 1e-6 else min(v)
  mx <- max(v)
  width <- (mx - mn) / bins
  idx <- pmin(floor((v - mn) / width), bins - 1) + 1
  cnt <- tabulate(idx, nbins = bins)
  cum <- cumsum(cnt) / length(v)
  q <- function(p) {
    j <- which(cum >= p)[1]
    p_prev <- if (j == 1) 0 else cum[j - 1]
    frac <- (p - p_prev) / (cnt[j] / length(v))
    mn + (j - 1) * width + frac * width
  }
  c(q(lo), q(hi), width)
}

test_that("C5: intensity truncation follows the histogram-quantile spec", {
  ph <- n4_phantom()
  v2r <- n4_phantom_vox2ras(ph$dim)
  x <- ph$biased

  ref <- ref_truncation_bounds(x, 0.025, 0.975, 256)
  q <- ravetools:::n4_truncate_quantiles(as.double(x), 0.025, 0.975, 256L)
  expect_length(q, 2L)
  expect_lte(abs(q[1] - ref[1]), 1e-9 * max(abs(x)))
  expect_lte(abs(q[2] - ref[2]), 1e-9 * max(abs(x)))
  # the bounds really are the 2.5 / 97.5 percent points of the counted voxels
  # (here every voxel: the phantom's background noise is partly negative)
  pos <- ref_truncation_counted(x)
  expect_equal(length(pos), length(x))
  expect_lte(abs(mean(pos < q[1]) - 0.025), 0.005)
  expect_lte(abs(mean(pos > q[2]) - 0.025), 0.005)

  # an exact-zero background is left out of the histogram (the minimum 0 is
  # raised to 1e-6), as are values below 1e-6, so the lower bound lies in the
  # tissue and the zeros are lifted to it
  set.seed(7)
  z <- array(0, c(20, 20, 20))
  z[6:15, 6:15, 6:15] <- stats::rnorm(1000, 100, 15)
  z[1, 1, 1:5] <- 5e-7
  qz <- ravetools:::n4_truncate_quantiles(as.double(z), 0.025, 0.975, 256L)
  refz <- ref_truncation_bounds(z, 0.025, 0.975, 256)
  expect_lte(max(abs(qz - refz[1:2])), 1e-9 * max(z))
  expect_gt(qz[1], 50)
  # with a negative voxel present the zeros do count, so the lower bound drops
  zn <- z
  zn[2, 2, 2] <- -1
  qn <- ravetools:::n4_truncate_quantiles(as.double(zn), 0.025, 0.975, 256L)
  expect_lt(qn[1], 1)
  # nothing to count (all zero): no truncation
  expect_true(all(is.na(ravetools:::n4_truncate_quantiles(numeric(27), 0.025, 0.975, 256L))))

  # iterations = 0 => no bias: the output is exactly the truncated input
  out <- bias_correction_n4(
    x, vox2ras = v2r, mask = ph$mask, iterations = 0L, shrink_factor = 2,
    intensity_truncation = c(0.025, 0.975, 256))
  expect_equal(min(out), q[1], tolerance = 1e-12)
  expect_equal(max(out), q[2], tolerance = 1e-12)
  expect_equal(as.vector(out), pmin(pmax(as.vector(x), q[1]), q[2]), tolerance = 1e-12)

  # NULL => no truncation at all
  out0 <- bias_correction_n4(
    x, vox2ras = v2r, mask = ph$mask, iterations = 0L, shrink_factor = 2,
    intensity_truncation = NULL)
  expect_equal(as.vector(out0), as.vector(x), tolerance = 1e-12)
  expect_equal(range(out0), range(x))

  # c(0, 1, bins) clamps only to the range of the counted voxels (here the
  # whole range, since every voxel is counted)
  out1 <- bias_correction_n4(
    x, vox2ras = v2r, mask = ph$mask, iterations = 0L, shrink_factor = 2,
    intensity_truncation = c(0, 1, 64))
  expect_equal(min(out1), min(pos), tolerance = 1e-12)
  expect_equal(max(out1), max(pos), tolerance = 1e-12)
})

test_that("C6: outside-mask intensities and zero-weight regions have no effect", {
  ph <- n4_phantom()
  v2r <- n4_phantom_vox2ras(ph$dim)
  args <- list(vox2ras = v2r, intensity_truncation = NULL, shrink_factor = 2)

  # (a) junk outside the mask
  set.seed(7)
  x1 <- ph$biased
  x2 <- x1
  x2[!ph$mask] <- stats::runif(sum(!ph$mask), -100, 500)
  b1 <- do.call(bias_correction_n4,
                c(list(x1, mask = ph$mask, return_bias_field = TRUE), args))
  b2 <- do.call(bias_correction_n4,
                c(list(x2, mask = ph$mask, return_bias_field = TRUE), args))
  expect_lte(max(abs(log(b1) - log(b2))), 1e-6)
  c1 <- do.call(bias_correction_n4, c(list(x1, mask = ph$mask), args))
  c2 <- do.call(bias_correction_n4, c(list(x2, mask = ph$mask), args))
  expect_lte(max(abs(c1[ph$mask] - c2[ph$mask])), 1e-6)

  # (b) zero weights behave like masking the region out
  slab <- array(FALSE, ph$dim)
  slab[1:24, , ] <- TRUE
  w <- array(1, ph$dim)
  w[slab] <- 0
  bw <- do.call(bias_correction_n4,
                c(list(x1, mask = ph$mask, weight_mask = w, return_bias_field = TRUE), args))
  bm <- do.call(bias_correction_n4,
                c(list(x1, mask = ph$mask & !slab, return_bias_field = TRUE), args))
  expect_lte(max(abs(log(bw) - log(bm))), 1e-6)
  # and the masking really changes something versus the full mask
  expect_gt(max(abs(log(bw) - log(b1))), 1e-4)
})

test_that("C7: results are identical with 1 thread and several threads", {
  # the multi-thread arm uses 2 threads: tests/testthat.R pins the suite to
  # 2 cores (CRAN policy), and the fit is deterministic by construction for
  # any count (fixed chunking, fixed-order reduction)
  ph <- n4_phantom()
  v2r <- n4_phantom_vox2ras(ph$dim)
  run <- function(...) {
    bias_correction_n4(ph$biased, vox2ras = v2r, shrink_factor = 2, ...)
  }
  b1 <- with_threads(1L, run(return_bias_field = TRUE))
  b2 <- with_threads(2L, run(return_bias_field = TRUE))
  expect_identical(as.vector(b1), as.vector(b2))
  c1 <- with_threads(1L, run())
  c2 <- with_threads(2L, run())
  expect_identical(as.vector(c1), as.vector(c2))
  m1 <- with_threads(1L, ravetools:::n4_default_mask(as.double(ph$biased), ph$dim))
  m2 <- with_threads(2L, ravetools:::n4_default_mask(as.double(ph$biased), ph$dim))
  expect_identical(m1, m2)
})

test_that("defaults follow ANTsPy 0.6.3: whole-image mask, one spline span per axis", {
  ph <- n4_phantom(nd = c(20L, 18L, 16L))
  spacing <- c(1, 1.5, 2)
  v2r <- diag(c(spacing, 1))
  base <- list(ph$biased, vox2ras = v2r, shrink_factor = 2, iterations = c(5, 5),
               return_bias_field = TRUE)
  run <- function(...) do.call(bias_correction_n4, c(base, list(...)))
  b_default <- run()
  # mask = NULL lets every voxel drive the fit
  expect_identical(b_default, run(mask = array(TRUE, ph$dim)))
  # spline_distance = NULL is one span over the image itself along each axis
  expect_identical(b_default, run(spline_distance = (ph$dim - 1) * spacing))
  # mask = "auto" is the get_mask recipe applied to the truncated image
  x <- as.double(ph$biased)
  q <- ravetools:::n4_truncate_quantiles(x, 0.025, 0.975, 256L)
  m <- ravetools:::n4_default_mask(pmin(pmax(x, q[1]), q[2]), ph$dim)
  b_auto <- run(mask = "auto")
  expect_identical(b_auto, run(mask = array(m, ph$dim)))
  expect_false(isTRUE(all.equal(b_auto, b_default)))
  # the older ANTsPy 0.3 recipe (automatic mask, 200 mm) still runs
  b_03 <- run(mask = "auto", spline_distance = 200)
  expect_true(all(is.finite(b_03)) && all(b_03 > 0))
  # other values are rejected
  expect_error(run(mask = "head"), "auto")
  expect_error(run(mask = c("auto", "auto")), "auto")
  expect_error(run(spline_distance = -1), "spline_distance")
})

# call bias_correction_n4 with `base` arguments, overriding (or
# adding) the named ones in `...` (a NULL value is passed through as NULL)
n4_call_with <- function(base, ...) {
  extra <- list(...)
  for (nm in names(extra)) base[nm] <- list(extra[[nm]])
  do.call(bias_correction_n4, base)
}

test_that("grid geometry: limits raise informative errors, extremes stay well defined", {
  ph <- n4_phantom(nd = c(16L, 16L, 16L))
  x <- ph$biased
  base <- list(x, mask = ph$mask, intensity_truncation = NULL, shrink_factor = 2,
               iterations = c(2, 2), return_bias_field = TRUE)
  run <- function(...) n4_call_with(base, ...)

  # far larger than the field of view: the virtual padding would exceed
  # 2^31 - 1 voxels (no silent integer overflow, no platform-dependent result)
  expect_error(run(spline_distance = 1e300), "too large")
  expect_error(run(spline_distance = 1e19), "too large")
  expect_error(run(spline_distance = 3e9), "too large")
  expect_error(run(spline_distance = c(200, 200, 1e12)), "dimension 3")
  # far smaller than the voxel size: more than 2^27 spans
  expect_error(run(spline_distance = 1e-12), "too small")
  # the lattice cap counts the doubling across levels and names the remedy
  expect_error(run(iterations = rep(1, 12), spline_distance = 1), "control points")
  expect_error(run(iterations = rep(1, 40), spline_distance = 200), "control points")
  # below the cap, a lattice far denser than the fitting grid (here 123^3
  # control points for 8^3 samples) warns about the units of spline_distance
  expect_warning(bw <- run(iterations = c(1, 1, 1, 1), spline_distance = 1), "millimeters")
  expect_true(all(is.finite(bw)))
  # sane configurations stay silent (default distance; or 4 levels at 30 mm)
  expect_silent(run())
  expect_silent(run(iterations = c(1, 1, 1, 1), spline_distance = 30))

  # a huge but representable distance works: every shrunk sample inside the
  # image drives the fit (full 8^3 box) and the bias is a constant
  geo <- function(dist, spacing = c(1, 1, 1), iters = c(2L, 2L)) {
    ravetools:::n4_bias_field(as.double(x), ph$dim, as.vector(ph$mask), numeric(0),
                              spacing, 2L, iters, 1e-7, rep_len(dist, 3L), 3L, 200L,
                              0.15, 0.01, FALSE)
  }
  res <- geo(1e9)
  expect_equal(res$box_dim, c(8L, 8L, 8L))
  expect_equal(as.numeric(res$padded_dim), rep(1e9 + 1, 3))
  expect_equal(res$control_points, c(5L, 5L, 5L))
  b <- exp(res$log_bias)
  expect_true(all(is.finite(b)))
  expect_lt(diff(range(b)), 1e-6)
  # same geometry whatever the magnitude beyond the field of view
  expect_equal(geo(1e6)$box_dim, c(8L, 8L, 8L))

  # one span over the unpadded image along each axis (the 'ANTsPy' 0.6
  # default mesh): no padding, 4 cubic control points per axis
  sp <- c(1, 1.5, 2)
  res1 <- geo((ph$dim - 1) * sp, spacing = sp, iters = 1L)
  expect_equal(res1$padded_dim, ph$dim)
  expect_equal(res1$control_points, c(4L, 4L, 4L))
  expect_equal(res1$box_dim, c(8L, 8L, 8L))
  # the exact multiple survives a 1e-9 relative round-off of the spacing
  res2 <- geo((ph$dim - 1) * sp, spacing = sp * (1 + 2e-10), iters = 1L)
  expect_equal(res2$padded_dim, ph$dim)
  # tiny voxel spacings: the geometry is still exact (one 1e-6 span over a
  # 1.5e-8 field of view needs 985 voxels of padding), and the default
  # 200-unit distance is rejected instead of overflowing the padded size
  res3 <- geo(1e-6, spacing = rep(1e-9, 3), iters = 1L)
  expect_equal(res3$padded_dim, c(1001L, 1001L, 1001L))
  expect_equal(res3$control_points, c(4L, 4L, 4L))
  expect_equal(res3$box_dim, c(8L, 8L, 8L))
  expect_error(geo(200, spacing = rep(1e-9, 3), iters = 1L), "too large")
})

test_that("logical switches and integer settings are validated", {
  ph <- n4_phantom(nd = c(16L, 16L, 16L))
  x <- ph$biased
  args <- list(x, mask = ph$mask, intensity_truncation = NULL, shrink_factor = 2,
               iterations = c(2, 2))
  run <- function(...) n4_call_with(args, ...)

  # numbers behave like base-R truth values; the bias field is never silently
  # swapped for the corrected image
  b <- run(return_bias_field = TRUE)
  expect_identical(run(return_bias_field = 1), b)
  expect_identical(run(return_bias_field = 0), run())
  expect_identical(run(rescale_intensities = 1), run(rescale_intensities = TRUE))
  for (bad in list(NA, "TRUE", c(TRUE, TRUE), NULL, logical(0))) {
    expect_error(run(return_bias_field = bad), "return_bias_field")
    expect_error(run(verbose = bad), "verbose")
    expect_error(run(rescale_intensities = bad), "rescale_intensities")
  }
  expect_silent(run(verbose = FALSE))

  # no silent truncation of fractional settings, no misleading type errors
  expect_error(run(iterations = c(3.9, 3.9)), "iterations")
  expect_error(run(iterations = c(3, NA)), "iterations")
  expect_error(run(iterations = c(3, Inf)), "iterations")
  expect_error(run(iterations = "3"), "iterations")
  expect_error(bias_correction_n4(x, intensity_truncation = c(0.025, 0.975, 256.7)),
               "truncation")
  expect_error(bias_correction_n4(x, intensity_truncation = c(0.025, NA, 256)),
               "truncation")
  expect_error(bias_correction_n4(array("a", dim(x))), "numeric")
  expect_error(run(spline_distance = "200"), "spline_distance")
  expect_error(run(histogram_bins = 2.5), "histogram_bins")
  expect_error(run(tolerance = -1), "tolerance")

  # a logical volume is accepted (converted to 0 / 1)
  lg <- array(x > 60, dim(x))
  out <- bias_correction_n4(lg, mask = ph$mask, shrink_factor = 2,
                                       iterations = c(2, 2), intensity_truncation = NULL)
  expect_equal(dim(out), dim(x))
  expect_true(all(is.finite(out)))
})

test_that("default mask reproduces the get_mask recipe on the phantom", {
  ph <- n4_phantom()
  m <- ravetools:::n4_default_mask(as.double(ph$biased), ph$dim)
  expect_type(m, "logical")
  expect_length(m, prod(ph$dim))
  agree <- mean(m == as.vector(ph$mask))
  expect_gte(agree, 0.99)
  # the mask is a single filled blob: no holes, one component
  expect_gt(sum(m), 0.9 * sum(ph$mask))
})

# brute-force R reference for binary erosion / dilation with the ITK ball and
# boundary rules (erosion: outside counts as foreground; dilation: background)
ref_morph <- function(m, radius, erode) {
  d <- dim(m)
  offs <- expand.grid(dx = -radius:radius, dy = -radius:radius, dz = -radius:radius)
  offs <- offs[offs$dx^2 + offs$dy^2 + offs$dz^2 <= (radius + 0.5)^2, ]
  out <- array(FALSE, d)
  for (k in seq_len(d[3])) for (j in seq_len(d[2])) for (i in seq_len(d[1])) {
    if (erode && !m[i, j, k]) next
    if (!erode && m[i, j, k]) {
      out[i, j, k] <- TRUE
      next
    }
    ii <- i + offs$dx
    jj <- j + offs$dy
    kk <- k + offs$dz
    ok <- ii >= 1 & ii <= d[1] & jj >= 1 & jj <= d[2] & kk >= 1 & kk <= d[3]
    vals <- m[cbind(ii[ok], jj[ok], kk[ok])]
    out[i, j, k] <- if (erode) all(vals) else any(vals)
  }
  out
}

test_that("mask cleanup primitives match brute-force references", {
  set.seed(3)
  d <- c(11L, 10L, 9L)
  m <- array(stats::runif(prod(d)) > 0.35, d)
  for (r in 1:2) {
    er <- ravetools:::n4_mask_morph(as.vector(m), d, r, TRUE)
    expect_identical(er, as.vector(ref_morph(m, r, TRUE)))
    di <- ravetools:::n4_mask_morph(as.vector(m), d, r, FALSE)
    expect_identical(di, as.vector(ref_morph(m, r, FALSE)))
  }

  # largest component: two blobs (sizes 60 and 1000) + a tiny speck;
  # hole filling: a cavity inside the big blob and a notch open to the outside
  d <- c(20L, 20L, 20L)
  big <- array(FALSE, d)
  big[3:12, 3:12, 3:12] <- TRUE               # 1000 voxels
  small <- array(FALSE, d)
  small[15:19, 15:17, 15:18] <- TRUE          # 60 voxels
  speck <- array(FALSE, d)
  speck[18, 2, 2] <- TRUE
  m <- big | small | speck
  lc <- ravetools:::n4_mask_cleanup(as.vector(m), d, 0L, TRUE, 50L)
  expect_identical(lc, as.vector(big))
  lc2 <- ravetools:::n4_mask_cleanup(as.vector(m), d, 0L, TRUE, 2000L)
  expect_false(any(lc2))

  cav <- big
  cav[6:8, 6:8, 6:8] <- FALSE                 # enclosed cavity -> filled
  cav[3:12, 7, 3] <- FALSE                    # groove reaching the face -> stays
  filled <- ravetools:::n4_mask_cleanup(as.vector(cav), d, 0L, FALSE, 50L)
  expected <- cav
  expected[6:8, 6:8, 6:8] <- TRUE
  expect_identical(filled, as.vector(expected))
})

test_that("B-spline lattice refinement reproduces the same function", {
  set.seed(11)
  for (order in 1:3) {
    ld <- c(5L, 6L, 7L) + order - 3L
    ld <- pmax(ld, order + 1L)
    lat <- stats::rnorm(prod(ld))
    grid <- c(17L, 19L, 23L)
    f0 <- ravetools:::n4_bspline_evaluate(lat, ld, order, grid)
    ref <- ravetools:::n4_bspline_refine(lat, ld, order)
    expect_length(ref, prod(2L * ld - order))
    f1 <- ravetools:::n4_bspline_evaluate(ref, 2L * ld - order, order, grid)
    expect_lte(max(abs(f0 - f1)), 1e-12)
  }
  # evaluation of a constant lattice is the constant (partition of unity)
  f <- ravetools:::n4_bspline_evaluate(rep(2.5, 4 * 4 * 4), c(4L, 4L, 4L), 3L, c(9L, 9L, 9L))
  expect_lte(max(abs(f - 2.5)), 1e-12)
})

test_that("rescaling, verbose output, other spline orders and per-axis distances", {
  ph <- n4_phantom(nd = c(24L, 24L, 24L))
  v2r <- n4_phantom_vox2ras(ph$dim)
  base <- list(ph$biased, vox2ras = v2r, mask = ph$mask, intensity_truncation = NULL,
               shrink_factor = 2, iterations = c(10, 10), spline_distance = 30)

  plain <- do.call(bias_correction_n4, base)
  resc <- do.call(bias_correction_n4, c(base, rescale_intensities = TRUE))
  # inside the mask the rescaled image spans exactly the input range and is an
  # affine function of the plain correction; outside it is untouched
  expect_equal(range(resc[ph$mask]), range(ph$biased[ph$mask]), tolerance = 1e-10)
  fit <- stats::lm(resc[ph$mask] ~ plain[ph$mask])
  expect_lt(max(abs(stats::residuals(fit))), 1e-8)
  expect_identical(resc[!ph$mask], plain[!ph$mask])

  # verbose prints the geometry and one line per iteration
  out <- utils::capture.output(
    do.call(bias_correction_n4, c(base, verbose = TRUE)),
    type = "output")
  expect_true(any(grepl("control points", out)))
  expect_gte(sum(grepl("convergence", out)), 20L)

  # linear and quadratic splines, and anisotropic spline distances, run and
  # still recover the bias direction
  for (ord in 1:2) {
    b <- do.call(bias_correction_n4,
                 c(base, spline_order = ord, return_bias_field = TRUE))
    expect_true(all(is.finite(b)) && all(b > 0))
    expect_gt(stats::cor(log(b[ph$mask]), log(ph$bias[ph$mask])), 0.8)
  }
  b3 <- do.call(bias_correction_n4,
                c(base[-which(names(base) == "spline_distance")],
                  spline_distance = list(c(30, 40, 50)), return_bias_field = TRUE))
  expect_gt(stats::cor(log(b3[ph$mask]), log(ph$bias[ph$mask])), 0.8)

  # NA voxels (here inside the mask) do not break the fit, stay NA in the
  # output, and every other voxel is still corrected
  xna <- ph$biased
  xna[1:3, 1:3, 1:3] <- NA
  cna <- do.call(bias_correction_n4, c(list(xna), base[-1]))
  expect_true(all(is.na(cna[1:3, 1:3, 1:3])))
  expect_true(all(is.finite(cna[!is.na(xna)])))
  # (dropping 27 driving voxels changes the fit by a small amount only)
  expect_lt(max(abs(cna[!is.na(xna)] - plain[!is.na(xna)])), 1e-2 * max(abs(plain)))
})

test_that("argument validation and attribute handling", {
  ph <- n4_phantom(nd = c(16L, 16L, 16L))
  x <- ph$biased
  expect_error(bias_correction_n4(1:10), "3D")
  expect_error(bias_correction_n4(x, mask = array(TRUE, c(8, 8, 8))), "mask")
  expect_error(bias_correction_n4(x, weight_mask = array(-1, dim(x))), "weight")
  expect_error(bias_correction_n4(x, intensity_truncation = c(0.5, 0.2, 10)), "truncation")
  expect_error(bias_correction_n4(x, shrink_factor = 0), "shrink")
  expect_error(bias_correction_n4(x, spline_order = 5), "spline_order")
  expect_error(bias_correction_n4(x, iterations = integer(0)), "iterations")
  expect_error(bias_correction_n4(x, mask = array(FALSE, dim(x))), "mask")
  expect_error(bias_correction_n4(x, vox2ras = diag(3)), "vox2ras")

  # vox2ras attribute is picked up and propagated; integer input is accepted
  v2r <- n4_phantom_vox2ras(ph$dim)
  xi <- array(as.integer(round(x)), dim(x))
  attr(xi, "vox2ras") <- v2r
  out <- bias_correction_n4(xi, shrink_factor = 2, iterations = c(5, 5))
  expect_equal(attr(out, "vox2ras"), v2r)
  expect_type(out, "double")
  expect_equal(dim(out), dim(x))
  # without any geometry the output carries no vox2ras
  out2 <- bias_correction_n4(x, shrink_factor = 2, iterations = c(5, 5))
  expect_null(attr(out2, "vox2ras"))
})
