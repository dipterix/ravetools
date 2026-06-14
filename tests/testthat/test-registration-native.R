test_that("trilinear resampling: values and nearest-neighbor default unchanged", {
  set.seed(1)
  nd <- c(12L, 11L, 10L)
  x <- array(rnorm(prod(nd)), nd)
  v2r <- diag(4)
  v2r[1:3, 4] <- c(-3, -2, -1)

  # nearest-neighbor (default) must equal the explicit nearest call
  a <- resample3D(nd, x, t(v2r), t(v2r), na = 0)
  b <- resample3D(nd, x, t(v2r), t(v2r), na = 0, interpolation = 0L)
  expect_identical(a[[1]], b[[1]])
  # identity transform reproduces the input
  expect_equal(as.vector(a[[1]]), as.vector(x), tolerance = 1e-9)

  # trilinear at integer coordinates equals the sampled voxel values
  tl <- resample3D(nd, x, t(v2r), t(v2r), na = 0, interpolation = 1L)
  expect_equal(as.vector(tl[[1]]), as.vector(x), tolerance = 1e-9)

  # trilinear midpoint between two voxels is their average
  half <- v2r
  half[1, 4] <- half[1, 4] + 0.5            # shift sampling by half a voxel in x
  mid <- resample3D(nd, x, t(half), t(v2r), na = NA_real_, interpolation = 1L)[[1]]
  expect_equal(mid[3, 4, 5], (x[3, 4, 5] + x[4, 4, 5]) / 2, tolerance = 1e-9)
})

# helper: smooth multi-blob phantom
phantom <- function(nd, seed = 7, nblob = 6) {
  set.seed(seed)
  g <- expand.grid(x = 0:(nd[1] - 1),
                   y = 0:(nd[2] - 1),
                   z = 0:(nd[3] - 1))
  val <- numeric(nrow(g))
  for (k in seq_len(nblob)) {
    c0 <- runif(3, nd / 3, 2 * nd / 3)
    val <- val + runif(1, .5, 1) *
      exp(-((g$x - c0[1])^2 + (g$y - c0[2])^2 + (g$z - c0[3])^2) / (2 * runif(1, 3, 6)^2))
  }
  array(val, nd)
}

recenter <- function(M, nd, v2r) {
  ctr <- v2r %*% c((nd - 1) / 2, 1)
  Tc <- diag(4)
  Tc[1:3, 4] <- ctr[1:3]
  Tc %*% M %*% solve(Tc)
}

rotz <- function(deg) {
  th <- deg * pi / 180
  R <- diag(4)
  R[1:2, 1:2] <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2)
  R
}

test_that("rigid registration recovers a known translation (all metrics)", {
  nd <- c(40L, 40L, 40L)
  v2r <- diag(4)
  v2r[1:3, 4] <- -20
  target <- phantom(nd)
  shift <- c(3, -2, 1)
  Ttrue <- diag(4)
  Ttrue[1:3, 4] <- shift
  source <- apply_transform3d(
    target,
    v2r,
    solve(Ttrue),
    reference_dim = nd,
    reference_vox2ras = v2r,
    interpolation = "trilinear"
  )

  for (m in c("meansquares", "cc", "mattes")) {
    res <- register_volume3d(
      source, target, v2r, v2r, type = "rigid", metric = m,
      shrink_factors = c(4, 2, 1), smoothing_sigmas = c(2, 1, 0),
      iterations = c(200, 150, 80), sampling_rate = 0.4, verbose = FALSE)
    expect_equal(res$transform[1:3, 4], shift, tolerance = 0.25,
                 info = paste("metric:", m))
    expect_lt(max(abs(res$transform[1:3, 1:3] - diag(3))), 0.03)
  }
})

test_that("rigid registration recovers a known rotation", {
  skip_on_cran()
  nd <- c(48L, 48L, 48L)
  v2r <- diag(4)
  v2r[1:3, 4] <- -24
  target <- phantom(nd)
  Ttrue <- recenter(rotz(12), nd, v2r)
  Ttrue[1:3, 4] <- Ttrue[1:3, 4] + c(2, -3, 1)
  source <- apply_transform3d(target,
                              v2r,
                              solve(Ttrue),
                              reference_dim = nd,
                              reference_vox2ras = v2r)

  res <- register_volume3d(
    source, target, v2r, v2r, type = "rigid", metric = "meansquares",
    shrink_factors = c(4, 2, 1), smoothing_sigmas = c(2, 1, 0),
    iterations = c(300, 200, 100), sampling_rate = 0.4, verbose = FALSE)
  expect_lt(max(abs(res$transform - Ttrue)), 0.1)
})

test_that("affine registration recovers scale + shear + rotation", {
  skip_on_cran()
  nd <- c(48L, 48L, 48L)
  v2r <- diag(4)
  v2r[1:3, 4] <- -24
  target <- phantom(nd)
  A <- diag(4)
  A[1, 1] <- 1.08
  A[2, 2] <- 0.93
  A[1, 2] <- 0.06
  A[2, 3] <- -0.04
  Ttrue <- recenter(A %*% rotz(5), nd, v2r)
  Ttrue[1:3, 4] <- Ttrue[1:3, 4] + c(1, 1, -1)
  source <- apply_transform3d(target,
                              v2r,
                              solve(Ttrue),
                              reference_dim = nd,
                              reference_vox2ras = v2r)

  res <- register_volume3d(
    source, target, v2r, v2r, type = "affine", metric = "cc",
    shrink_factors = c(4, 2, 1), smoothing_sigmas = c(2, 1, 0),
    iterations = c(300, 200, 120), sampling_rate = 0.4, verbose = FALSE)
  expect_lt(max(abs(res$transform - Ttrue)), 0.06)
})

test_that("registration handles different vox2ras (spacing) on source and target", {
  skip_on_cran()
  nd_t <- c(40L, 40L, 40L)
  v2r_t <- diag(4)
  v2r_t[1:3, 4] <- -20
  target <- phantom(nd_t)

  # source on a finer grid (0.5 mm) with an origin shift, holding the same anatomy
  nd_s <- c(80L, 80L, 80L)
  v2r_s <- diag(c(0.5, 0.5, 0.5, 1))
  v2r_s[1:3, 4] <- -20
  source <- apply_transform3d(target,
                              v2r_t,
                              diag(4),
                              reference_dim = nd_s,
                              reference_vox2ras = v2r_s)

  shift <- c(2.5, -1.5, 1)
  Ttrue <- diag(4)
  Ttrue[1:3, 4] <- shift
  source <- apply_transform3d(source, v2r_s, solve(Ttrue),
                              reference_dim = nd_s, reference_vox2ras = v2r_s)

  res <- register_volume3d(
    source, target, v2r_s, v2r_t, type = "rigid", metric = "cc",
    shrink_factors = c(4, 2, 1), smoothing_sigmas = c(2, 1, 0),
    iterations = c(200, 150, 100), sampling_rate = 0.5, verbose = FALSE)
  expect_equal(res$transform[1:3, 4], shift, tolerance = 0.4)

  # warped source lives on the target grid and correlates with the target
  expect_identical(dim(res$image), dim(target))
  ok <- target > 0.05
  expect_gt(cor(res$image[ok], target[ok]), 0.9)
})

test_that("SyN recovers a known smooth deformation", {
  skip_on_cran()
  nd <- c(32L, 32L, 32L)
  v2r <- diag(4)
  v2r[1:3, 4] <- -16
  target <- phantom(nd, seed = 3, nblob = 8)

  # vectorized trilinear sampler (0-based continuous coords)
  tril <- function(vol, cx, cy, cz) {
    d <- dim(vol)
    inb <- cx >= 0 & cy >= 0 & cz >= 0 & cx <= d[1] - 1 & cy <= d[2] - 1 & cz <= d[3] - 1
    cx <- pmin(pmax(cx, 0), d[1] - 1)
    cy <- pmin(pmax(cy, 0), d[2] - 1)
    cz <- pmin(pmax(cz, 0), d[3] - 1)
    x0 <- floor(cx)
    y0 <- floor(cy)
    z0 <- floor(cz)
    x1 <- pmin(x0 + 1, d[1] - 1)
    y1 <- pmin(y0 + 1, d[2] - 1)
    z1 <- pmin(z0 + 1, d[3] - 1)
    fx <- cx - x0
    fy <- cy - y0
    fz <- cz - z0
    idx <- function(x, y, z) {
      x + d[1] * (y + d[2] * z) + 1
    }
    v <- function(x, y, z) {
      vol[idx(x, y, z)]
    }
    c00 <- v(x0, y0, z0) * (1 - fx) + v(x1, y0, z0) * fx
    c10 <- v(x0, y1, z0) * (1 - fx) + v(x1, y1, z0) * fx
    c01 <- v(x0, y0, z1) * (1 - fx) + v(x1, y0, z1) * fx
    c11 <- v(x0, y1, z1) * (1 - fx) + v(x1, y1, z1) * fx
    c0 <- c00 * (1 - fy) + c10 * fy
    c1 <- c01 * (1 - fy) + c11 * fy
    out <- c0 * (1 - fz) + c1 * fz
    out[!inb] <- 0
    out
  }

  g <- expand.grid(x = 0:(nd[1] - 1),
                   y = 0:(nd[2] - 1),
                   z = 0:(nd[3] - 1))
  amp <- 1.4
  L <- 16
  moving <- array(tril(
    target,
    g$x + amp * sin(2 * pi * g$y / L),
    g$y + amp * sin(2 * pi * g$z / L),
    g$z + amp * sin(2 * pi * g$x / L)
  ), nd)

  base_cor <- cor(as.vector(moving), as.vector(target))
  res <- register_volume3d(
    moving, target, v2r, v2r, type = "syn", metric = "cc",
    shrink_factors = c(4, 2, 1), smoothing_sigmas = c(2, 1, 0),
    iterations = c(200, 150, 80), sampling_rate = 1,
    syn_iterations = c(50, 40, 20), syn_sigma = 3, verbose = FALSE)

  expect_true(all(is.finite(res$forward_field)))
  expect_identical(dim(res$forward_field), c(nd, 3L))
  expect_gt(cor(as.vector(res$image), as.vector(target)), base_cor + 0.01)
  expect_lt(mean((res$image - target)^2), mean((moving - target)^2) / 2)
})

test_that("metric improves over the identity transform", {
  nd <- c(40L, 40L, 40L)
  v2r <- diag(4)
  v2r[1:3, 4] <- -20
  target <- phantom(nd)
  Ttrue <- recenter(rotz(8), nd, v2r)
  Ttrue[1:3, 4] <- Ttrue[1:3, 4] + c(2, -2, 1)
  source <- apply_transform3d(target, v2r, solve(Ttrue),
                              reference_dim = nd, reference_vox2ras = v2r)

  res <- register_volume3d(
    source, target, v2r, v2r, type = "rigid", metric = "mattes",
    shrink_factors = c(4, 2, 1), smoothing_sigmas = c(2, 1, 0),
    iterations = c(200, 150, 80), sampling_rate = 0.4, verbose = FALSE)
  # metric is a cost (negative MI); the final value should be below the first
  expect_lt(tail(res$metric_trace, 1), res$metric_trace[1])
})

test_that("masks: all-ones equals unmasked, partial mask runs, bad dims error", {
  nd <- c(32L, 32L, 32L)
  v2r <- diag(4)
  v2r[1:3, 4] <- -16
  target <- phantom(nd, seed = 5)
  shift <- c(2, -1, 1.5)
  Ttrue <- diag(4)
  Ttrue[1:3, 4] <- shift
  source <- apply_transform3d(target, v2r, solve(Ttrue),
                              reference_dim = nd, reference_vox2ras = v2r)

  lin <- function(tm = NULL) {
    register_volume3d(
      source,
      target,
      v2r,
      v2r,
      target_mask = tm,
      type = "rigid",
      metric = "cc",
      shrink_factors = c(4, 2),
      smoothing_sigmas = c(2, 1),
      iterations = c(100, 60),
      sampling_rate = 1,
      verbose = FALSE
    )$transform
  }

  # an all-ones mask includes every voxel, so it must reproduce the unmasked
  # result bit-for-bit (the sampler skips nothing and the RNG stream is the same)
  expect_identical(lin(array(1, nd)), lin())

  # a partial mask still runs and recovers roughly the right translation
  m <- array(0, nd)
  m[6:26, 6:26, 6:26] <- 1
  expect_equal(lin(m)[1:3, 4], shift, tolerance = 0.5)

  # SyN with a mask (target + source) yields a finite field of the right shape
  syn <- register_volume3d(
    source, target, v2r, v2r, target_mask = m, source_mask = m, type = "syn",
    shrink_factors = c(2, 1), smoothing_sigmas = c(1, 0), iterations = c(50, 30),
    syn_iterations = c(10, 5), metric = "cc", verbose = FALSE)
  expect_true(all(is.finite(syn$forward_field)))
  expect_identical(dim(syn$forward_field), c(nd, 3L))

  # a mask must match its image dimensions
  expect_error(
    register_volume3d(source, target, v2r, v2r, target_mask = array(1, c(8, 8, 8)),
                      type = "rigid", verbose = FALSE),
    "must match its image dimensions")
})

test_that("cortical landmark term pulls the SyN warp onto correspondences", {
  nd <- c(32L, 32L, 32L)
  v2r <- diag(4)
  v2r[1:3, 4] <- -16
  set.seed(11)
  img <- array(rnorm(prod(nd)), nd)   # identical fixed = moving

  # corresponding points at integer voxels; target = fixed + a known +2.5 mm in x
  set.seed(12)
  ijk <- unique(cbind(sample(8:24, 90, TRUE), sample(8:24, 90, TRUE), sample(8:24, 90, TRUE)))
  nP <- nrow(ijk)
  ras_fix <- t(v2r %*% rbind(t(ijk), 1))[, 1:3]
  off <- 2.5
  ras_mov <- ras_fix
  ras_mov[, 1] <- ras_mov[, 1] + off

  # points_weight = 0 disables the term -> bit-identical to supplying no points
  pars <- list(type = "syn_only", init_transform = diag(4), metric = "cc",
               shrink_factors = c(2, 1), smoothing_sigmas = c(1, 0),
               syn_iterations = c(8, 8), verbose = FALSE)
  a <- do.call(register_volume3d, c(list(img, img, v2r, v2r), pars))$forward_field
  b <- do.call(register_volume3d, c(list(img, img, v2r, v2r,
        target_points = ras_fix, source_points = ras_mov, points_weight = 0), pars))$forward_field
  expect_identical(a, b)

  # identical images => the image force is ~0, so the landmarks drive the field
  res <- register_volume3d(
    img, img, v2r, v2r, type = "syn_only", init_transform = diag(4),
    target_points = ras_fix, source_points = ras_mov, points_weight = 2,
    shrink_factors = c(2, 1), smoothing_sigmas = c(1, 0),
    syn_iterations = c(40, 40), syn_sigma = 2, metric = "cc", verbose = FALSE)
  ff <- res$forward_field
  disp <- t(vapply(seq_len(nP), function(p) {
    ff[ijk[p, 1] + 1, ijk[p, 2] + 1, ijk[p, 3] + 1, ]
  }, numeric(3)))
  # the recovered displacement at the points should approach the known offset
  expect_lt(mean(sqrt(rowSums((
    matrix(c(off, 0, 0), nP, 3, byrow = TRUE) - disp
  )^2))), 0.8)
  expect_gt(mean(disp[, 1]), 0.7 * off)

  # points must be supplied together and with matching row counts
  expect_error(
    register_volume3d(img, img, v2r, v2r, type = "syn",
                      target_points = ras_fix, verbose = FALSE),
    "supplied together")
  expect_error(
    register_volume3d(img, img, v2r, v2r, type = "syn", target_points = ras_fix,
                      source_points = ras_mov[1:4, ], verbose = FALSE),
    "same number of rows")
})
