# apply_transform3d_volume / apply_transform3d_points, the displacement-field
# inverse of register_volume3d, and the exposed SyN controls.
# All data are small synthetic volumes (<= 56^3).

# ---- helpers ----------------------------------------------------------------

ta_phantom <- function(nd, seed = 7, nblob = 6) {
  set.seed(seed)
  g <- expand.grid(x = 0:(nd[1] - 1), y = 0:(nd[2] - 1), z = 0:(nd[3] - 1))
  val <- numeric(nrow(g))
  for (k in seq_len(nblob)) {
    c0 <- runif(3, nd / 3, 2 * nd / 3)
    val <- val + runif(1, .5, 1) *
      exp(-((g$x - c0[1])^2 + (g$y - c0[2])^2 + (g$z - c0[3])^2) / (2 * runif(1, 3, 6)^2))
  }
  array(val, nd)
}

# vectorized trilinear sampler at 0-based continuous voxel coordinates; zero
# outside the [0, n-1] box (the same rule as the C++ samplers)
ta_tril <- function(vol, cx, cy, cz) {
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
  v <- function(x, y, z) vol[x + d[1] * (y + d[2] * z) + 1]
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

# RAS coordinates (N x 3) of every voxel centre of a grid
ta_grid_ras <- function(nd, v2r) {
  g <- expand.grid(x = 0:(nd[1] - 1), y = 0:(nd[2] - 1), z = 0:(nd[3] - 1))
  t(v2r %*% rbind(g$x, g$y, g$z, 1))[, 1:3]
}

# sample a (nx, ny, nz, 3) RAS field (with vox2ras attribute) at RAS points
ta_sample_field <- function(field, ras) {
  v2r <- attr(field, "vox2ras")
  vox <- t(solve(v2r) %*% rbind(t(ras), 1))[, 1:3, drop = FALSE]
  d <- dim(field)[1:3]
  comp <- function(k) ta_tril(array(field[, , , k], d), vox[, 1], vox[, 2], vox[, 3])
  cbind(comp(1), comp(2), comp(3))
}

ta_apply_affine <- function(M, p) t(M %*% rbind(t(p), 1))[, 1:3, drop = FALSE]

ta_rotz <- function(deg) {
  th <- deg * pi / 180
  R <- diag(4)
  R[1:2, 1:2] <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2)
  R
}

ta_recenter <- function(M, nd, v2r) {
  ctr <- v2r %*% c((nd - 1) / 2, 1)
  Tc <- diag(4)
  Tc[1:3, 4] <- ctr[1:3]
  Tc %*% M %*% solve(Tc)
}

ta_inner <- function(nd, border = 3L) {
  m <- array(FALSE, nd)
  m[(border + 1):(nd[1] - border), (border + 1):(nd[2] - border), (border + 1):(nd[3] - border)] <- TRUE
  m
}

ta_norm <- function(m) sqrt(rowSums(m^2))

# smooth radial window (1 inside r0, cosine taper to 0 at r1, in units of the
# half-extent) so a phantom vanishes at the volume border
ta_window <- function(nd, r0 = 0.75, r1 = 1) {
  g <- expand.grid(x = 0:(nd[1] - 1), y = 0:(nd[2] - 1), z = 0:(nd[3] - 1))
  c0 <- (nd - 1) / 2
  r <- sqrt(((g$x - c0[1]) / c0[1])^2 + ((g$y - c0[2]) / c0[2])^2 + ((g$z - c0[3]) / c0[3])^2)
  w <- ifelse(r <= r0, 1, ifelse(r >= r1, 0, cos(pi / 2 * (r - r0) / (r1 - r0))^2))
  array(w, nd)
}

# one toy SyN registration shared by several tests (computed once per file):
# a windowed blob phantom deformed by a smooth one-period sinusoid of 2-voxel
# amplitude (smooth relative to the grid, so the inverse field interpolates
# accurately between nodes, and zero at the border, so the registration has
# no edge structure to chase)
ta_cache <- new.env(parent = emptyenv())
ta_syn_toy <- function() {
  if (!is.null(ta_cache$toy)) return(ta_cache$toy)
  nd <- c(48L, 48L, 48L)
  v2r <- diag(4)
  v2r[1:3, 4] <- -24
  target <- ta_phantom(nd, seed = 3, nblob = 8) * ta_window(nd)
  g <- expand.grid(x = 0:(nd[1] - 1), y = 0:(nd[2] - 1), z = 0:(nd[3] - 1))
  amp <- 2
  L <- 48
  moving <- array(ta_tril(target,
                          g$x + amp * sin(2 * pi * g$y / L),
                          g$y + amp * sin(2 * pi * g$z / L),
                          g$z + amp * sin(2 * pi * g$x / L)), nd)
  args <- list(moving, target, v2r, v2r, type = "syn", metric = "cc",
               shrink_factors = c(4, 2, 1), smoothing_sigmas = c(2, 1, 0),
               iterations = c(100, 60, 30), sampling_rate = 1,
               syn_iterations = c(40, 30, 15), syn_sigma = 3, verbose = FALSE)
  res <- do.call(register_volume3d, args)
  lab_src <- (moving > 0.3) + (moving > 0.6)
  lab_tgt <- (target > 0.3) + (target > 0.6)
  res_multi <- register_volume3d(
    list(moving, lab_src), list(target, lab_tgt), v2r, v2r, type = "syn",
    metric = c("cc", "meansquares"), weights = c(1, 0.3),
    interpolation = c("trilinear", "nearest"),
    shrink_factors = c(4, 2, 1), smoothing_sigmas = c(2, 1, 0),
    iterations = c(60, 40, 20), sampling_rate = 1,
    syn_iterations = c(20, 15, 5), syn_sigma = 3, verbose = FALSE)
  ta_cache$toy <- list(nd = nd, v2r = v2r, target = target, moving = moving,
                       lab_src = lab_src, args = args, res = res, res_multi = res_multi,
                       base_cor = cor(as.vector(moving), as.vector(target)))
  ta_cache$toy
}

# ---- A1 ---------------------------------------------------------------------

test_that("A1: a single affine matches apply_transform3d on a rotated anisotropic grid", {
  set.seed(101)
  nd <- c(24L, 20L, 22L)
  vol <- ta_phantom(nd, seed = 1)
  v2r <- ta_rotz(20) %*% diag(c(1.2, 0.8, 1.0, 1))
  v2r[1:3, 4] <- c(-14, -9, -11)
  rd <- c(20L, 24L, 18L)
  rv2r <- ta_rotz(-10) %*% diag(c(0.9, 1.1, 1.3, 1))
  rv2r[1:3, 4] <- c(-10, -12, -12)
  M <- diag(4)
  M[1:3, 1:3] <- diag(3) + matrix(rnorm(9, sd = 0.05), 3)
  M[1:3, 4] <- c(1.5, -2, 0.7)

  ref <- apply_transform3d(vol, v2r, M, reference_dim = rd, reference_vox2ras = rv2r,
                           interpolation = "trilinear")
  out <- apply_transform3d_volume(vol, list(M), vox2ras = v2r, reference_dim = rd,
                                  reference_vox2ras = rv2r, interpolation = "trilinear")
  expect_identical(dim(out), rd)
  expect_equal(attr(out, "vox2ras"), rv2r)
  expect_lte(max(abs(out - ref)), 1e-10)

  ref_nn <- apply_transform3d(vol, v2r, M, reference_dim = rd, reference_vox2ras = rv2r,
                              interpolation = "nearest")
  out_nn <- apply_transform3d_volume(vol, list(M), vox2ras = v2r, reference_dim = rd,
                                     reference_vox2ras = rv2r, interpolation = "nearest")
  expect_identical(as.vector(out_nn), as.vector(ref_nn))

  # a bare matrix is a one-element list
  expect_identical(apply_transform3d_volume(vol, M, vox2ras = v2r, reference_dim = rd,
                                            reference_vox2ras = rv2r), out)

  # an all-zero field forces the streaming composite evaluator (no pure-affine
  # shortcut); it must agree with the resampler to the same tolerances
  zf <- array(0, c(rd, 3))
  attr(zf, "vox2ras") <- rv2r
  out_c <- apply_transform3d_volume(vol, list(zf, M), vox2ras = v2r, reference_dim = rd,
                                    reference_vox2ras = rv2r, interpolation = "trilinear")
  expect_lte(max(abs(out_c - ref)), 1e-10)
  out_cn <- apply_transform3d_volume(vol, list(zf, M), vox2ras = v2r, reference_dim = rd,
                                     reference_vox2ras = rv2r, interpolation = "nearest")
  expect_identical(as.vector(out_cn), as.vector(ref_nn))
})

# ---- A2 ---------------------------------------------------------------------

test_that("A2: the forward volume warp replays register_volume3d's warped images", {
  toy <- ta_syn_toy()
  res <- toy$res
  out <- apply_transform3d_volume(toy$moving, res, direction = "forward")
  expect_identical(dim(out), toy$nd)
  expect_equal(attr(out, "vox2ras"), toy$v2r)
  expect_lte(max(abs(out - res$image)), 1e-5 * max(abs(toy$moving)))

  # a nearest-neighbour label channel: <= 0.1% mismatching voxels
  rm <- toy$res_multi
  out_lab <- apply_transform3d_volume(toy$lab_src, rm, direction = "forward",
                                      interpolation = "nearest")
  expect_lte(mean(out_lab != rm$images[[2]]), 0.001)
  expect_setequal(unique(as.vector(out_lab)), c(0, 1, 2))
})

# ---- A3 ---------------------------------------------------------------------

test_that("A3: affine + sinusoidal field composite matches the analytic map", {
  rd <- c(40L, 40L, 40L)
  rv2r <- diag(4)
  rv2r[1:3, 4] <- -20
  sd_nd <- c(56L, 56L, 56L)
  sv2r <- diag(4)
  sv2r[1:3, 4] <- -28
  A <- ta_recenter(ta_rotz(6), rd, rv2r)
  A[1, 1] <- A[1, 1] * 1.05
  A[1:3, 4] <- A[1:3, 4] + c(1.5, -1, 0.5)

  ras <- ta_grid_ras(rd, rv2r)
  amp <- 2
  L <- 25
  ufun <- function(p) {
    cbind(
      amp * sin(2 * pi * p[, 2] / L),
      amp * sin(2 * pi * p[, 3] / L),
      amp * sin(2 * pi * p[, 1] / L)
    )
  }
  u <- ufun(ras)
  field <- array(u, c(rd, 3))
  attr(field, "vox2ras") <- rv2r

  centers <- rbind(c(0, 0, 0), c(6, -4, 3), c(-5, 5, -6), c(3, 7, -2))
  sds <- c(6, 5, 5, 4.5)
  amps <- c(1, .8, .9, .7)
  f <- function(p) {
    val <- numeric(nrow(p))
    for (k in seq_len(nrow(centers))) {
      val <- val + amps[k] * exp(-rowSums(sweep(p, 2, centers[k, ])^2) / (2 * sds[k]^2))
    }
    val
  }
  src <- array(f(ta_grid_ras(sd_nd, sv2r)), sd_nd)

  out <- apply_transform3d_volume(src, list(field, A), vox2ras = sv2r,
                                  reference_dim = rd, reference_vox2ras = rv2r)
  q <- ta_apply_affine(A, ras + u)              # analytic composite A(r + u(r))
  expected <- array(f(q), rd)
  inner <- ta_inner(rd, 3L)
  err <- (out - expected)[inner]
  rng <- diff(range(expected[inner]))
  expect_lte(sqrt(mean(err^2)), 0.01 * rng)
  expect_lte(max(abs(err)), 0.03 * rng)

  # points: exact at the field nodes, <= 0.05 mm between nodes
  set.seed(5)
  idx <- sample(which(inner), 400)
  p <- ras[idx, ]
  expect_lte(max(abs(apply_transform3d_points(p, list(field, A)) - q[idx, ])), 1e-10)
  p_off <- p + matrix(runif(length(p), -0.5, 0.5), ncol = 3)
  expected_off <- ta_apply_affine(A, p_off + ufun(p_off))
  got_off <- apply_transform3d_points(p_off, list(field, A))
  expect_lte(max(ta_norm(got_off - expected_off)), 0.05)

  # list order is ANTs order: the first element is applied first to the point
  got_rev <- apply_transform3d_points(p, list(A, field))
  ap <- ta_apply_affine(A, p)
  expect_lte(max(abs(got_rev - (ap + ta_sample_field(field, ap)))), 1e-10)
  expect_gt(max(abs(got_rev - q[idx, ])), 0.1)
})

# ---- A4 ---------------------------------------------------------------------

test_that("A4: 'inverse' point mapping equals the sample positions of the 'forward' volume warp", {
  toy <- ta_syn_toy()
  res <- toy$res
  ras <- ta_grid_ras(toy$nd, toy$v2r)
  p_src <- apply_transform3d_points(ras, res, direction = "inverse")
  expect_identical(dim(p_src), c(nrow(ras), 3L))
  expect_type(p_src, "double")
  u <- matrix(res$forward_field, ncol = 3)
  expect_lte(max(abs(p_src - ta_apply_affine(res$transform, ras + u))), 1e-10)

  vox <- t(solve(toy$v2r) %*% rbind(t(p_src), 1))[, 1:3]
  manual <- array(ta_tril(toy$moving, vox[, 1], vox[, 2], vox[, 3]), toy$nd)
  out <- apply_transform3d_volume(toy$moving, res, direction = "forward")
  expect_lte(max(abs(out - manual)), 1e-10)
})

# ---- A5 ---------------------------------------------------------------------

# Residuals of both compositions at the grid nodes, in voxel units (the grids
# below have unit spacing): forward o inverse, v(y) + u(y + v(y)), and
# inverse o forward, u(x) + v(x + u(x)). A node counts as interior when it and
# the point it is carried to both lie at least `border` voxels inside the grid,
# so the check never reads the field's border zone.
ta_inverse_residuals <- function(fwd, inv, border = 3L) {
  d <- dim(fwd)[1:3]
  v2r <- attr(fwd, "vox2ras")
  ras <- ta_grid_ras(d, v2r)
  u <- matrix(fwd, ncol = 3)
  v <- matrix(inv, ncol = 3)
  inside <- function(p) {
    vox <- t(solve(v2r) %*% rbind(t(p), 1))[, 1:3, drop = FALSE]
    ok <- rep(TRUE, nrow(vox))
    for (k in 1:3) ok <- ok & vox[, k] >= border & vox[, k] <= d[k] - 1 - border
    ok
  }
  node_in <- inside(ras)
  r1 <- v + ta_sample_field(fwd, ras + v)
  r2 <- u + ta_sample_field(inv, ras + u)
  n1 <- ta_norm(r1)[node_in & inside(ras + v)]
  n2 <- ta_norm(r2)[node_in & inside(ras + u)]
  c(mean_fi = mean(n1), max_fi = max(n1), mean_if = mean(n2), max_if = max(n2),
    n_fi = length(n1), n_if = length(n2))
}

test_that("A5: the fixed-point inverse is accurate (and the old negated field is not)", {
  # a smooth sinusoidal field (one period across the grid) with a 3-voxel
  # amplitude per component
  nd <- c(64L, 64L, 64L)
  v2r <- diag(4)
  v2r[1:3, 4] <- -32
  ras <- ta_grid_ras(nd, v2r)
  L <- 64
  amp <- 3
  u <- cbind(amp * sin(2 * pi * ras[, 2] / L),
             amp * sin(2 * pi * ras[, 3] / L),
             amp * sin(2 * pi * ras[, 1] / L))
  expect_gte(max(ta_norm(u)), 3)              # max displacement >= 3 voxels
  field <- array(u, c(nd, 3))
  attr(field, "vox2ras") <- v2r

  inv <- ravetools:::invert_displacement_field(field)
  expect_identical(dim(inv), c(nd, 3L))
  expect_equal(attr(inv, "vox2ras"), v2r)
  r <- ta_inverse_residuals(field, inv, 3L)
  expect_gt(r[["n_fi"]], 0.8 * prod(nd - 6))
  expect_lte(r[["mean_fi"]], 0.01)
  expect_lte(r[["max_fi"]], 0.1)
  expect_lte(r[["mean_if"]], 0.01)
  expect_lte(r[["max_if"]], 0.1)

  # the former first-order inverse (-u) fails the same check
  neg <- -field
  attr(neg, "vox2ras") <- v2r
  r_old <- ta_inverse_residuals(field, neg, 3L)
  expect_gt(r_old[["mean_fi"]], 0.3)
  expect_gt(r_old[["max_fi"]], 0.3)

  # the inverse is deterministic across thread counts
  ravetools_threads(1)
  on.exit(ravetools_threads("auto"), add = TRUE)
  inv1 <- ravetools:::invert_displacement_field(field)
  ravetools_threads("auto")
  expect_identical(as.vector(inv1), as.vector(inv))

  # register_volume3d's inverse_field is such an inverse of forward_field
  toy <- ta_syn_toy()
  res <- toy$res
  rs <- ta_inverse_residuals(res$forward_field, res$inverse_field, 3L)
  expect_gt(rs[["n_if"]], 0.5 * prod(toy$nd - 6))
  expect_lte(rs[["mean_fi"]], 0.01)
  expect_lte(rs[["max_fi"]], 0.1)
  expect_lte(rs[["mean_if"]], 0.01)
  expect_lte(rs[["max_if"]], 0.1)
  # and the negated forward field would not pass it
  neg_syn <- -res$forward_field
  attr(neg_syn, "vox2ras") <- toy$v2r
  expect_gt(ta_inverse_residuals(res$forward_field, neg_syn, 3L)[["max_fi"]], 0.1)
})

# ---- A6 ---------------------------------------------------------------------

test_that("A6: forward then inverse round-trips a smooth image and points", {
  toy <- ta_syn_toy()
  res <- toy$res
  img <- toy$moving
  fwd <- apply_transform3d_volume(img, res, direction = "forward")
  back <- apply_transform3d_volume(fwd, res, direction = "inverse")
  expect_identical(dim(back), toy$nd)
  expect_equal(attr(back, "vox2ras"), toy$v2r)
  inner <- ta_inner(toy$nd, 3L)
  expect_gte(cor(back[inner], img[inner]), 0.99)

  set.seed(8)
  ras <- ta_grid_ras(toy$nd, toy$v2r)
  p <- ras[sample(which(inner), 300), ]
  pt <- apply_transform3d_points(p, res, direction = "forward")   # source -> target
  pb <- apply_transform3d_points(pt, res, direction = "inverse")  # target -> source
  expect_lte(mean(ta_norm(pb - p)), 0.05)
})

# ---- A7 ---------------------------------------------------------------------

test_that("A7: save/load keeps both appliers consistent and persists SourceDim", {
  skip_if_not_installed("freesurferformats")
  toy <- ta_syn_toy()
  res <- toy$res
  d <- tempfile()
  man <- save_registration(res, d)
  expect_true(any(grepl("^SourceDim:", readLines(man))))
  lr <- load_registration(man)
  expect_equal(lr$geometry$source_dim, toy$nd)
  expect_equal(lr$geometry$target_dim, toy$nd)

  same <- function(a, b, tol = 1e-5) {
    expect_identical(dim(a), dim(b))
    expect_lte(max(abs(as.vector(a) - as.vector(b))), tol)
  }
  same(apply_transform3d_volume(toy$moving, lr, direction = "forward"),
       apply_transform3d_volume(toy$moving, res, direction = "forward"))
  same(apply_transform3d_volume(toy$target, lr, direction = "inverse"),
       apply_transform3d_volume(toy$target, res, direction = "inverse"))
  set.seed(9)
  p <- ta_grid_ras(toy$nd, toy$v2r)[sample(prod(toy$nd), 200), ]
  same(apply_transform3d_points(p, lr, direction = "forward"),
       apply_transform3d_points(p, res, direction = "forward"))
  same(apply_transform3d_points(p, lr, direction = "inverse"),
       apply_transform3d_points(p, res, direction = "inverse"))

  # default output grids: target grid for 'forward', source grid for 'inverse'
  nd_s <- c(20L, 24L, 22L)
  v2r_s <- diag(c(1.2, 1, 0.9, 1))
  v2r_s[1:3, 4] <- c(-12, -12, -10)
  src <- apply_transform3d(toy$target, toy$v2r, diag(4), reference_dim = nd_s,
                           reference_vox2ras = v2r_s)
  rig <- register_volume3d(src, toy$target, v2r_s, toy$v2r, type = "rigid", metric = "cc",
                           shrink_factors = c(4, 2), smoothing_sigmas = c(2, 1),
                           iterations = c(30, 20), sampling_rate = 0.5, verbose = FALSE)
  expect_identical(rig$geometry$source_dim, nd_s)
  fw <- apply_transform3d_volume(src, rig, direction = "forward")
  expect_identical(dim(fw), toy$nd)
  expect_equal(attr(fw, "vox2ras"), toy$v2r)
  iv <- apply_transform3d_volume(toy$target, rig, direction = "inverse")
  expect_identical(dim(iv), nd_s)
  expect_equal(attr(iv, "vox2ras"), v2r_s)
  man2 <- save_registration(rig, tempfile())
  lr2 <- load_registration(man2)
  expect_identical(dim(apply_transform3d_volume(toy$target, lr2, direction = "inverse")), nd_s)

  # an older manifest without SourceDim: the inverse direction needs reference_dim
  txt <- readLines(man2)
  man_old <- tempfile(fileext = ".dcf")
  writeLines(txt[!grepl("^SourceDim:", txt)], man_old)
  file.copy(file.path(dirname(man2), setdiff(list.files(dirname(man2)), basename(man2))),
            dirname(man_old))
  lr_old <- load_registration(man_old)
  expect_null(lr_old$geometry$source_dim)
  expect_error(apply_transform3d_volume(toy$target, lr_old, direction = "inverse"),
               "reference_dim")
  iv_old <- apply_transform3d_volume(toy$target, lr_old, direction = "inverse",
                                     reference_dim = nd_s)
  expect_identical(dim(iv_old), nd_s)
  expect_lte(max(abs(iv_old - iv)), 1e-5)
  expect_identical(dim(apply_transform3d_volume(src, lr_old, direction = "forward")), toy$nd)
})

# ---- A8 ---------------------------------------------------------------------

test_that("A8: edge cases (NA points, outside the field, errors, storage modes, 4D)", {
  rd <- c(16L, 14L, 12L)
  rv2r <- diag(4)
  rv2r[1:3, 4] <- c(-8, -7, -6)
  set.seed(21)
  field <- array(rnorm(prod(rd) * 3, sd = 0.5), c(rd, 3))
  attr(field, "vox2ras") <- rv2r
  A <- diag(4)
  A[1:3, 1:3] <- diag(3) + matrix(rnorm(9, sd = 0.05), 3)
  A[1:3, 4] <- c(2, -1, 0.5)

  # NA rows stay NA, the others are mapped
  p <- rbind(c(1, 2, 3), c(NA, 0, 0), c(100, 100, 100), c(-2, 1, 0.5))
  out <- apply_transform3d_points(p, list(field, A))
  expect_true(all(is.na(out[2, ])))
  expect_false(anyNA(out[-2, ]))
  # points outside the field get the affine only
  expect_equal(out[3, ], ta_apply_affine(A, p[3, , drop = FALSE])[1, ])
  # inside the field: affine of the displaced point
  expect_equal(out[1, ], ta_apply_affine(A, p[1, , drop = FALSE] +
                                          ta_sample_field(field, p[1, , drop = FALSE]))[1, ])
  # a data.frame and a length-3 vector are accepted; the output is N x 3 double
  expect_equal(apply_transform3d_points(as.data.frame(p), list(field, A)), out)
  expect_equal(apply_transform3d_points(c(1, 2, 3), list(field, A)), out[1, , drop = FALSE])
  # invert on a matrix uses its inverse; bare matrix == one-element list
  expect_equal(apply_transform3d_points(p[1, ], A, invert = TRUE),
               ta_apply_affine(solve(A), p[1, , drop = FALSE]))
  expect_equal(apply_transform3d_points(p[1, ], list(A, A), invert = c(TRUE, FALSE)),
               p[1, , drop = FALSE])

  # errors
  expect_error(apply_transform3d_points(p, list(field), invert = TRUE), "only be TRUE for")
  expect_error(apply_transform3d_points(p, list(field, A), invert = c(TRUE, FALSE, TRUE)), "length")
  expect_error(apply_transform3d_points(p[, 1:2], list(A)), "3 columns")
  expect_error(apply_transform3d_points(p, list(A), direction = "forward"), "registration object")
  expect_error(apply_transform3d_points(p, list("not a transform")), "mat")
  nf <- field
  attr(nf, "vox2ras") <- NULL
  expect_error(apply_transform3d_points(p, list(nf)), "vox2ras")
  vol <- ta_phantom(rd, seed = 2)
  expect_error(apply_transform3d_volume(vol, list(A), vox2ras = rv2r), "reference_dim")
  expect_error(apply_transform3d_volume(vol, list(A), reference_dim = rd, reference_vox2ras = rv2r),
               "vox2ras")
  expect_error(apply_transform3d_volume(vol, list(A), vox2ras = rv2r, reference_dim = rd,
                                        reference_vox2ras = rv2r, direction = "inverse"),
               "registration object")
  toy <- ta_syn_toy()
  expect_error(apply_transform3d_volume(toy$moving, toy$res, invert = TRUE), "direction")
  expect_error(apply_transform3d_points(p, toy$res, invert = TRUE), "direction")
  expect_error(apply_transform3d_volume(vol, list(field, A), vox2ras = rv2r, reference_dim = rd,
                                        reference_vox2ras = rv2r, invert = TRUE), "only be TRUE for")
  expect_error(apply_transform3d_volume(array(0, c(4, 4)), list(A), vox2ras = rv2r,
                                        reference_dim = rd, reference_vox2ras = rv2r), "3D")

  # integer and logical volumes are accepted (nearest keeps labels intact)
  lab <- array(sample(0:3, prod(rd), TRUE), rd)
  out_d <- apply_transform3d_volume(lab, list(field, A), vox2ras = rv2r, reference_dim = rd,
                                    reference_vox2ras = rv2r, interpolation = "nearest")
  lab_i <- lab
  storage.mode(lab_i) <- "integer"
  out_i <- apply_transform3d_volume(lab_i, list(field, A), vox2ras = rv2r, reference_dim = rd,
                                    reference_vox2ras = rv2r, interpolation = "nearest")
  expect_identical(out_i, out_d)
  expect_true(all(out_i %in% 0:3))
  msk <- lab > 1
  out_l <- apply_transform3d_volume(msk, list(field, A), vox2ras = rv2r, reference_dim = rd,
                                    reference_vox2ras = rv2r, interpolation = "nearest")
  expect_identical(out_l, apply_transform3d_volume(array(as.double(msk), rd), list(field, A),
                                                   vox2ras = rv2r, reference_dim = rd,
                                                   reference_vox2ras = rv2r,
                                                   interpolation = "nearest"))

  # vox2ras may travel as an attribute; na_fill fills out-of-bounds voxels
  attr(vol, "vox2ras") <- rv2r
  big <- diag(4)
  big[1:3, 4] <- c(50, 0, 0)
  out_na <- apply_transform3d_volume(vol, list(big), reference_dim = rd, reference_vox2ras = rv2r,
                                     na_fill = -1)
  expect_true(all(out_na == -1))

  # 4D input: every frame is transformed identically
  v4 <- array(c(vol, vol * 2 + 1), c(rd, 2))
  out4 <- apply_transform3d_volume(v4, list(field, A), vox2ras = rv2r, reference_dim = rd,
                                   reference_vox2ras = rv2r)
  expect_identical(dim(out4), c(rd, 2L))
  out3 <- apply_transform3d_volume(array(vol, rd), list(field, A), vox2ras = rv2r,
                                   reference_dim = rd, reference_vox2ras = rv2r)
  expect_equal(out4[, , , 1], array(out3, rd))
  # the second frame is 2 * frame + 1 inside, na_fill (0) outside
  inside <- apply_transform3d_volume(array(1, rd), list(field, A), vox2ras = rv2r,
                                     reference_dim = rd, reference_vox2ras = rv2r,
                                     interpolation = "trilinear")
  expect_equal(out4[, , , 2], array((2 * out3 + 1) * (inside == 1), rd))
  out4a <- apply_transform3d_volume(v4, list(A), vox2ras = rv2r, reference_dim = rd,
                                    reference_vox2ras = rv2r, interpolation = "nearest")
  expect_identical(dim(out4a), c(rd, 2L))
  expect_identical(out4a[, , , 1], array(apply_transform3d(array(vol, rd), rv2r, A, rd, rv2r,
                                                            interpolation = "nearest"), rd))

  # results do not depend on the thread count
  ref_vol <- apply_transform3d_volume(vol, list(field, A), reference_dim = rd, reference_vox2ras = rv2r)
  ref_pts <- apply_transform3d_points(p, list(field, A))
  ravetools_threads(1)
  on.exit(ravetools_threads("auto"), add = TRUE)
  expect_identical(apply_transform3d_volume(vol, list(field, A), reference_dim = rd,
                                            reference_vox2ras = rv2r), ref_vol)
  expect_identical(apply_transform3d_points(p, list(field, A)), ref_pts)
  ravetools_threads("auto")
})

test_that("A8: file paths are accepted in a transform list", {
  skip_if_not_installed("freesurferformats")
  rd <- c(12L, 10L, 11L)
  rv2r <- diag(c(1.1, 0.9, 1, 1))
  rv2r[1:3, 4] <- c(-6, -5, -5)
  set.seed(31)
  field <- array(rnorm(prod(rd) * 3, sd = 0.4), c(rd, 3))
  attr(field, "vox2ras") <- rv2r
  A <- diag(4)
  A[1:3, 1:3] <- diag(3) + matrix(rnorm(9, sd = 0.05), 3)
  A[1:3, 4] <- c(1, -2, 0.5)
  fm <- tempfile(fileext = ".mat")
  fw <- tempfile(fileext = ".nii.gz")
  write_ants_transform(A, fm)
  write_ants_warp(field, fw)
  p <- ta_grid_ras(rd, rv2r)[seq(1, prod(rd), by = 7), ]
  mem <- apply_transform3d_points(p, list(field, A))
  disk <- apply_transform3d_points(p, c(fw, fm))
  expect_lte(max(abs(disk - mem)), 1e-5)
  expect_lte(max(abs(apply_transform3d_points(p, list(fw, A), invert = c(FALSE, TRUE)) -
                       apply_transform3d_points(p, list(field, solve(A))))), 1e-5)
  vol <- ta_phantom(rd, seed = 4)
  expect_lte(max(abs(
    apply_transform3d_volume(vol, list(fw, fm), vox2ras = rv2r, reference_dim = rd,
                             reference_vox2ras = rv2r) -
      apply_transform3d_volume(vol, list(field, A), vox2ras = rv2r, reference_dim = rd,
                               reference_vox2ras = rv2r))), 1e-5)
})

# ---- A9 ---------------------------------------------------------------------

test_that("A9: register_volume3d no longer calls apply_transform3d and keeps its outputs", {
  body_txt <- paste(deparse(body(register_volume3d)), collapse = "\n")
  expect_false(grepl("apply_transform3d(", body_txt, fixed = TRUE))
  expect_true(grepl("apply_transform3d_volume(", body_txt, fixed = TRUE))

  # rigid / affine warped outputs still come from the same resampler (bit-identical
  # to apply_transform3d), per channel with its own interpolation
  nd <- c(24L, 24L, 24L)
  v2r <- diag(4)
  v2r[1:3, 4] <- -12
  target <- ta_phantom(nd, seed = 6)
  Ttrue <- diag(4)
  Ttrue[1:3, 4] <- c(1.5, -1, 0.5)
  source <- apply_transform3d(target, v2r, solve(Ttrue), reference_dim = nd, reference_vox2ras = v2r)
  expect_warning(res <- register_volume3d(
    list(source, source^2, source > 0.4), list(target, target^2, target > 0.4), v2r, v2r,
    type = "rigid", metric = "cc", interpolation = c("trilinear", "bspline", "nearest"),
    shrink_factors = c(2, 1), smoothing_sigmas = c(1, 0), iterations = c(30, 20),
    sampling_rate = 0.5, verbose = FALSE), "multiple channels")
  expect_length(res$images, 3L)
  for (k in 1:3) {
    ref <- apply_transform3d(list(source, source^2, source > 0.4)[[k]], v2r, res$transform,
                             reference_dim = nd, reference_vox2ras = v2r,
                             interpolation = c("trilinear", "bspline", "nearest")[k])
    expect_identical(as.vector(res$images[[k]]), as.vector(ref))
  }
  expect_identical(res$image, res$images[[1]])
  expect_identical(res$geometry$source_dim, nd)
  expect_identical(res$geometry$target_dim, nd)
})

# ---- B1-B3 ------------------------------------------------------------------

test_that("B1-B3: SyN controls (explicit defaults, non-default values, validation)", {
  toy <- ta_syn_toy()
  args <- toy$args
  # B1: explicit defaults are bit-identical to the implicit call
  res_b1 <- do.call(register_volume3d, c(args, list(syn_grad_step = 0.2, syn_total_sigma = 0,
                                                    syn_cc_radius = 2L)))
  expect_identical(res_b1$forward_field, toy$res$forward_field)
  expect_identical(res_b1$inverse_field, toy$res$inverse_field)
  expect_identical(res_b1$image, toy$res$image)
  expect_identical(res_b1$transform, toy$res$transform)

  # both fields are identical with a single thread
  ravetools_threads(1)
  on.exit(ravetools_threads("auto"), add = TRUE)
  res_t1 <- do.call(register_volume3d, args)
  ravetools_threads("auto")
  expect_identical(res_t1$forward_field, toy$res$forward_field)
  expect_identical(res_t1$inverse_field, toy$res$inverse_field)

  # B2: non-default values change the field and still recover the deformation
  res_b2 <- do.call(register_volume3d, c(args, list(syn_grad_step = 0.25, syn_total_sigma = 0.5,
                                                    syn_cc_radius = 3L)))
  expect_false(identical(res_b2$forward_field, toy$res$forward_field))
  expect_true(all(is.finite(res_b2$forward_field)))
  expect_gt(cor(as.vector(res_b2$image), as.vector(toy$target)), toy$base_cor + 0.01)
  expect_lt(mean((res_b2$image - toy$target)^2), mean((toy$moving - toy$target)^2) / 2)

  # B3: invalid values error
  bad <- function(...) do.call(register_volume3d, c(args, list(...)))
  expect_error(bad(syn_grad_step = 0), "syn_grad_step")
  expect_error(bad(syn_grad_step = -0.1), "syn_grad_step")
  expect_error(bad(syn_grad_step = NA), "syn_grad_step")
  expect_error(bad(syn_grad_step = c(0.1, 0.2)), "syn_grad_step")
  expect_error(bad(syn_total_sigma = -1), "syn_total_sigma")
  expect_error(bad(syn_total_sigma = NA), "syn_total_sigma")
  expect_error(bad(syn_cc_radius = 0), "syn_cc_radius")
  expect_error(bad(syn_cc_radius = 1.5), "syn_cc_radius")
  expect_error(bad(syn_cc_radius = "a"), "syn_cc_radius")
})

# ---- review fixes ------------------------------------------------------------

test_that("fields use the edge value in the half-voxel margin on both sides (ITK rule)", {
  nd <- c(8L, 4L, 4L)
  v2r <- diag(4)
  field <- array(0, c(nd, 3))
  field[, , , 1] <- 1 + 2 * (0:(nd[1] - 1))            # u_x = 1 + 2 i along x
  attr(field, "vox2ras") <- v2r
  xs <- c(-0.5, -0.49, -0.25, -0.1, 0, 0.5, 6.5, 7, 7.25, 7.49)
  u <- apply_transform3d_points(cbind(xs, 1, 1), list(field))[, 1] - xs
  expect_equal(u[xs < 0], rep(1, sum(xs < 0)), tolerance = 1e-12)      # node 0, no extrapolation
  expect_equal(u[xs > 7], rep(15, sum(xs > 7)), tolerance = 1e-12)     # node 7
  expect_equal(u[xs %in% c(0, 0.5, 6.5, 7)], c(1, 2, 14, 15), tolerance = 1e-12)
  # beyond the margin the field contributes nothing
  q <- apply_transform3d_points(rbind(c(-0.51, 1, 1), c(7.5, 1, 1)), list(field))
  expect_equal(q[, 1], c(-0.51, 7.5))
})

test_that("singular or non-finite matrices are rejected, never silently ignored", {
  nd <- c(8L, 6L, 5L)
  v2r <- diag(4)
  field <- array(0.3, c(nd, 3))
  sing <- diag(4)
  sing[1, 1] <- 0
  attr(field, "vox2ras") <- sing
  p <- rbind(c(1, 1, 1), c(2, 3, 1))
  expect_error(apply_transform3d_points(p, list(field)), "invertible")
  vol <- array(seq_len(prod(nd)) / prod(nd), nd)
  expect_error(apply_transform3d_volume(vol, list(diag(4)), vox2ras = sing,
                                        reference_dim = nd, reference_vox2ras = v2r),
               "invertible")
  expect_error(apply_transform3d_volume(vol, list(diag(4)), vox2ras = v2r,
                                        reference_dim = nd, reference_vox2ras = sing),
               "invertible")
  expect_error(apply_transform3d_points(p, list(sing)), "invertible")
  bad <- diag(4)
  bad[1, 4] <- Inf
  expect_error(apply_transform3d_points(p, list(bad)), "finite")
})

test_that("non-finite point rows give NA rows on the affine and the field path", {
  nd <- c(8L, 6L, 5L)
  field <- array(0.3, c(nd, 3))
  attr(field, "vox2ras") <- diag(4)
  A <- diag(4)
  A[1:3, 4] <- c(1, 2, 3)
  pin <- rbind(c(Inf, 1, 1), c(1, NA, 1), c(1, 1, -Inf), c(2, 1, 1))
  for (tr in list(list(A), list(field, A))) {
    out <- apply_transform3d_points(pin, tr)
    expect_true(all(is.na(out[1:3, ])))
    expect_true(all(is.finite(out[4, ])))
  }
})

test_that("a 5D stack with a trailing singleton is treated as 4D frames", {
  nd <- c(8L, 6L, 5L)
  v2r <- diag(4)
  A <- diag(4)
  A[1:3, 4] <- c(0.5, -0.25, 0.75)
  set.seed(5)
  vol4 <- array(stats::runif(prod(nd) * 2), c(nd, 2))
  vol5 <- vol4
  dim(vol5) <- c(nd, 2L, 1L)
  run <- function(v) {
    apply_transform3d_volume(
      v,
      list(A),
      vox2ras = v2r,
      reference_dim = nd,
      reference_vox2ras = v2r
    )
  }
  expect_identical(run(vol5), run(vol4))
  expect_identical(dim(run(vol5)), c(nd, 2L))
  vol6 <- vol4
  dim(vol6) <- c(nd, 1L, 2L)
  expect_error(run(vol6), "4D array of frames")
})

test_that("register_volume3d inverts its SyN field with at most 20 sweeps", {
  nd <- c(16L, 16L, 16L)
  v2r <- diag(4)
  v2r[1:3, 4] <- -8
  target <- ta_phantom(nd, seed = 5, nblob = 3)
  moving <- ta_phantom(nd, seed = 6, nblob = 3)
  out <- capture.output(res <- register_volume3d(
    moving, target, v2r, v2r, type = "syn", metric = "cc",
    shrink_factors = c(2, 1), smoothing_sigmas = c(1, 0), iterations = c(20, 10),
    syn_iterations = c(10, 5), verbose = TRUE))
  line <- grep("inverse field:", out, value = TRUE)
  expect_length(line, 1L)
  sweeps <- as.integer(sub(".*inverse field: ([0-9]+) fixed-point.*", "\\1", line))
  expect_lte(sweeps, 20L)
  expect_true(all(is.finite(res$inverse_field)))
})
