
# Closed box surfaces from a binary volume, 1 mm apart. Each element of
# `boxes` lists the occupied voxel indices along x, y, and z; the surface is
# shifted by 0.3 mm so no face lies exactly halfway between lattice points
box_surface <- function(boxes, dim = c(80L, 80L, 80L)) {
  volume <- array(0, dim)
  for (b in boxes) {
    volume[b$x, b$y, b$z] <- 1
  }
  mesh <- vcg_isosurface(volume, threshold_lb = 0.5, vox_to_ras = diag(4))
  mesh$vb[1:3, ] <- mesh$vb[1:3, ] + 0.3
  mesh
}

# Brute-force binary closing with the digital disk/ball of `radius`; points
# outside the array count as outside for the dilation, inside for the erosion
brute_force_close <- function(mask, radius, slice_axis) {
  d <- dim(mask)
  r <- ceiling(radius)
  offsets <- as.matrix(expand.grid(-r:r, -r:r, -r:r))
  if (slice_axis >= 0) {
    offsets <- offsets[offsets[, slice_axis + 1] == 0, , drop = FALSE]
  }
  offsets <- offsets[rowSums(offsets^2) <= radius^2, , drop = FALSE]
  points <- as.matrix(expand.grid(seq_len(d[1]), seq_len(d[2]), seq_len(d[3])))
  value_at <- function(m, offset, outside) {
    p <- sweep(points, 2, offset, "+")
    ok <- p[, 1] >= 1 & p[, 1] <= d[1] & p[, 2] >= 1 & p[, 2] <= d[2] &
      p[, 3] >= 1 & p[, 3] <= d[3]
    v <- rep(outside, nrow(points))
    v[ok] <- m[p[ok, , drop = FALSE]]
    v
  }
  dilated <- array(FALSE, d)
  for (i in seq_len(nrow(offsets))) {
    dilated <- dilated | value_at(mask, offsets[i, ], FALSE)
  }
  closed <- array(TRUE, d)
  for (i in seq_len(nrow(offsets))) {
    closed <- closed & value_at(dilated, offsets[i, ], TRUE)
  }
  closed
}

test_that("mrisCloseBinary matches a brute-force closing", {
  set.seed(42)
  dims <- c(14L, 12L, 10L)
  mask <- array(stats::runif(prod(dims)) > 0.75, dims)
  for (slice_axis in c(-1L, 0L, 1L, 2L)) {
    for (radius in c(0, 1, 1.5, 2.5)) {
      expect_identical(
        c(ravetools:::mrisCloseBinary(mask, dims, radius, slice_axis)),
        c(brute_force_close(mask, radius, slice_axis)),
        info = sprintf("slice_axis = %d, radius = %g", slice_axis, radius)
      )
    }
  }
})

test_that("mrisFillInteriorParity fills the lattice points inside a sphere", {
  sphere <- vcg_sphere(sub_division = 5L)
  vb <- sphere$vb[1:3, ] * 10 + 12.3
  dims <- c(26L, 26L, 26L)
  filled <- ravetools:::mrisFillInteriorParity(vb, sphere$it, dims)
  expect_identical(attr(filled, "odd_columns"), 0L)
  expect_identical(dim(filled), dims)

  grid <- as.matrix(expand.grid(0:25, 0:25, 0:25))
  r <- sqrt(rowSums(sweep(grid, 2, 12.3)^2))
  # the polyhedron is inscribed in the sphere, its faces less than 0.01 inside
  expect_true(all(filled[r < 9.95]))
  expect_false(any(filled[r > 10]))
})

test_that("mrisLargestComponent keeps the biggest connected piece", {
  two <- box_surface(list(
    list(x = 10:19, y = 10:19, z = 10:19),
    list(x = 30:49, y = 10:29, z = 10:29)
  ))
  components <- ravetools:::mrisLargestComponent(two$it, ncol(two$vb))
  expect_identical(components$n_components, 2L)
  # only the larger box (x from 28.8) is kept
  expect_gt(min(two$vb[1, components$keep]), 25)
  expect_lt(max(two$vb[1, !components$keep]), 25)
})

test_that("mris_outer_surface leaves a convex surface nearly unchanged", {
  sphere <- vcg_sphere(sub_division = 4L)
  sphere$vb[1:3, ] <- sphere$vb[1:3, ] * 40
  for (method in c("freesurfer", "isotropic")) {
    envelope <- mris_outer_surface(sphere, method = method)
    info <- attr(envelope, "info")
    radius <- sqrt(colSums(envelope$vb[1:3, ]^2))
    expect_lt(abs(mean(radius) - 40), 1)
    expect_lt(max(abs(radius - 40)), 2)
    expect_true(info$is_closed_manifold)
    expect_equal(info$euler, 2)
    expect_equal(info$n_components, 1)
    expect_equal(info$odd_columns, 0)
  }
})

test_that("mris_outer_surface bridges gaps narrower than the closing element", {
  # 6 mm apart along x, which lies in the coronal plane: bridged by both methods
  near <- box_surface(list(
    list(x = 10:29, y = 10:29, z = 10:29),
    list(x = 36:55, y = 10:29, z = 10:29)
  ))
  # 40 mm apart: wider than the 31-voxel disk, so not bridged
  far <- box_surface(list(
    list(x = 3:18, y = 10:29, z = 10:29),
    list(x = 59:74, y = 10:29, z = 10:29)
  ))
  for (method in c("freesurfer", "isotropic")) {
    expect_equal(attr(mris_outer_surface(near, method = method), "info")$n_components, 1)
    expect_equal(attr(mris_outer_surface(far, method = method), "info")$n_components, 2)
  }
  # no closing: the boxes stay apart
  expect_equal(
    attr(mris_outer_surface(near, close_radius = 0), "info")$n_components, 2)
})

test_that("the freesurfer method closes within coronal slices only", {
  # boxes stacked along y (anterior-posterior), 3 mm apart: the slit lies in a
  # coronal plane, which the per-slice disk never crosses
  slit <- box_surface(list(
    list(x = 10:39, y = 10:24, z = 10:39),
    list(x = 10:39, y = 28:42, z = 10:39)
  ))
  freesurfer <- mris_outer_surface(slit, method = "freesurfer")
  isotropic <- mris_outer_surface(slit, method = "isotropic")
  expect_equal(attr(freesurfer, "info")$n_components, 2)
  expect_equal(attr(isotropic, "info")$n_components, 1)
})

test_that("the freesurfer one-voxel growth goes toward +x and +z only", {
  mask <- array(FALSE, c(5L, 5L, 5L))
  mask[3, 3, 3] <- TRUE
  grown <- ravetools:::outer_surface_grow_right_superior(mask)
  expected <- array(FALSE, c(5L, 5L, 5L))
  expected[3:4, 3, 3:4] <- TRUE
  expect_identical(grown, expected)
})

test_that("mris_outer_surface is deterministic and keeps the input coordinates", {
  box <- box_surface(list(list(x = 20:49, y = 20:49, z = 20:49)))
  a <- mris_outer_surface(box)
  b <- mris_outer_surface(box)
  expect_identical(a, b)
  expect_s3_class(a, "mesh3d")
  # within the one-voxel growth of the freesurfer method
  expect_equal(rowMeans(a$vb[1:3, ]), rowMeans(box$vb[1:3, ]), tolerance = 1,
               scale = 1)
})

test_that("mris_outer_surface validates its arguments", {
  sphere <- vcg_sphere()
  sphere$vb[1:3, ] <- sphere$vb[1:3, ] * 20
  expect_error(mris_outer_surface(sphere, close_radius = -1), "close_radius")
  expect_error(mris_outer_surface(sphere, close_radius = NA), "close_radius")
  expect_error(mris_outer_surface(sphere, smooth_iterations = -1), "smooth_iterations")
  expect_error(mris_outer_surface(sphere, method = "convex"))
  expect_error(mris_outer_surface(list(vb = sphere$vb)), "face")
})
