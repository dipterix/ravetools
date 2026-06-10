
# ---- vcg_max_edge_length -----------------------------------------------

test_that("vcg_max_edge_length returns a positive scalar", {
  sphere <- vcg_sphere()
  mel <- vcg_max_edge_length(sphere)
  expect_length(mel, 1L)
  expect_true(is.finite(mel))
  expect_gt(mel, 0)
})

test_that("vcg_max_edge_length >= vcg_average_edge_length", {
  sphere <- vcg_sphere()
  expect_gte(vcg_max_edge_length(sphere), vcg_average_edge_length(sphere))
})

test_that("vcg_max_edge_length agrees with brute-force calculation", {
  sphere <- vcg_sphere()
  vb <- sphere$vb[1:3, ]
  it <- sphere$it
  e1 <- vb[, it[1, ]] - vb[, it[2, ]]
  e2 <- vb[, it[2, ]] - vb[, it[3, ]]
  e3 <- vb[, it[3, ]] - vb[, it[1, ]]
  expected <- max(sqrt(colSums(e1^2)), sqrt(colSums(e2^2)), sqrt(colSums(e3^2)))
  expect_equal(vcg_max_edge_length(sphere), expected, tolerance = 1e-5)
})

# ---- vcg_subdivide_max_edge_length -------------------------------------

test_that("vcg_subdivide_max_edge_length satisfies the threshold", {
  sphere    <- vcg_sphere()
  cur_max   <- vcg_max_edge_length(sphere)
  threshold <- cur_max * 0.4
  result    <- vcg_subdivide_max_edge_length(sphere, max_edge_len = threshold)
  expect_lte(vcg_max_edge_length(result), threshold * 1.001)
})

test_that("vcg_subdivide_max_edge_length is a no-op when threshold exceeds current max", {
  sphere  <- vcg_sphere()
  cur_max <- vcg_max_edge_length(sphere)
  result  <- vcg_subdivide_max_edge_length(sphere, max_edge_len = cur_max * 2)
  expect_equal(sphere$it, result$it)
  expect_equal(sphere$vb, result$vb)
})

test_that("vcg_subdivide_max_edge_length returns a valid mesh3d", {
  sphere <- vcg_sphere()
  result <- vcg_subdivide_max_edge_length(sphere, max_edge_len = 0.1)
  expect_true(inherits(result, "mesh3d"))
  expect_equal(nrow(result$vb), 4L)
  expect_equal(nrow(result$it), 3L)
  expect_gt(ncol(result$it), ncol(sphere$it))
})

test_that("vcg_subdivide_max_edge_length explicit max_iter=1 may leave some long edges", {
  sphere    <- vcg_sphere()
  cur_max   <- vcg_max_edge_length(sphere)
  threshold <- cur_max * 0.3
  # one pass may not fully satisfy the threshold
  result1   <- vcg_subdivide_max_edge_length(sphere, max_edge_len = threshold, max_iter = 1L)
  result_full <- vcg_subdivide_max_edge_length(sphere, max_edge_len = threshold)
  # full should have at least as many faces as single-pass
  expect_gte(ncol(result_full$it), ncol(result1$it))
})

# ---- vcg_mesh_patch ----------------------------------------------------
# waypoints are now an n x 3 coordinate matrix (x, y, z per row).
# Construct by transposing columns of sphere$vb.

test_that("vcg_mesh_patch returns a length-2 list of mesh3d objects", {
  sphere <- vcg_sphere()
  n      <- ncol(sphere$vb)
  wp_idx <- c(1L, as.integer(n * 0.25), as.integer(n * 0.5), as.integer(n * 0.75))
  wp     <- t(sphere$vb[1:3, wp_idx, drop = FALSE])
  result <- vcg_mesh_patch(sphere, waypoints = wp)
  expect_length(result, 2L)
  expect_true(inherits(result[[1]], "mesh3d"))
  expect_true(inherits(result[[2]], "mesh3d"))
})

test_that("vcg_mesh_patch face counts partition the original mesh", {
  sphere <- vcg_sphere()
  n      <- ncol(sphere$vb)
  wp_idx <- c(1L, as.integer(n * 0.25), as.integer(n * 0.5), as.integer(n * 0.75))
  wp     <- t(sphere$vb[1:3, wp_idx, drop = FALSE])
  result <- vcg_mesh_patch(sphere, waypoints = wp)
  expect_equal(ncol(result[[1]]$it) + ncol(result[[2]]$it), ncol(sphere$it))
})

test_that("vcg_mesh_patch orig_vertex maps correctly to original vertices", {
  sphere <- vcg_sphere()
  n      <- ncol(sphere$vb)
  wp_idx <- c(1L, as.integer(n * 0.25), as.integer(n * 0.5), as.integer(n * 0.75))
  wp     <- t(sphere$vb[1:3, wp_idx, drop = FALSE])
  result <- vcg_mesh_patch(sphere, waypoints = wp)
  p1 <- result[[1]]
  p2 <- result[[2]]
  expect_equal(p1$vb[1:3, ], sphere$vb[1:3, p1$orig_vertex], tolerance = 1e-6)
  expect_equal(p2$vb[1:3, ], sphere$vb[1:3, p2$orig_vertex], tolerance = 1e-6)
})

test_that("vcg_mesh_patch seed_vertex selects the correct side", {
  sphere <- vcg_sphere()
  n      <- ncol(sphere$vb)
  wp_idx <- c(1L, as.integer(n * 0.25), as.integer(n * 0.5), as.integer(n * 0.75))
  wp     <- t(sphere$vb[1:3, wp_idx, drop = FALSE])

  result_auto <- vcg_mesh_patch(sphere, waypoints = wp)
  p1_auto <- result_auto[[1]]

  # A vertex from p1_auto should stay in the first patch
  inside_v <- p1_auto$orig_vertex[ceiling(length(p1_auto$orig_vertex) / 2L)]
  result_seed <- vcg_mesh_patch(sphere, waypoints = wp, seed_vertex = inside_v)
  p1_seed <- result_seed[[1]]

  expect_setequal(p1_auto$orig_vertex, p1_seed$orig_vertex)
})

test_that("vcg_mesh_patch p1 centroid is closer to mean waypoint than p2", {
  sphere <- vcg_sphere()
  n      <- ncol(sphere$vb)
  wp_idx <- c(1L, as.integer(n * 0.25), as.integer(n * 0.5), as.integer(n * 0.75))
  wp     <- t(sphere$vb[1:3, wp_idx, drop = FALSE])
  result  <- vcg_mesh_patch(sphere, waypoints = wp)
  p1 <- result[[1]]
  p2 <- result[[2]]
  mean_wp <- colMeans(wp)   # mean of the 4 waypoint coordinates
  d1 <- sqrt(sum((rowMeans(p1$vb[1:3, , drop = FALSE]) - mean_wp)^2))
  d2 <- sqrt(sum((rowMeans(p2$vb[1:3, , drop = FALSE]) - mean_wp)^2))
  expect_lte(d1, d2)
})

test_that("vcg_mesh_patch errors with fewer than 3 waypoints", {
  sphere <- vcg_sphere()
  wp2 <- t(sphere$vb[1:3, 1:2, drop = FALSE])   # 2-row matrix
  expect_error(
    vcg_mesh_patch(sphere, waypoints = wp2),
    "at least 3"
  )
})

test_that("vcg_mesh_patch errors with wrong waypoints column count", {
  sphere <- vcg_sphere()
  expect_error(
    vcg_mesh_patch(sphere, waypoints = matrix(1, nrow = 4, ncol = 2)),
    "3 columns"
  )
})

test_that("vcg_mesh_patch works with 3 waypoints", {
  sphere <- vcg_sphere()
  n      <- ncol(sphere$vb)
  wp_idx <- c(1L, as.integer(n * 0.33), as.integer(n * 0.66))
  wp     <- t(sphere$vb[1:3, wp_idx, drop = FALSE])
  result <- vcg_mesh_patch(sphere, waypoints = wp)
  expect_equal(ncol(result[[1]]$it) + ncol(result[[2]]$it), ncol(sphere$it))
})

test_that("vcg_mesh_patch max_edge_length refines the output patches", {
  sphere <- vcg_sphere()
  n      <- ncol(sphere$vb)
  wp_idx <- c(1L, as.integer(n * 0.25), as.integer(n * 0.5), as.integer(n * 0.75))
  wp     <- t(sphere$vb[1:3, wp_idx, drop = FALSE])
  mel    <- vcg_max_edge_length(sphere) * 0.5

  result <- vcg_mesh_patch(sphere, waypoints = wp, max_edge_length = mel)
  p1 <- result[[1]]
  p2 <- result[[2]]

  # Both patches should respect the max_edge_length
  expect_lte(vcg_max_edge_length(p1), mel * 1.001)
  expect_lte(vcg_max_edge_length(p2), mel * 1.001)

  # More faces than without refinement
  result_plain <- vcg_mesh_patch(sphere, waypoints = wp)
  expect_gt(ncol(p1$it) + ncol(p2$it),
            ncol(result_plain[[1]]$it) + ncol(result_plain[[2]]$it))

  # orig_vertex values are valid indices
  expect_true(all(p1$orig_vertex >= 1L))
  expect_true(all(p2$orig_vertex >= 1L))
})
