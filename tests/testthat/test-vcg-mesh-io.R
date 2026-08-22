
# Regression coverage for the functions that load a mesh through
# `IOMesh::vcgReadR` (src/vcgCommon.h). `vcgReadR` zeroes every vertex normal
# it does not read from R; before that fix the normals of vertices belonging to
# no face were left uninitialised, which valgrind flagged inside
# tri::MidPoint during vcgEdgeLengthSubdivision.
#
# Vertices referenced by a face always had their normals recomputed by
# UpdateNormal::PerVertexClear, so isolated vertices are the only place where
# behaviour could differ - the tests below pin that down, and give the rest of
# the vcgReadR call sites enough coverage to catch a future change there.

# Sphere plus `n` vertices that belong to no face.
sphere_with_isolated <- function(n = 2L, sub_division = 2L) {
  mesh <- vcg_sphere(sub_division = sub_division)
  extra <- matrix(c(5, 5, 5, -5, -5, -5), nrow = 3L)[, seq_len(n), drop = FALSE]
  mesh$vb <- cbind(mesh$vb, rbind(extra, 1))
  mesh$normals <- NULL
  mesh
}

# Column indices of the isolated vertices appended by sphere_with_isolated()
isolated_columns <- function(mesh, n = 2L) {
  seq.int(ncol(mesh$vb) - n + 1L, ncol(mesh$vb))
}

# ---- vcg_update_normals ------------------------------------------------

test_that("vcg_update_normals leaves isolated vertices with exactly zero normals", {
  mesh <- sphere_with_isolated()
  iso  <- isolated_columns(mesh)

  for (weight in c("area", "angle")) {
    normals <- vcg_update_normals(mesh, weight = weight)$normals
    expect_equal(ncol(normals), ncol(mesh$vb))
    expect_true(all(is.finite(normals)))
    # vertices belonging to no face are never touched by PerVertexClear
    expect_true(all(normals[1:3, iso] == 0))
  }
})

test_that("vcg_update_normals returns unit normals for face-referenced vertices", {
  mesh <- sphere_with_isolated()
  iso  <- isolated_columns(mesh)
  normals <- vcg_update_normals(mesh)$normals
  lengths <- sqrt(colSums(normals[1:3, -iso, drop = FALSE]^2))
  expect_equal(lengths, rep(1, length(lengths)), tolerance = 1e-5)
})

test_that("vcg_update_normals is deterministic", {
  mesh <- sphere_with_isolated()
  expect_identical(
    vcg_update_normals(mesh)$normals,
    vcg_update_normals(mesh)$normals
  )
})

test_that("vcg_update_normals handles a point cloud with no faces", {
  # every vertex is isolated here, so this takes the PointCloudNormal branch
  set.seed(1)
  points  <- matrix(rnorm(60), ncol = 3)
  normals <- vcg_update_normals(points)$normals

  expect_equal(ncol(normals), nrow(points))
  expect_true(all(is.finite(normals)))
  lengths <- sqrt(colSums(normals[1:3, , drop = FALSE]^2))
  expect_equal(lengths, rep(1, length(lengths)), tolerance = 1e-5)
  expect_identical(normals, vcg_update_normals(points)$normals)
})

# ---- vcg_smooth_explicit -----------------------------------------------

test_that("vcg_smooth_explicit keeps isolated vertices intact for every type", {
  mesh <- sphere_with_isolated()
  iso  <- isolated_columns(mesh)
  types <- c("taubin", "laplace", "HClaplace", "fujiLaplace",
             "angWeight", "surfPreserveLaplace")

  for (type in types) {
    smoothed <- vcg_smooth_explicit(mesh, type = type, iteration = 2)
    expect_equal(ncol(smoothed$vb), ncol(mesh$vb))
    expect_true(all(is.finite(smoothed$normals)))
    expect_true(all(smoothed$normals[1:3, iso] == 0))
    # smoothing moves vertices along incident faces; isolated ones have none
    expect_equal(smoothed$vb[1:3, iso], mesh$vb[1:3, iso])
  }
})

# ---- vcg_smooth_implicit -----------------------------------------------

# NOTE: vcg_smooth_implicit is deliberately exercised without isolated
# vertices. An unreferenced vertex makes the implicit solve singular and the
# whole mesh comes back as garbage - all zeros, denormals or NaN depending on
# heap state. That is a pre-existing defect (it reproduces before the vcgReadR
# normal fix and implicit_smooth.h never reads a vertex normal), so there is no
# defined behaviour to assert on. Do not "fix" this test by pinning down
# whatever the solver happens to return for that input.

test_that("vcg_smooth_implicit shrinks a sphere without distorting it", {
  sphere   <- vcg_sphere(sub_division = 2)
  smoothed <- vcg_smooth_implicit(sphere)

  expect_equal(ncol(smoothed$vb), ncol(sphere$vb))
  expect_equal(ncol(smoothed$it), ncol(sphere$it))
  expect_true(all(is.finite(smoothed$vb)))
  expect_true(all(is.finite(smoothed$normals)))

  # smoothing a sphere keeps it a sphere, slightly contracted
  radii <- sqrt(colSums(smoothed$vb[1:3, ]^2))
  expect_true(all(radii > 0.5 & radii < 1))
  expect_lt(diff(range(radii)), 0.05)

  lengths <- sqrt(colSums(smoothed$normals[1:3, , drop = FALSE]^2))
  expect_equal(lengths, rep(1, length(lengths)), tolerance = 1e-5)
})

# ---- vcg_fix_defects ---------------------------------------------------

test_that("vcg_fix_defects fills a hole and drops unreferenced vertices", {
  mesh <- sphere_with_isolated()
  mesh$it <- mesh$it[, -1, drop = FALSE]   # punch a triangular hole

  repaired <- vcg_fix_defects(mesh)
  info     <- attr(repaired, "info")

  expect_s3_class(repaired, "mesh3d")
  expect_gte(info$holes_filled, 1L)
  expect_equal(info$boundary_edges_after, 0L)
  expect_equal(info$nonmanifold_edges_after, 0L)

  # unlike the smoothers, vcg_fix_defects compacts unreferenced vertices away
  expect_equal(ncol(repaired$vb), ncol(mesh$vb) - 2L)
  expect_true(all(is.finite(repaired$normals)))
  lengths <- sqrt(colSums(repaired$normals[1:3, , drop = FALSE]^2))
  expect_equal(lengths, rep(1, length(lengths)), tolerance = 1e-5)
})

# ---- vcg_subdivision ---------------------------------------------------

test_that("vcg_subdivision edge method splits every face into four", {
  sphere <- vcg_sphere(sub_division = 2)
  result <- vcg_subdivision(sphere, method = "edge")

  expect_s3_class(result, "mesh3d")
  expect_equal(ncol(result$it), ncol(sphere$it) * 4L)
  expect_equal(nrow(result$vb), 4L)
  expect_true(all(is.finite(result$vb)))
  expect_true(all(result$it >= 1L & result$it <= ncol(result$vb)))
})

test_that("vcg_subdivision edge method only adds original edge midpoints", {
  sphere <- vcg_sphere(sub_division = 2)
  result <- vcg_subdivision(sphere, method = "edge")

  vb <- sphere$vb[1:3, ]
  it <- sphere$it
  midpoints <- cbind(
    (vb[, it[1, ]] + vb[, it[2, ]]) / 2,
    (vb[, it[2, ]] + vb[, it[3, ]]) / 2,
    (vb[, it[3, ]] + vb[, it[1, ]]) / 2
  )
  candidates <- cbind(vb, midpoints)

  nearest <- vcg_kdtree_nearest(
    target = t(candidates), query = t(result$vb[1:3, ]), k = 1
  )
  expect_lt(max(nearest$distance), 1e-6)
})

test_that("vcg_subdivision barycenter method adds one vertex per face", {
  sphere <- vcg_sphere(sub_division = 2)
  result <- vcg_subdivision(sphere, method = "barycenter")

  expect_equal(ncol(result$vb), ncol(sphere$vb) + ncol(sphere$it))
  expect_equal(ncol(result$it), ncol(sphere$it) * 3L)
  expect_true(all(is.finite(result$vb)))
})

# ---- vcg_subdivide_max_edge_length -------------------------------------

test_that("vcg_subdivide_max_edge_length produces finite coordinates", {
  sphere <- vcg_sphere(sub_division = 2)
  result <- vcg_subdivide_max_edge_length(
    sphere, max_edge_len = vcg_max_edge_length(sphere) * 0.4
  )
  expect_true(all(is.finite(result$vb)))
  expect_true(all(result$it >= 1L & result$it <= ncol(result$vb)))
})

test_that("vcg_subdivide_max_edge_length is idempotent at the same threshold", {
  sphere    <- vcg_sphere(sub_division = 2)
  threshold <- vcg_max_edge_length(sphere) * 0.4
  once  <- vcg_subdivide_max_edge_length(sphere, max_edge_len = threshold)
  twice <- vcg_subdivide_max_edge_length(once, max_edge_len = threshold)
  expect_equal(twice$vb, once$vb)
  expect_equal(twice$it, once$it)
})

# ---- vcg_mesh_volume ---------------------------------------------------

test_that("vcg_mesh_volume approximates the volume of a unit sphere", {
  volume <- vcg_mesh_volume(vcg_sphere())
  expect_length(volume, 1L)
  expect_gt(volume, 0)
  # a 642-vertex icosphere is inscribed, so it slightly under-estimates
  expect_equal(volume, 4 / 3 * pi, tolerance = 0.05)
})

# ---- vcg_uniform_remesh ------------------------------------------------

test_that("vcg_uniform_remesh returns a mesh with unit normals", {
  sphere <- vcg_sphere(sub_division = 2)
  result <- vcg_uniform_remesh(sphere, voxel_size = 0.15, verbose = FALSE)

  expect_s3_class(result, "mesh3d")
  expect_gt(ncol(result$it), 0L)
  expect_true(all(is.finite(result$vb)))
  expect_true(all(is.finite(result$normals)))
  lengths <- sqrt(colSums(result$normals[1:3, , drop = FALSE]^2))
  expect_equal(lengths, rep(1, length(lengths)), tolerance = 1e-5)
})

# ---- vcg_raycaster -----------------------------------------------------

test_that("vcg_raycaster hits a unit sphere along the principal axes", {
  sphere <- vcg_sphere()
  origin <- cbind(c(0, 0, 3), c(0, 3, 0), c(3, 0, 0),
                  c(0, 0, -3), c(0, -3, 0), c(-3, 0, 0))
  direction <- cbind(c(0, 0, -1), c(0, -1, 0), c(-1, 0, 0),
                     c(0, 0, 1), c(0, 1, 0), c(1, 0, 0))

  result <- vcg_raycaster(sphere, origin, direction)

  expect_true(all(result$has_intersection))
  # travel from radius 3 to the sphere surface at radius 1
  expect_equal(result$distance, rep(2, 6), tolerance = 1e-3)
  expect_equal(sqrt(colSums(result$intersection^2)), rep(1, 6), tolerance = 1e-3)
  expect_equal(sqrt(colSums(result$normals^2)), rep(1, 6), tolerance = 1e-5)
  # surface normals point outward, against the inbound rays
  expect_true(all(colSums(result$normals * direction) < 0))
  expect_true(all(result$face_index >= 1L & result$face_index <= ncol(sphere$it)))
})

# ---- vcg_subset_vertex -------------------------------------------------

test_that("vcg_subset_vertex keeps only the selected vertices", {
  sphere   <- vcg_sphere(sub_division = 2)
  selector <- rep(FALSE, ncol(sphere$vb))
  selector[sphere$it[1, 1:50]] <- TRUE

  result <- vcg_subset_vertex(sphere, selector)

  expect_equal(ncol(result$vb), sum(selector))
  expect_true(all(result$it >= 1L & result$it <= ncol(result$vb)))
  # every surviving vertex is one of the selected input vertices
  nearest <- vcg_kdtree_nearest(
    target = t(sphere$vb[1:3, selector, drop = FALSE]),
    query  = t(result$vb[1:3, ]),
    k      = 1
  )
  expect_lt(max(nearest$distance), 1e-6)
})

# ---- vcg_count_edge_defects --------------------------------------------

test_that("vcg_count_edge_defects reports a closed sphere as defect free", {
  defects <- vcg_count_edge_defects(vcg_sphere(sub_division = 2))
  expect_equal(defects$boundary_edges, 0L)
  expect_equal(defects$nonmanifold_edges, 0L)
  expect_true(defects$is_closed_manifold)
})

test_that("vcg_count_edge_defects counts the boundary of a single triangle", {
  triangle <- structure(
    list(
      vb = rbind(cbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0)), 1),
      it = matrix(1:3, nrow = 3L)
    ),
    class = c("mesh3d", "shape3d")
  )
  defects <- vcg_count_edge_defects(triangle)
  expect_equal(defects$boundary_edges, 3L)
  expect_equal(defects$nonmanifold_edges, 0L)
  expect_false(defects$is_closed_manifold)
})

# ---- dijkstras_surface_distance ----------------------------------------

test_that("dijkstras_surface_distance measures geodesics on a sphere", {
  sphere <- vcg_sphere(sub_division = 2)
  result <- dijkstras_surface_distance(
    positions  = t(sphere$vb[1:3, ]),
    faces      = t(sphere$it),
    start_node = 1L
  )
  distance <- result$paths$distance

  expect_equal(nrow(result$paths), ncol(sphere$vb))
  expect_equal(distance[1], 0)
  expect_true(all(is.finite(distance)))
  expect_true(all(distance >= 0))
  # the antipode of a unit sphere is pi away; a polyhedral path is never shorter
  expect_lte(max(distance), pi * 1.05)
})
