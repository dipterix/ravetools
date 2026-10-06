
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

test_that("vcg_smooth_implicit keeps isolated vertices intact", {
  mesh <- sphere_with_isolated()
  iso  <- isolated_columns(mesh)
  smoothed <- vcg_smooth_implicit(mesh)

  expect_equal(ncol(smoothed$vb), ncol(mesh$vb))
  expect_true(all(is.finite(smoothed$vb)))
  expect_true(all(is.finite(smoothed$normals)))
  expect_true(all(smoothed$normals[1:3, iso] == 0))
  # a vertex with no incident face has nothing to smooth against
  expect_equal(smoothed$vb[1:3, iso], mesh$vb[1:3, iso])
  expect_identical(smoothed$it, mesh$it)
})

test_that("vcg_smooth_implicit ignores unreferenced vertices entirely", {
  # an unreferenced vertex used to make the implicit solve singular, which
  # destroyed the whole mesh rather than just that vertex
  clean  <- vcg_sphere(sub_division = 2)
  padded <- sphere_with_isolated()
  kept   <- seq_len(ncol(clean$vb))

  expect_equal(
    vcg_smooth_implicit(padded)$vb[1:3, kept],
    vcg_smooth_implicit(clean)$vb[1:3, ]
  )
})

test_that("vcg_smooth_implicit is deterministic with unreferenced vertices", {
  mesh <- sphere_with_isolated()
  expect_identical(vcg_smooth_implicit(mesh)$vb, vcg_smooth_implicit(mesh)$vb)
})

test_that("vcg_smooth_implicit returns a face-less mesh unchanged", {
  cloud <- vcg_sphere(sub_division = 2)
  cloud$it <- NULL
  cloud$normals <- NULL
  smoothed <- vcg_smooth_implicit(cloud)

  expect_equal(smoothed$vb[1:3, ], cloud$vb[1:3, ])
  expect_true(all(smoothed$normals[1:3, ] == 0))
})

# The tests below derive the expected result straight from the definition of
# implicit smoothing, independently of the solver in src/vcgCommon.cpp:
# solve (M + lambda * L^(2^(degree - 1))) X = M V, with M the per-vertex sum
# of doubled incident face areas scaled by its maximum, and L the Laplacian
# accumulated per face edge (weight `laplacian_weight`, or half the cotangent
# of the angle opposite the edge).

# Cap of a sphere: the faces whose centroid lies above z = 0.5 are dropped,
# leaving one open boundary loop. Scaling by 1.1 makes the coordinates
# doubles that single precision cannot hold, so "kept exactly" means it.
open_sphere <- function(sub_division = 3L) {
  mesh <- vcg_sphere(sub_division = sub_division)
  cz <- colMeans(matrix(mesh$vb[3, mesh$it], nrow = 3L))
  mesh$it <- mesh$it[, cz < 0.5, drop = FALSE]
  mesh$vb[1:3, ] <- mesh$vb[1:3, ] * 1.1
  mesh$normals <- NULL
  mesh
}

# Vertices on an edge that belongs to exactly one face
border_vertices <- function(it) {
  e <- cbind(it[1:2, ], it[2:3, ], it[c(3, 1), ])
  key <- paste(pmin(e[1, ], e[2, ]), pmax(e[1, ], e[2, ]))
  once <- names(which(table(key) == 1L))
  sort(unique(as.integer(unlist(strsplit(once, " ", fixed = TRUE)))))
}

cross3 <- function(u, w) {
  c(u[2] * w[3] - u[3] * w[2], u[3] * w[1] - u[1] * w[3], u[1] * w[2] - u[2] * w[1])
}

# Face-edge Laplacian as a dense matrix
dense_laplacian <- function(vb, it, laplacian_weight = 1, cot = FALSE) {
  n <- ncol(vb)
  L <- matrix(0, n, n)
  for (f in seq_len(ncol(it))) {
    idx <- it[, f]
    for (e in 1:3) {
      a <- idx[e]
      b <- idx[e %% 3L + 1L]
      if (cot) {
        o  <- idx[(e + 1L) %% 3L + 1L]
        ca <- vb[, a] - vb[, o]
        cb <- vb[, b] - vb[, o]
        wt <- sum(ca * cb) / sqrt(sum(cross3(ca, cb)^2)) / 2
      } else {
        wt <- laplacian_weight
      }
      L[a, a] <- L[a, a] + wt
      L[b, b] <- L[b, b] + wt
      L[a, b] <- L[a, b] - wt
      L[b, a] <- L[b, a] - wt
    }
  }
  L
}

dense_implicit_reference <- function(mesh, lambda = 0.2, degree = 1L,
                                     laplacian_weight = 1, cot = FALSE) {
  vb <- mesh$vb[1:3, , drop = FALSE]
  it <- mesh$it
  h  <- numeric(ncol(vb))
  for (f in seq_len(ncol(it))) {
    idx <- it[, f]
    da  <- sqrt(sum(cross3(vb[, idx[2]] - vb[, idx[1]], vb[, idx[3]] - vb[, idx[1]])^2))
    h[idx] <- h[idx] + da
  }
  L <- dense_laplacian(vb, it, laplacian_weight = laplacian_weight, cot = cot)
  for (i in seq_len(degree - 1L)) L <- L %*% L
  M <- diag(h / max(h))
  t(solve(M + lambda * L, M %*% t(vb)))
}

# An ellipsoid, so that face areas and therefore the mass matrix vary
ellipsoid <- function(sub_division = 2L) {
  mesh <- vcg_sphere(sub_division = sub_division)
  mesh$vb[1, ] <- mesh$vb[1, ] * 1.5
  mesh$vb[2, ] <- mesh$vb[2, ] * 0.8
  mesh$normals <- NULL
  mesh
}

test_that("vcg_smooth_implicit solves the implicit system it documents", {
  mesh <- ellipsoid()
  for (degree in 1:3) {
    expect_equal(
      vcg_smooth_implicit(mesh, lambda = 0.2, degree = degree)$vb[1:3, ],
      dense_implicit_reference(mesh, lambda = 0.2, degree = degree),
      tolerance = 1e-5
    )
  }
  expect_equal(
    vcg_smooth_implicit(mesh, lambda = 0.5, laplacian_weight = 2)$vb[1:3, ],
    dense_implicit_reference(mesh, lambda = 0.5, laplacian_weight = 2),
    tolerance = 1e-5
  )
  expect_equal(
    vcg_smooth_implicit(mesh, lambda = 0.2, use_cot_weight = TRUE)$vb[1:3, ],
    dense_implicit_reference(mesh, lambda = 0.2, cot = TRUE),
    tolerance = 1e-5
  )
})

test_that("vcg_smooth_implicit holds border vertices exactly when fix_border = TRUE", {
  mesh     <- open_sphere()
  bnd      <- border_vertices(mesh$it)
  interior <- setdiff(sort(unique(as.vector(mesh$it))), bnd)

  for (degree in 1:3) {
    fixed <- vcg_smooth_implicit(mesh, degree = degree, fix_border = TRUE)
    expect_identical(fixed$vb[1:3, bnd], mesh$vb[1:3, bnd])
    # the interior is still smoothed
    expect_gt(max(abs(fixed$vb[1:3, interior] - mesh$vb[1:3, interior])), 1e-3)
  }

  # without the flag the rim is free to move
  free <- vcg_smooth_implicit(mesh, degree = 2, fix_border = FALSE)
  expect_gt(max(abs(free$vb[1:3, bnd] - mesh$vb[1:3, bnd])), 1e-3)
})

test_that("vcg_smooth_implicit without the mass matrix fairs the interior from the fixed border", {
  mesh     <- open_sphere()
  bnd      <- border_vertices(mesh$it)
  interior <- setdiff(sort(unique(as.vector(mesh$it))), bnd)

  faired <- vcg_smooth_implicit(mesh, use_mass_matrix = FALSE,
                                fix_border = TRUE, degree = 1)
  vb <- faired$vb[1:3, ]

  expect_true(all(is.finite(vb)))
  expect_identical(vb[, bnd], mesh$vb[1:3, bnd])

  # degree-1 fairing is harmonic: L x vanishes on every free vertex
  L <- dense_laplacian(mesh$vb[1:3, ], mesh$it)
  residual <- vb %*% L
  expect_lt(max(abs(residual[, interior])), 1e-6)

  # with no data term, lambda only scales the system
  expect_equal(
    vcg_smooth_implicit(mesh, lambda = 5, use_mass_matrix = FALSE,
                        fix_border = TRUE, degree = 1)$vb,
    faired$vb, tolerance = 1e-8
  )
})

test_that("vcg_smooth_implicit refuses to fair a mesh with nothing held fixed", {
  # a closed mesh has no border to hold
  expect_error(
    vcg_smooth_implicit(vcg_sphere(sub_division = 2L),
                        use_mass_matrix = FALSE, fix_border = TRUE),
    "use_mass_matrix"
  )
  # an open mesh whose border is not held
  expect_error(
    vcg_smooth_implicit(open_sphere(), use_mass_matrix = FALSE,
                        fix_border = FALSE),
    "use_mass_matrix"
  )
})

test_that("vcg_smooth_implicit stops before a solve that needs more than max_memory", {
  expect_error(
    vcg_smooth_implicit(vcg_sphere(sub_division = 2L), max_memory = 1e-6),
    "GiB"
  )
})

test_that("vcg_smooth_implicit reports the whole memory need, not the first stage over the limit", {
  sphere <- vcg_sphere(sub_division = 3L)
  reported <- function(max_memory) {
    msg <- tryCatch({
      vcg_smooth_implicit(sphere, degree = 2, max_memory = max_memory)
      ""
    }, error = conditionMessage)
    as.numeric(sub(".*needs (about|at least) ([0-9.eE+-]+) GiB.*", "\\2", msg))
  }
  # a limit below even the first stage
  floor_need <- reported(1e-9)
  expect_true(is.finite(floor_need))

  # past the first stage the error names the total, wherever the limit falls
  total <- reported(floor_need * 1.01)
  expect_gt(total, floor_need)
  expect_equal(reported(total * 0.99), total, tolerance = 1e-2)

  # and a limit just above the total is enough
  smoothed <- vcg_smooth_implicit(sphere, degree = 2, max_memory = total * 1.01)
  expect_true(all(is.finite(smoothed$vb)))
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

test_that("vcg_mesh_volume reports a stable non-manifold edge count", {
  # vcglib's CountNonManifoldEdgeFF used to leak three process-lifetime user
  # bit flags per call; once the allocator ran past the sign bit it handed out
  # 0, the de-duplication stopped working, and the count drifted 3 -> 4 -> 9.
  broken <- vcg_sphere()
  broken$it <- cbind(broken$it, broken$it[, 1, drop = FALSE])

  counts <- vapply(seq_len(8L), function(i) {
    msg <- gsub("\n", " ", tryCatch(vcg_mesh_volume(broken), error = conditionMessage))
    as.integer(sub(".*Non-manifold edges:\\s*(\\d+).*", "\\1", msg))
  }, integer(1L))

  # the duplicated face makes exactly its three edges non-manifold
  expect_equal(counts, rep(3L, 8L))
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
