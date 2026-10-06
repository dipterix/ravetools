
# Quadric edge-collapse decimation (src/vcgDecimate.cpp). The expectations are
# geometric facts checked independently of vcglib: face counts, distance to
# the unit sphere, edge multiplicities and the Euler characteristic.

# Edge multiplicities of a face matrix, keyed "a b" with a < b
edge_counts <- function(it) {
  e <- cbind(it[1:2, ], it[2:3, ], it[c(3, 1), ])
  table(paste(pmin(e[1, ], e[2, ]), pmax(e[1, ], e[2, ])))
}

# V - E + F over the vertices that faces reference
euler_characteristic <- function(it) {
  length(unique(as.vector(it))) - length(edge_counts(it)) + ncol(it)
}

test_that("vcg_decimate reaches the requested face count", {
  sphere <- vcg_sphere(sub_division = 4L)            # 5120 faces
  target <- ncol(sphere$it) / 4

  decimated <- vcg_decimate(sphere, ratio = 0.25)
  expect_lte(ncol(decimated$it), target)
  expect_gte(ncol(decimated$it), target - 2)

  by_count <- vcg_decimate(sphere, target_faces = 1000)
  expect_lte(ncol(by_count$it), 1000)
  expect_gte(ncol(by_count$it), 998)
})

test_that("vcg_decimate keeps the vertices on the surface it simplifies", {
  sphere    <- vcg_sphere(sub_division = 4L)
  decimated <- vcg_decimate(sphere, ratio = 0.25)

  radii <- sqrt(colSums(decimated$vb[1:3, ]^2))
  expect_true(all(radii > 0.98 & radii < 1.02))
  expect_true(all(is.finite(decimated$normals)))
  expect_equal(ncol(decimated$normals), ncol(decimated$vb))
})

test_that("vcg_decimate keeps a closed manifold mesh closed and manifold", {
  decimated <- vcg_decimate(vcg_sphere(sub_division = 4L), ratio = 0.1)
  counts <- edge_counts(decimated$it)
  expect_true(all(counts == 2L))
})

test_that("vcg_decimate preserves topology when asked", {
  # two disjoint spheres: Euler characteristic 2 + 2
  one <- vcg_sphere(sub_division = 3L)
  two <- one
  two$vb[1, ] <- two$vb[1, ] + 5
  pair <- one
  pair$vb <- cbind(one$vb, two$vb)
  pair$it <- cbind(one$it, two$it + ncol(one$vb))
  pair$normals <- NULL
  expect_equal(euler_characteristic(pair$it), 4L)

  decimated <- vcg_decimate(pair, ratio = 0.05, preserve_topology = TRUE)
  expect_equal(euler_characteristic(decimated$it), 4L)
  # both spheres are still there
  expect_true(any(decimated$vb[1, ] > 4) && any(decimated$vb[1, ] < 1))
})

test_that("vcg_decimate leaves the boundary where it was", {
  sphere <- vcg_sphere(sub_division = 4L)
  cz <- colMeans(matrix(sphere$vb[3, sphere$it], nrow = 3L))
  cap <- sphere
  cap$it <- sphere$it[, cz < 0.5, drop = FALSE]
  cap$normals <- NULL

  border_positions <- function(mesh) {
    counts <- edge_counts(mesh$it)
    ids <- unique(as.integer(unlist(strsplit(names(counts)[counts == 1L], " "))))
    mesh$vb[1:3, ids, drop = FALSE]
  }
  before <- border_positions(cap)
  after  <- border_positions(vcg_decimate(cap, ratio = 0.3, preserve_boundary = TRUE))

  # every boundary vertex left is one of the original boundary vertices
  gaps <- apply(after, 2, function(p) min(colSums((before - p)^2)))
  expect_true(all(gaps < 1e-12))
  expect_equal(ncol(after), ncol(before))
})

test_that("vcg_decimate validates its targets", {
  sphere <- vcg_sphere(sub_division = 2L)
  expect_error(vcg_decimate(sphere, ratio = 0), "ratio")
  expect_error(vcg_decimate(sphere, ratio = 1.5), "ratio")
  expect_error(vcg_decimate(sphere$vb[1:3, ], ratio = 0.5))

  # nothing to remove: the mesh comes back as it was
  same <- vcg_decimate(sphere, target_faces = ncol(sphere$it))
  expect_equal(same$vb[1:3, ], sphere$vb[1:3, ])
  expect_equal(same$it, sphere$it)
})
