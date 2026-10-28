
# Convex hulls (src/vcgConvexHull.cpp, R/vcg-convex-hull.R). The expectations
# are geometric facts checked independently of either backend: containment of
# every input point, a closed and consistently oriented surface, and volumes
# from the divergence theorem.

hull_methods <- c("vcg", "quickhull3d")

# Edge multiplicities of a face matrix, keyed "a b" with a < b
edge_counts <- function(it) {
  e <- cbind(it[1:2, ], it[2:3, ], it[c(3, 1), ])
  table(paste(pmin(e[1, ], e[2, ]), pmax(e[1, ], e[2, ])))
}

# Closed (every edge in two faces), consistently oriented (each directed edge
# once) and genus 0 (V - E + F = 2)
expect_closed_surface <- function(hull) {
  counts <- edge_counts(hull$it)
  expect_true(all(counts == 2L))
  e <- cbind(hull$it[1:2, ], hull$it[2:3, ], hull$it[c(3, 1), ])
  expect_equal(anyDuplicated(paste(e[1, ], e[2, ])), 0L)
  n_vertex <- length(unique(as.vector(hull$it)))
  expect_equal(n_vertex - length(counts) + ncol(hull$it), 2)
}

# Largest signed distance from a point (rows of `points`) to the plane of a
# hull face, positive outside; faces too thin to have a plane are skipped
max_outside_distance <- function(hull, points) {
  v <- hull$vb[1:3, , drop = FALSE]
  a <- v[, hull$it[1, ], drop = FALSE]
  u <- v[, hull$it[2, ], drop = FALSE] - a
  w <- v[, hull$it[3, ], drop = FALSE] - a
  n <- rbind(u[2, ] * w[3, ] - u[3, ] * w[2, ],
             u[3, ] * w[1, ] - u[1, ] * w[3, ],
             u[1, ] * w[2, ] - u[2, ] * w[1, ])
  len <- sqrt(colSums(n^2))
  keep <- len > 1e-12 * max(len)
  n <- n[, keep, drop = FALSE] / rep(len[keep], each = 3L)
  offset <- colSums(n * a[, keep, drop = FALSE])
  max(crossprod(n, t(points)) - offset)
}

# Volume enclosed by a closed mesh; positive when the faces are ordered
# counter-clockwise seen from outside
signed_volume <- function(hull) {
  v <- hull$vb[1:3, , drop = FALSE]
  a <- v[, hull$it[1, ], drop = FALSE]
  b <- v[, hull$it[2, ], drop = FALSE]
  c <- v[, hull$it[3, ], drop = FALSE]
  sum(a[1, ] * (b[2, ] * c[3, ] - b[3, ] * c[2, ]) +
      a[2, ] * (b[3, ] * c[1, ] - b[1, ] * c[3, ]) +
      a[3, ] * (b[1, ] * c[2, ] - b[2, ] * c[1, ])) / 6
}

test_that("vcg_convex_hull keeps exactly the corners of a filled cube", {
  set.seed(1)
  x <- rbind(
    as.matrix(expand.grid(0:1, 0:1, 0:1)),
    matrix(runif(600, 0.05, 0.95), ncol = 3)
  )
  for (method in hull_methods) {
    hull <- vcg_convex_hull(x, method = method)
    expect_equal(attr(hull, "vertex_index"), 1:8)
    expect_equal(signed_volume(hull), 1)
    expect_closed_surface(hull)
  }
})

test_that("vcg_convex_hull returns a closed outward surface around every point", {
  set.seed(2)
  x <- matrix(rnorm(3000), ncol = 3)
  for (method in hull_methods) {
    hull <- vcg_convex_hull(x, method = method)
    index <- attr(hull, "vertex_index")
    expect_false(is.unsorted(index, strictly = TRUE))
    # hull vertices are the input points themselves, not rounded copies
    expect_identical(hull$vb[1:3, ], t(x[index, ]))
    expect_lt(max_outside_distance(hull, x), 1e-9)
    expect_closed_surface(hull)
    expect_gt(signed_volume(hull), 0)
    expect_no_warning(volume <- vcg_mesh_volume(hull))
    expect_equal(volume, signed_volume(hull), tolerance = 1e-5)
    expect_equal(dim(hull$normals), dim(hull$vb))
  }
})

test_that("both methods find the same hull for points in general position", {
  set.seed(2)
  x <- matrix(rnorm(3000), ncol = 3)
  by_vcg <- vcg_convex_hull(x, method = "vcg")
  by_quickhull3d <- vcg_convex_hull(x, method = "quickhull3d")
  expect_identical(attr(by_quickhull3d, "vertex_index"), attr(by_vcg, "vertex_index"))
  expect_equal(signed_volume(by_quickhull3d), signed_volume(by_vcg))
})

test_that("the quickhull3d port reproduces quickhull3d 3.1.2", {
  # faces that quickhull3d 3.1.2 (Node.js) returns for the same points, as
  # 1-based row indices; the tetrahedron is the example in its README
  tet <- rbind(c(0, 1, 0), c(1, -1, 1), c(-1, -1, 1), c(0, -1, -1))
  expect_identical(
    quickhull3dConvexHull(t(tet), -1),
    matrix(c(3L, 2L, 1L, 4L, 2L, 3L, 4L, 1L, 2L, 4L, 3L, 1L), nrow = 3L)
  )

  set.seed(5)
  random20 <- matrix(rnorm(60), ncol = 3)
  expect_identical(
    quickhull3dConvexHull(t(random20), -1),
    matrix(c(11L, 5L, 13L, 4L, 11L, 13L, 19L, 5L, 11L, 1L, 18L, 4L, 1L,
             4L, 13L, 1L, 13L, 18L, 3L, 13L, 7L, 3L, 7L, 18L, 3L, 18L, 13L,
             17L, 4L, 18L, 17L, 18L, 19L, 17L, 19L, 11L, 17L, 11L, 4L, 9L,
             19L, 18L, 2L, 5L, 19L, 2L, 19L, 9L, 2L, 7L, 13L, 2L, 13L, 5L,
             16L, 7L, 2L, 16L, 2L, 9L, 16L, 9L, 18L, 16L, 18L, 7L),
           nrow = 3L)
  )

  # coplanar faces merged into squares, then split into fans; exact ties make
  # this one sensitive to how the compiler rounds
  skip_on_cran()
  grid <- as.matrix(expand.grid(0:2, 0:2, 0:2))
  storage.mode(grid) <- "double"
  expect_identical(
    quickhull3dConvexHull(t(grid), -1),
    matrix(c(19L, 25L, 7L, 19L, 7L, 1L, 9L, 7L, 25L, 9L, 25L, 27L, 9L, 3L,
             1L, 9L, 1L, 7L, 21L, 27L, 25L, 21L, 25L, 19L, 21L, 19L, 1L, 21L,
             1L, 3L, 21L, 3L, 9L, 21L, 9L, 27L),
           nrow = 3L)
  )
})

test_that("the quickhull3d port rejects degenerate points instead of crashing", {
  # vcg_convex_hull screens these out first; the port must not rely on it
  expect_error(quickhull3dConvexHull(diag(3), -1), "<4 points")
  expect_error(quickhull3dConvexHull(matrix(1, 3, 5), -1), "identical")
  expect_error(quickhull3dConvexHull(rbind(1:10, 2 * (1:10), 3 * (1:10)) + 0, -1),
               "collinear")
  set.seed(6)
  expect_error(quickhull3dConvexHull(rbind(runif(10), runif(10), 0), -1),
               "coplanar")
})

test_that("vcg_convex_hull keeps every vertex of a sphere and ignores faces", {
  sphere <- vcg_sphere(sub_division = 3L)
  for (method in hull_methods) {
    hull <- vcg_convex_hull(sphere, method = method)
    expect_equal(attr(hull, "vertex_index"), seq_len(ncol(sphere$vb)))
    from_points <- vcg_convex_hull(t(sphere$vb[1:3, ]), method = method)
    expect_identical(from_points$it, hull$it)
  }
})

test_that("vcg_convex_hull treats points within `tolerance` of a face as on it", {
  # a unit corner tetrahedron and a point 0.001 outside its slanted face
  x <- rbind(c(0, 0, 0), diag(3), rep(1 / 3, 3) + 1e-3 / sqrt(3))
  for (method in hull_methods) {
    expect_equal(attr(vcg_convex_hull(x, method = method), "vertex_index"), 1:5)
    expect_equal(
      attr(vcg_convex_hull(x, method = method, tolerance = 0.01), "vertex_index"),
      1:4
    )
  }
})

test_that("vcg_convex_hull handles grids full of coplanar points", {
  grid <- as.matrix(expand.grid(0:5, 0:5, 0:5))
  corners <- which(rowSums(grid == 0 | grid == 5) == 3)
  for (method in hull_methods) {
    hull <- vcg_convex_hull(grid, method = method)
    expect_closed_surface(hull)
    expect_lt(max_outside_distance(hull, grid), 1e-9)
    expect_equal(signed_volume(hull), 125)
    expect_true(all(corners %in% attr(hull, "vertex_index")))
  }

  mesh <- vcg_isosurface(left_hippocampus_mask)
  points <- t(mesh$vb[1:3, ])
  for (method in hull_methods) {
    hull <- vcg_convex_hull(mesh, method = method)
    expect_closed_surface(hull)
    expect_lt(max_outside_distance(hull, points), 1e-6)
  }
})

test_that("vcg_convex_hull raises a tolerance below rounding error", {
  # A rotated grid: many points share planes, but only up to rounding. Below
  # the rounding error the visibility tests contradict each other; the vcg
  # backend then lost face links and crashed R, and the quickhull3d port (like
  # quickhull3d itself) merged faces forever on this very grid. Such
  # tolerances are raised to the rounding error, 3 times the machine
  # precision times the sum over the axes of the largest absolute coordinate.
  set.seed(704)
  rotation <- qr.Q(qr(matrix(rnorm(9), 3)))
  points <- as.matrix(expand.grid(0:7, 0:7, 0:7)) %*% rotation
  rounding <- 3 * .Machine$double.eps * sum(apply(abs(points), 2, max))
  for (method in hull_methods) {
    at_rounding <- tryCatch(
      vcg_convex_hull(points, method = method, tolerance = rounding),
      error = function(e) e)
    for (tolerance in c(0, rounding / 100)) {
      below <- tryCatch(
        vcg_convex_hull(points, method = method, tolerance = tolerance),
        error = function(e) e)
      if (inherits(at_rounding, "error")) {
        # a numerical failure at the rounding error is reported, not returned
        expect_s3_class(below, "error")
        expect_match(conditionMessage(below), "larger `tolerance`")
      } else {
        expect_closed_surface(below)
        expect_identical(below$it, at_rounding$it)
      }
    }
  }
})

test_that("vcg_convex_hull does not depend on the position or scale of the points", {
  set.seed(3)
  x <- matrix(rnorm(600), ncol = 3)
  for (method in hull_methods) {
    index <- attr(vcg_convex_hull(x, method = method), "vertex_index")
    expect_identical(
      attr(vcg_convex_hull(x * 1e-6, method = method), "vertex_index"), index)
    expect_identical(
      attr(vcg_convex_hull(x * 1000 + 5000, method = method), "vertex_index"), index)

    # a repeated point is kept once, through one of its copies
    doubled <- vcg_convex_hull(rbind(x, x), method = method)
    expect_closed_surface(doubled)
    expect_setequal((attr(doubled, "vertex_index") - 1L) %% nrow(x) + 1L, index)
  }
})

test_that("is_closed_surface accepts only closed, oriented genus-0 meshes", {
  # a consistently oriented tetrahedron on vertices 1 to 4
  tet <- cbind(c(1, 3, 2), c(1, 2, 4), c(2, 3, 4), c(3, 1, 4))
  expect_true(is_closed_surface(tet, 4))
  # a hole
  expect_false(is_closed_surface(tet[, -4], 4))
  # one face reversed
  flipped <- tet
  flipped[, 2] <- flipped[3:1, 2]
  expect_false(is_closed_surface(flipped, 4))
  # two tetrahedra pinched at vertex 4: every edge is fine, V - E + F = 3
  second <- matrix(c(5, 6, 7, 4)[tet], nrow = 3)
  expect_false(is_closed_surface(cbind(tet, second), 7))
})

test_that("vcg_convex_hull explains why it cannot build a hull", {
  set.seed(4)
  points <- matrix(rnorm(30), ncol = 3)
  for (method in hull_methods) {
    expect_error(vcg_convex_hull(diag(3), method = method), "at least 4 points")
    expect_error(vcg_convex_hull(cbind(runif(10), runif(10), 0), method = method),
                 "enclose no volume")
    expect_error(vcg_convex_hull(cbind(1:10, 2 * (1:10), 3 * (1:10)), method = method),
                 "enclose no volume")
    expect_error(vcg_convex_hull(matrix(1, 5, 3), method = method),
                 "enclose no volume")
    with_na <- points
    with_na[2, 1] <- NA
    expect_error(vcg_convex_hull(with_na, method = method), "finite")
    expect_error(vcg_convex_hull(points, method = method, tolerance = -1),
                 "tolerance")
  }
  expect_error(vcg_convex_hull(points, method = "qhull"))
})
