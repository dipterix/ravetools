test_that("vcg_detect_collision matches an independent oracle", {

  # ---- Independent oracle -------------------------------------------------
  # Densely sample each geometry and take the minimum pairwise distance. This
  # shares no algorithm with the implementation, so it catches the failure mode
  # that motivates the hand-written primitives: vcglib's SegmentSegmentDistance
  # clamps the two line-line closest points independently and overestimates by
  # up to ~4x on skew pairs.
  #
  # Sampling can only overestimate, and by at most half a sample spacing per
  # side, which gives a two-sided bound on the implementation.

  sample_point <- function(p) list(pts = matrix(p, nrow = 1L), h = 0)

  sample_segment <- function(a, b, n = 400L) {
    t <- seq(0, 1, length.out = n)
    list(
      pts = cbind(a[1] + t * (b[1] - a[1]),
                  a[2] + t * (b[2] - a[2]),
                  a[3] + t * (b[3] - a[3])),
      h = sqrt(sum((b - a)^2)) / (n - 1L)
    )
  }

  sample_triangle <- function(v0, v1, v2, k = 40L) {
    grid <- expand.grid(i = 0:k, j = 0:k)
    grid <- grid[grid$i + grid$j <= k, , drop = FALSE]
    u <- grid$i / k
    v <- grid$j / k
    w <- 1 - u - v
    edge <- max(sqrt(sum((v1 - v0)^2)), sqrt(sum((v2 - v1)^2)), sqrt(sum((v0 - v2)^2)))
    list(
      pts = cbind(w * v0[1] + u * v1[1] + v * v2[1],
                  w * v0[2] + u * v1[2] + v * v2[2],
                  w * v0[3] + u * v1[3] + v * v2[3]),
      h = edge / k
    )
  }

  min_pair_dist <- function(A, B) {
    best <- Inf
    for (i in seq_len(nrow(A))) {
      best <- min(best, min(sqrt(colSums((t(B) - A[i, ])^2))))
    }
    best
  }

  # `impl` is the implementation's answer, `sa`/`sb` the sampled geometries.
  expect_matches_oracle <- function(impl, sa, sb, label) {
    sampled <- min_pair_dist(sa$pts, sb$pts)
    tol <- (sa$h + sb$h) / 2 + 1e-5
    # Must never overestimate (the vcglib bug class) ...
    expect_lt(impl, sampled + 1e-5)
    # ... and must not underestimate beyond the sampling resolution.
    expect_lt(sampled - impl, tol)
    invisible(label)
  }

  measure <- function(x, y, mode_x, mode_y) {
    # `radius = Inf` is not allowed, so use a bound larger than any distance here
    res <- vcg_detect_collision(x, y, mode_x = mode_x, mode_y = mode_y,
                                radius = 1000)
    res$distance[[1L]]
  }

  set.seed(20260819)

  # -- point vs point -------------------------------------------------------
  for (i in 1:20) {
    p <- matrix(runif(6, -3, 3), 2L, 3L)
    d <- measure(p[1, , drop = FALSE], p[2, , drop = FALSE], "points", "points")
    expect_equal(d, sqrt(sum((p[1, ] - p[2, ])^2)), tolerance = 1e-5)
  }

  # -- point vs segment -----------------------------------------------------
  for (i in 1:20) {
    p <- matrix(runif(9, -3, 3), 3L, 3L)
    d <- measure(p[1, , drop = FALSE], p[2:3, , drop = FALSE], "points", "segments")
    expect_matches_oracle(d, sample_point(p[1, ]), sample_segment(p[2, ], p[3, ]),
                          "point-segment")
  }

  # -- segment vs segment ---------------------------------------------------
  for (i in 1:40) {
    p <- matrix(runif(12, -3, 3), 4L, 3L)
    d <- measure(p[1:2, , drop = FALSE], p[3:4, , drop = FALSE], "segments", "segments")
    expect_matches_oracle(d, sample_segment(p[1, ], p[2, ]),
                          sample_segment(p[3, ], p[4, ]), "segment-segment")
  }

  # -- point vs triangle ----------------------------------------------------
  one_triangle <- function(v) {
    list(vb = rbind(t(v), 1), it = matrix(1:3, nrow = 3L),
         class = "mesh3d")
  }
  as_mesh <- function(v) {
    structure(list(vb = rbind(t(v), 1), it = matrix(1:3, ncol = 1L)),
              class = c("ravetools_mesh3d", "mesh3d"))
  }

  for (i in 1:20) {
    v <- matrix(runif(9, -3, 3), 3L, 3L)
    p <- matrix(runif(3, -3, 3), 1L, 3L)
    d <- measure(p, as_mesh(v), "points", "mesh")
    expect_matches_oracle(d, sample_point(p[1, ]),
                          sample_triangle(v[1, ], v[2, ], v[3, ]), "point-triangle")
  }

  # -- segment vs triangle --------------------------------------------------
  for (i in 1:20) {
    v <- matrix(runif(9, -3, 3), 3L, 3L)
    s <- matrix(runif(6, -3, 3), 2L, 3L)
    d <- measure(s, as_mesh(v), "segments", "mesh")
    expect_matches_oracle(d, sample_segment(s[1, ], s[2, ]),
                          sample_triangle(v[1, ], v[2, ], v[3, ]), "segment-triangle")
  }

  # -- triangle vs triangle -------------------------------------------------
  for (i in 1:20) {
    a <- matrix(runif(9, -3, 3), 3L, 3L)
    b <- matrix(runif(9, -3, 3), 3L, 3L)
    d <- measure(as_mesh(a), as_mesh(b), "mesh", "mesh")
    expect_matches_oracle(d, sample_triangle(a[1, ], a[2, ], a[3, ]),
                          sample_triangle(b[1, ], b[2, ], b[3, ]), "triangle-triangle")
  }
})


test_that("vcg_detect_collision is exact where vcglib's own routine is not", {
  # vcg::SegmentSegmentDistance returns 4.123 for this pair; the truth is 3.0.
  x <- rbind(c(0, 0, 0), c(1, 0, 0))
  y <- rbind(c(5, 0, 1), c(1, 4, 1))
  res <- vcg_detect_collision(x, y, mode_x = "segments", mode_y = "segments",
                              radius = 10)
  expect_equal(res$distance[[1L]], 3, tolerance = 1e-5)

  # The case the tractography viewer gets wrong: a segment passing straight
  # through a single-point ROI, with both of its vertices offset. Vertex
  # sampling reports 0.781; the true segment distance is 0.
  roi <- rbind(c(0, 0, 0))
  seg <- rbind(c(-0.781, 0, 0), c(0.781, 0, 0))
  res <- vcg_detect_collision(roi, seg, mode_y = "segments", radius = 1)
  expect_true(res$hit[[1L]])
  expect_equal(res$distance[[1L]], 0, tolerance = 1e-6)
})


test_that("vcg_detect_collision agrees with vcg_kdtree_nearest on point clouds", {
  set.seed(7)
  target <- matrix(rnorm(300), ncol = 3L)
  query <- matrix(rnorm(600), ncol = 3L)

  for (radius in c(0.25, 0.5, 1, 2)) {
    res <- vcg_detect_collision(target, query, mode_x = "points",
                                mode_y = "points", radius = radius)
    knn <- vcg_kdtree_nearest(target = target, query = query, k = 1L)

    expect_identical(as.logical(res$hit), as.vector(knn$distance <= radius))

    # Where there is a hit, both the distance and the reported target row must
    # match the nearest neighbour.
    hit <- which(res$hit)
    expect_equal(res$distance[hit], as.vector(knn$distance)[hit], tolerance = 1e-5)
    expect_identical(res$index[hit], as.integer(knn$index)[hit])
  }
})


test_that("vcg_detect_collision is symmetric in its two sides", {
  set.seed(11)
  a <- matrix(rnorm(60), ncol = 3L)
  b <- matrix(rnorm(90), ncol = 3L)

  for (modes in list(c("points", "points"), c("points", "segments"),
                     c("segments", "segments"))) {
    forward <- vcg_detect_collision(a, b, mode_x = modes[[1L]],
                                    mode_y = modes[[2L]], radius = 1.5)
    reverse <- vcg_detect_collision(b, a, mode_x = modes[[2L]],
                                    mode_y = modes[[1L]], radius = 1.5)

    expect_identical(any(forward$hit, na.rm = TRUE),
                     any(reverse$hit, na.rm = TRUE))
    expect_equal(min(forward$distance, na.rm = TRUE),
                 min(reverse$distance, na.rm = TRUE), tolerance = 1e-5)
  }
})


test_that("vcg_detect_collision accepts NA rows as polyline separators", {
  set.seed(3)
  roi <- matrix(rnorm(30), ncol = 3L)

  tract1 <- cbind(seq(-2, 2, by = 0.4), 0, 0)
  tract2 <- cbind(0, seq(-2, 2, by = 0.4), 1)
  combined <- rbind(tract1, c(NA, NA, NA), tract2)

  joint <- vcg_detect_collision(roi, combined, mode_y = "segments", radius = 1)

  sep <- nrow(tract1) + 1L
  expect_true(is.na(joint$hit[[sep]]))            # the separator row
  expect_true(is.na(joint$hit[[nrow(tract1)]]))   # last vertex of tract 1
  expect_true(is.na(joint$hit[[nrow(combined)]])) # last vertex of tract 2

  # Splitting the matrix into one call per polyline must give the same answers.
  solo1 <- vcg_detect_collision(roi, tract1, mode_y = "segments", radius = 1)
  solo2 <- vcg_detect_collision(roi, tract2, mode_y = "segments", radius = 1)

  expect_identical(joint$hit[seq_len(nrow(tract1))], solo1$hit)
  expect_equal(joint$distance[seq_len(nrow(tract1))], solo1$distance)
  expect_identical(joint$hit[(sep + 1L):nrow(combined)], solo2$hit)
  expect_equal(joint$distance[(sep + 1L):nrow(combined)], solo2$distance)

  # A separator on the indexed side must not create a phantom segment that
  # bridges the end of one polyline to the start of the next. This probe sits
  # near the midpoint of that would-be bridge and far from either polyline.
  a <- rbind(c(0, 0, 0), c(1, 0, 0))     # polyline 1
  b <- rbind(c(1, 20, 0), c(0, 20, 0))   # polyline 2, well away
  probe <- rbind(c(1, 10, 0))            # midway along the phantom bridge

  separated <- vcg_detect_collision(rbind(a, c(NA, NA, NA), b), probe,
                                    mode_x = "segments", radius = 1)
  bridged <- vcg_detect_collision(rbind(a, b), probe,
                                  mode_x = "segments", radius = 1)
  expect_false(separated$hit)   # nearest real segment is 10 away
  expect_true(bridged$hit)      # without the separator, the bridge passes through

  # Rows that are only partly finite are treated as separators too.
  partial <- combined
  partial[3L, 2L] <- NA_real_
  expect_true(is.na(vcg_detect_collision(roi, partial, mode_y = "segments",
                                         radius = 1)$hit[[3L]]))
})


test_that("vcg_detect_collision early_stop agrees with a full scan", {
  set.seed(5)
  roi <- matrix(rnorm(30), ncol = 3L)
  tracts <- rbind(
    cbind(seq(-2, 2, by = 0.25), 0, 0),
    c(NA, NA, NA),
    cbind(0.2, seq(-2, 2, by = 0.25), 0)
  )

  full <- vcg_detect_collision(roi, tracts, mode_y = "segments", radius = 0.8)
  early <- vcg_detect_collision(roi, tracts, mode_y = "segments", radius = 0.8,
                                early_stop = TRUE)

  # Every element the early-stopping scan did evaluate must agree.
  evaluated <- !is.na(early$hit)
  expect_identical(early$hit[evaluated], full$hit[evaluated])
  expect_equal(early$distance[evaluated], full$distance[evaluated])

  # Once it stops, it stops: nothing after the first hit in a group is
  # evaluated, and at most one hit is reported per group.
  expect_lte(sum(early$hit, na.rm = TRUE), 2L)
  expect_true(any(full$hit, na.rm = TRUE))

  # The unit is the NA-delimited group, not the row. A point cloud carrying no
  # separators is therefore ONE group, and stops after a single hit overall
  # rather than reporting every near point.
  origin <- rbind(c(0, 0, 0))
  cloud <- rbind(c(0, 0, 0), c(0.1, 0, 0), c(9, 9, 9), c(0.2, 0, 0))

  flat <- vcg_detect_collision(origin, cloud, mode_y = "points", radius = 0.5,
                               early_stop = TRUE)
  expect_identical(flat$hit, c(TRUE, NA, NA, NA))

  # Splitting the same cloud with a separator yields one hit per group. The
  # mode is stated explicitly throughout because these assertions are about
  # point-cloud scoping specifically, not about what `auto` happens to guess.
  grouped <- vcg_detect_collision(origin,
                                  rbind(cloud[1:2, ], c(NA, NA, NA), cloud[3:4, ]),
                                  mode_y = "points", radius = 0.5,
                                  early_stop = TRUE)
  expect_identical(grouped$hit, c(TRUE, NA, NA, FALSE, TRUE))

  # A target mesh is likewise a single group: it stops at the first colliding
  # face, where a full scan finds many.
  a <- vcg_sphere()
  b <- vcg_sphere()
  b$vb[1L, ] <- b$vb[1L, ] + 1
  expect_gt(sum(vcg_detect_collision(a, b, radius = 0)$hit), 1L)
  expect_identical(sum(vcg_detect_collision(a, b, radius = 0,
                                            early_stop = TRUE)$hit,
                       na.rm = TRUE), 1L)

  # The hit reported is the first in row order, not the closest approach. This
  # chain grazes the origin at 0.5 before returning at 0.1.
  chain <- rbind(c(-2, 0.5, 0), c(2, 0.5, 0), c(2, 0.1, 0), c(-2, 0.1, 0))
  scan_all <- vcg_detect_collision(origin, chain, mode_y = "segments",
                                   radius = 0.6)
  scan_early <- vcg_detect_collision(origin, chain, mode_y = "segments",
                                     radius = 0.6, early_stop = TRUE)

  expect_equal(min(scan_all$distance, na.rm = TRUE), 0.1, tolerance = 1e-5)
  expect_equal(scan_early$distance[which(scan_early$hit)], 0.5, tolerance = 1e-5)
  expect_gt(scan_early$distance[which(scan_early$hit)],
            min(scan_all$distance, na.rm = TRUE))
})


test_that("vcg_detect_collision accepted set grows monotonically with radius", {
  set.seed(13)
  roi <- matrix(rnorm(60), ncol = 3L)
  query <- matrix(rnorm(300), ncol = 3L)

  previous <- integer(0)
  for (radius in c(0.1, 0.3, 0.6, 1, 1.5, 3)) {
    hit <- which(vcg_detect_collision(roi, query, radius = radius)$hit)
    expect_true(all(previous %in% hit))
    previous <- hit
  }
  expect_gt(length(previous), 0L)

  # A radius that spans the whole bounding box must accept everything.
  span <- vcg_detect_collision(roi, query, radius = 1e6)
  expect_true(all(span$hit))
})


test_that("vcg_detect_collision handles mesh-to-mesh proximity", {
  a <- vcg_sphere()
  b <- vcg_sphere()
  b$vb[1L, ] <- b$vb[1L, ] + 3

  # Two unit spheres 3 apart leave a surface gap of about 1.
  expect_false(any(vcg_detect_collision(a, b, radius = 0.9)$hit))
  expect_true(any(vcg_detect_collision(a, b, radius = 1.1)$hit))

  # Overlapping spheres genuinely intersect, so some faces are at distance 0.
  c <- vcg_sphere()
  c$vb[1L, ] <- c$vb[1L, ] + 1
  overlap <- vcg_detect_collision(a, c, radius = 0)
  expect_true(any(overlap$hit))
  expect_equal(min(overlap$distance, na.rm = TRUE), 0, tolerance = 1e-6)

  # The reported index must point at a real face of `a`.
  hit <- which(overlap$hit)
  expect_true(all(overlap$index[hit] >= 1L & overlap$index[hit] <= ncol(a$it)))
})


test_that("vcg_detect_collision can count the interior of a closed surface", {
  roi <- vcg_sphere()
  centre <- rbind(c(0, 0, 0))

  expect_false(vcg_detect_collision(roi, centre, radius = 0.1)$hit)

  inside <- vcg_detect_collision(roi, centre, radius = 0.1,
                                 include_interior = TRUE)
  expect_true(inside$hit)
  expect_equal(inside$distance, 0)

  # A point well outside stays a miss either way.
  outside <- rbind(c(10, 10, 10))
  expect_false(vcg_detect_collision(roi, outside, radius = 0.1,
                                    include_interior = TRUE)$hit)

  # A segment that lies wholly inside is only caught with the interior test.
  seg <- rbind(c(-0.2, 0, 0), c(0.2, 0, 0))
  expect_false(any(vcg_detect_collision(roi, seg, mode_y = "segments",
                                        radius = 0.1)$hit, na.rm = TRUE))
  expect_true(any(vcg_detect_collision(roi, seg, mode_y = "segments",
                                       radius = 0.1,
                                       include_interior = TRUE)$hit,
                  na.rm = TRUE))

  # It needs a mesh on the indexed side.
  expect_error(
    vcg_detect_collision(matrix(rnorm(30), ncol = 3L), centre, radius = 1,
                         include_interior = TRUE),
    "include_interior"
  )
})


test_that("vcg_detect_collision validates and resolves its arguments", {
  roi <- vcg_sphere()
  pts <- matrix(rnorm(30), ncol = 3L)

  # `auto` picks mesh for a surface with faces, points for any matrix. Separator
  # rows do not change that guess, so a matrix of connected segments must say so.
  expect_identical(
    vcg_detect_collision(roi, pts, radius = 1),
    vcg_detect_collision(roi, pts, mode_x = "mesh", mode_y = "points",
                         radius = 1)
  )
  chained <- rbind(pts[1:4, ], c(NA, NA, NA), pts[5:8, ])
  expect_identical(
    vcg_detect_collision(roi, chained, radius = 1),
    vcg_detect_collision(roi, chained, mode_x = "mesh", mode_y = "points",
                         radius = 1)
  )

  # A matrix cannot be a mesh.
  expect_error(vcg_detect_collision(pts, pts, mode_x = "mesh", radius = 1),
               "no faces")
  expect_error(vcg_detect_collision(roi, pts, mode_y = "mesh", radius = 1),
               "no faces")

  expect_error(vcg_detect_collision(roi, pts, radius = -1), "radius")
  expect_error(vcg_detect_collision(roi, pts, radius = NA), "radius")
  expect_error(vcg_detect_collision(roi, pts, radius = Inf), "radius")

  # 2-column input is padded to 3D, as elsewhere in the package.
  flat <- vcg_detect_collision(matrix(c(0, 0), ncol = 2L),
                               matrix(c(3, 4), ncol = 2L), radius = 10)
  expect_equal(flat$distance, 5)
})


test_that("vcg_detect_collision handles degenerate input", {
  roi <- matrix(rnorm(30), ncol = 3L)

  # Empty target.
  empty <- vcg_detect_collision(roi, matrix(numeric(0), ncol = 3L), radius = 1)
  expect_length(empty$hit, 0L)
  expect_length(empty$distance, 0L)

  # Empty ROI: nothing can collide.
  none <- vcg_detect_collision(matrix(numeric(0), ncol = 3L), roi, radius = 1e6)
  expect_false(any(none$hit))

  # A single row in `segments` mode carries no segment at all.
  lone <- vcg_detect_collision(roi, rbind(c(0, 0, 0)), mode_y = "segments",
                               radius = 10)
  expect_true(is.na(lone$hit))

  # A zero-length segment behaves like the point it degenerates to.
  degenerate <- vcg_detect_collision(rbind(c(0, 0, 0)),
                                     rbind(c(1, 0, 0), c(1, 0, 0)),
                                     mode_y = "segments", radius = 10)
  expect_equal(degenerate$distance[[1L]], 1, tolerance = 1e-6)

  # Coincident points are at distance 0, and radius 0 still counts as contact.
  touching <- vcg_detect_collision(rbind(c(1, 2, 3)), rbind(c(1, 2, 3)),
                                   radius = 0)
  expect_true(touching$hit)
  expect_equal(touching$distance, 0)
})
