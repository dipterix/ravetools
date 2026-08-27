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
    res$representation$distance[[1L]]
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
  expect_equal(res$representation$distance[[1L]], 3, tolerance = 1e-5)

  # The case the tractography viewer gets wrong: a segment passing straight
  # through a single-point ROI, with both of its vertices offset. Vertex
  # sampling reports 0.781; the true segment distance is 0.
  roi <- rbind(c(0, 0, 0))
  seg <- rbind(c(-0.781, 0, 0), c(0.781, 0, 0))
  res <- vcg_detect_collision(roi, seg, mode_y = "segments", radius = 1)
  expect_true(res$collide)
  expect_equal(res$representation$distance[[1L]], 0, tolerance = 1e-6)
})


test_that("vcg_detect_collision agrees with vcg_kdtree_nearest on point clouds", {
  set.seed(7)
  target <- matrix(rnorm(300), ncol = 3L)
  query <- matrix(rnorm(600), ncol = 3L)

  for (radius in c(0.25, 0.5, 1, 2)) {
    res <- vcg_detect_collision(target, query, mode_x = "points",
                                mode_y = "points", radius = radius)
    knn <- vcg_kdtree_nearest(target = target, query = query, k = 1L)

    # Every row of a point cloud is its own unit.
    expect_identical(as.logical(res$hit_unit), as.vector(knn$distance <= radius))
    expect_identical(res$collide, any(knn$distance <= radius))

    # The hit table carries one row per hit unit, and both the distance and the
    # reported target row must match the nearest neighbour.
    rep <- res$representation
    expect_identical(rep$unit, which(res$hit_unit))
    expect_equal(rep$distance, as.vector(knn$distance)[rep$unit], tolerance = 1e-5)
    expect_identical(rep$x_index, as.integer(knn$index)[rep$unit])

    # A point has exactly one element: itself.
    expect_true(all(rep$index == 1L))
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

    expect_identical(forward$collide, reverse$collide)
    expect_equal(min(forward$representation$distance),
                 min(reverse$representation$distance), tolerance = 1e-5)
  }
})


test_that("vcg_detect_collision accepts NA rows as polyline separators", {
  set.seed(3)
  roi <- matrix(rnorm(30), ncol = 3L)

  tract1 <- cbind(seq(-2, 2, by = 0.4), 0, 0)
  tract2 <- cbind(0, seq(-2, 2, by = 0.4), 1)
  combined <- rbind(tract1, c(NA, NA, NA), tract2)

  joint <- vcg_detect_collision(roi, combined, mode_y = "segments", radius = 1)

  # Two chains in one matrix are two units, whatever their length.
  expect_length(joint$hit_unit, 2L)

  # Splitting the matrix into one call per polyline must give the same answers,
  # and the segment index must be relative to each chain's own start.
  solo1 <- vcg_detect_collision(roi, tract1, mode_y = "segments", radius = 1)
  solo2 <- vcg_detect_collision(roi, tract2, mode_y = "segments", radius = 1)

  expect_identical(joint$hit_unit, c(solo1$hit_unit, solo2$hit_unit))

  for (k in seq_len(2L)) {
    solo <- list(solo1, solo2)[[k]]
    mine <- joint$representation[joint$representation$unit == k, , drop = FALSE]
    expect_identical(nrow(mine), nrow(solo$representation))
    if (nrow(mine)) {
      expect_identical(mine$index, solo$representation$index)
      expect_equal(mine$distance, solo$representation$distance)
      expect_identical(mine$x_index, solo$representation$x_index)
    }
  }

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
  expect_false(separated$collide)   # nearest real segment is 10 away
  expect_true(bridged$collide)      # without the separator, the bridge passes through

  # Rows that are only partly finite are treated as separators too, so this one
  # splits the first tract in half and yields a third unit.
  partial <- combined
  partial[3L, 2L] <- NA_real_
  expect_length(
    vcg_detect_collision(roi, partial, mode_y = "segments", radius = 1)$hit_unit,
    3L
  )
})


test_that("vcg_detect_collision counts every chain a separator delimits", {
  roi <- matrix(rnorm(30), ncol = 3L)
  A <- cbind(c(0, 1), 0, 0)
  B <- cbind(c(0, 1), 5, 0)
  sep <- c(NA, NA, NA)

  # Every slot is a chain, empty ones included. A leading separator opens an
  # empty chain; a trailing one closes the chain before it.
  cases <- list(
    list(y = A, n = 1L),
    list(y = rbind(A, sep, B), n = 2L),
    list(y = rbind(A, sep, sep, B), n = 3L),
    list(y = rbind(A, sep, sep, sep, B), n = 4L),
    list(y = rbind(sep, A), n = 2L),
    list(y = rbind(A, sep), n = 1L),
    list(y = rbind(sep), n = 1L),
    list(y = matrix(numeric(0), ncol = 3L), n = 0L)
  )

  for (case in cases) {
    res <- vcg_detect_collision(roi, case$y, mode_y = "segments", radius = 1)
    expect_length(res$hit_unit, case$n)
    expect_identical(res$summary$y$n_units, case$n)
  }

  # An empty chain must not shift the answers of the chains around it: here the
  # colliding chain is the third, not the second.
  origin <- rbind(c(0, 0, 0))
  away <- cbind(c(-5, -4), 20, 0)
  through <- cbind(c(-1, 1), 0, 0)

  shifted <- vcg_detect_collision(origin, rbind(away, sep, sep, through),
                                  mode_y = "segments", radius = 0.5)
  expect_identical(shifted$hit_unit, c(FALSE, NA, TRUE))
  expect_identical(shifted$representation$unit, 3L)

  # The rule exists so that appending a separator to every chain gives exactly
  # one unit per chain, wherever the empty ones fall.
  none <- matrix(numeric(0), ncol = 3L)
  bundle_of <- function(tracts) {
    do.call("rbind", lapply(tracts, function(line) rbind(line, sep)))
  }

  for (tracts in list(list(A, B, A), list(none, B, A),
                      list(A, none, A), list(A, B, none))) {
    res <- vcg_detect_collision(roi, bundle_of(tracts), mode_y = "segments",
                                radius = 1)
    expect_length(res$hit_unit, 3L)
    empty <- vapply(tracts, nrow, 0L) == 0L
    expect_true(all(is.na(res$hit_unit[empty])))
    expect_false(any(is.na(res$hit_unit[!empty])))
  }

  # Grouping applies to the indexed side too, where an empty chain simply
  # contributes no segment.
  indexed <- vcg_detect_collision(rbind(A, sep, sep, B), rbind(c(0.5, 0, 0)),
                                  mode_x = "segments", radius = 0.1)
  expect_true(indexed$collide)
  expect_identical(indexed$summary$x$n_units, 3L)
})


test_that("vcg_detect_collision keeps NA distinct from FALSE in hit_unit", {
  # NA means "no geometry here to test"; FALSE means "tested and missed". This
  # test pins that decision: flattening NA to FALSE anywhere in the pipeline
  # must fail here. See the "Why hit_unit reports NA" section of the docs.
  origin <- rbind(c(0, 0, 0))
  sep <- c(NA, NA, NA)
  near <- rbind(c(-0.1, 0, 0), c(0.1, 0, 0))   # collides with the origin
  far <- rbind(c(9, 9, 9), c(9, 9, 8))         # does not

  # points: a separator row is not a point, so it is neither a hit nor a miss.
  pts <- vcg_detect_collision(origin, rbind(c(0, 0, 0), sep, c(9, 9, 9)),
                              mode_y = "points", radius = 0.5)
  expect_identical(pts$hit_unit, c(TRUE, NA, FALSE))

  # segments: an empty chain and a single-vertex chain both carry no segment,
  # and sit alongside a genuine miss that must stay FALSE.
  chains <- vcg_detect_collision(
    origin,
    rbind(near, sep, sep, rbind(c(5, 5, 5)), sep, far),
    mode_y = "segments", radius = 0.5
  )
  expect_identical(chains$hit_unit, c(TRUE, NA, NA, FALSE))

  # mesh: a face is always testable, so no NA can appear.
  a <- vcg_sphere()
  b <- vcg_sphere()
  b$vb[1L, ] <- b$vb[1L, ] + 1
  expect_false(anyNA(vcg_detect_collision(a, b, radius = 0)$hit_unit))

  # The documented escape hatch produces the flattened mask on demand.
  expect_identical(chains$hit_unit %in% TRUE, c(TRUE, FALSE, FALSE, FALSE))

  # NA never leaks into the summaries: `collide` and `representation` already
  # give the NA-free view, and `which` drops NA on its own.
  expect_true(chains$collide)
  expect_identical(chains$representation$unit, which(chains$hit_unit))
  expect_identical(nrow(chains$representation),
                   sum(chains$hit_unit, na.rm = TRUE))
})


test_that("vcg_detect_collision summarises the query", {
  set.seed(17)
  roi <- vcg_sphere()
  pts <- matrix(rnorm(30), ncol = 3L)

  res <- vcg_detect_collision(roi, pts, radius = 1)
  s <- res$summary

  expect_named(s, c("x", "y", "radius", "test_level", "include_interior"))

  # The resolved mode, never "auto".
  expect_identical(s$x$mode, "mesh")
  expect_identical(s$y$mode, "points")
  expect_identical(s$x$unit_type, "face")
  expect_identical(s$y$unit_type, "point")
  expect_identical(s$x$n_units, ncol(roi$it))
  expect_identical(s$y$n_units, nrow(pts))
  expect_identical(s$y$n_units, length(res$hit_unit))

  expect_identical(s$radius, 1)
  expect_identical(s$test_level, "element")
  expect_false(s$include_interior)

  # A chain is a "line" on either side, and a separator row is still a unit in
  # points mode, so the count keeps matching `hit_unit`.
  chained <- rbind(pts[1:4, ], c(NA, NA, NA), pts[5:8, ])
  seg <- vcg_detect_collision(chained, chained, mode_x = "segments",
                              mode_y = "segments", radius = 1)
  expect_identical(seg$summary$x$unit_type, "line")
  expect_identical(seg$summary$x$n_units, 2L)
  expect_identical(seg$summary$y$n_units, length(seg$hit_unit))

  cloud <- vcg_detect_collision(chained, chained, radius = 1)
  expect_identical(cloud$summary$y$n_units, nrow(chained))
  expect_identical(cloud$summary$y$n_units, length(cloud$hit_unit))

  # `whole` drops `hit_unit`, so the summary is the only record of the size.
  whole <- vcg_detect_collision(roi, pts, radius = 1, test_level = "whole")
  expect_null(whole$hit_unit)
  expect_identical(whole$summary$y$n_units, s$y$n_units)
  expect_identical(whole$summary$test_level, "whole")

  expect_true(
    vcg_detect_collision(roi, pts, radius = 1,
                         include_interior = TRUE)$summary$include_interior
  )
})


test_that("vcg_detect_collision reports every element at test_level = 'element'", {
  set.seed(5)
  roi <- matrix(rnorm(30), ncol = 3L)
  tracts <- rbind(
    cbind(seq(-2, 2, by = 0.25), 0, 0),
    c(NA, NA, NA),
    cbind(0.2, seq(-2, 2, by = 0.25), 0)
  )

  res <- vcg_detect_collision(roi, tracts, mode_y = "segments", radius = 0.8)

  # `element` is the default, and `collide` summarises `hit_unit`.
  expect_identical(res, vcg_detect_collision(roi, tracts, mode_y = "segments",
                                             radius = 0.8,
                                             test_level = "element"))
  expect_identical(res$collide, any(res$hit_unit, na.rm = TRUE))
  expect_identical(nrow(res$representation), sum(res$hit_unit, na.rm = TRUE))
  expect_named(res, c("collide", "hit_unit", "representation", "summary"))
  expect_named(res$representation, c("unit", "index", "distance", "x_index"))

  # A point cloud reports one unit per row, separators included, so `hit_unit`
  # stays aligned with the input matrix.
  origin <- rbind(c(0, 0, 0))
  cloud <- rbind(c(0, 0, 0), c(0.1, 0, 0), c(9, 9, 9), c(0.2, 0, 0))

  flat <- vcg_detect_collision(origin, cloud, mode_y = "points", radius = 0.5)
  expect_identical(flat$hit_unit, c(TRUE, TRUE, FALSE, TRUE))
  expect_identical(flat$representation$unit, c(1L, 2L, 4L))

  grouped <- vcg_detect_collision(origin,
                                  rbind(cloud[1:2, ], c(NA, NA, NA), cloud[3:4, ]),
                                  mode_y = "points", radius = 0.5)
  expect_identical(grouped$hit_unit, c(TRUE, TRUE, NA, FALSE, TRUE))

  # A mesh reports one unit per face, and names the face vertex nearest `x`.
  a <- vcg_sphere()
  b <- vcg_sphere()
  b$vb[1L, ] <- b$vb[1L, ] + 1
  mesh <- vcg_detect_collision(a, b, radius = 0)
  expect_length(mesh$hit_unit, ncol(b$it))
  expect_gt(nrow(mesh$representation), 1L)
  expect_true(all(mesh$representation$index %in% 1:3))
  expect_true(all(mesh$representation$x_index %in% seq_len(ncol(a$it))))

  # Nothing within reach: `collide` is FALSE and the table is empty.
  far <- vcg_detect_collision(origin, rbind(c(50, 50, 50)), radius = 0.5)
  expect_false(far$collide)
  expect_identical(nrow(far$representation), 0L)
})


test_that("vcg_detect_collision test_level = 'unit' stops inside each unit", {
  origin <- rbind(c(0, 0, 0))

  # A point and a face each ARE their own unit, so there is nothing to stop
  # short of: the two levels must agree exactly. Only `summary$test_level`,
  # which records the question rather than the answer, is allowed to differ.
  answer <- function(res) res[c("collide", "hit_unit", "representation")]

  set.seed(5)
  cloud <- matrix(rnorm(300), ncol = 3L)
  expect_identical(answer(vcg_detect_collision(origin, cloud, radius = 1)),
                   answer(vcg_detect_collision(origin, cloud, radius = 1,
                                               test_level = "unit")))

  a <- vcg_sphere()
  b <- vcg_sphere()
  b$vb[1L, ] <- b$vb[1L, ] + 1
  expect_identical(answer(vcg_detect_collision(a, b, radius = 0)),
                   answer(vcg_detect_collision(a, b, radius = 0,
                                               test_level = "unit")))

  # A polyline does have several elements. This chain grazes the origin at 0.5
  # on its first segment before returning at 0.1 on its third, so the two levels
  # report the same unit but a different representative element.
  chain <- rbind(c(-2, 0.5, 0), c(2, 0.5, 0), c(2, 0.1, 0), c(-2, 0.1, 0))

  by_element <- vcg_detect_collision(origin, chain, mode_y = "segments",
                                     radius = 0.6)
  by_unit <- vcg_detect_collision(origin, chain, mode_y = "segments",
                                  radius = 0.6, test_level = "unit")

  expect_identical(by_element$hit_unit, by_unit$hit_unit)
  expect_identical(by_element$collide, by_unit$collide)

  # `element` finds the closest segment ...
  expect_identical(by_element$representation$index, 3L)
  expect_equal(by_element$representation$distance, 0.1, tolerance = 1e-5)
  # ... `unit` stops at the first one that qualifies.
  expect_identical(by_unit$representation$index, 1L)
  expect_equal(by_unit$representation$distance, 0.5, tolerance = 1e-5)

  # Segment indices are relative to the start of their own chain, not the
  # matrix, so the same chain repeated twice gives the same index twice.
  bundle <- rbind(chain, c(NA, NA, NA), chain)
  both <- vcg_detect_collision(origin, bundle, mode_y = "segments",
                               radius = 0.6, test_level = "unit")
  expect_identical(both$representation$unit, c(1L, 2L))
  expect_identical(both$representation$index, c(1L, 1L))

  # A chain of one vertex carries no segment at all: nothing to test.
  lone <- vcg_detect_collision(origin, rbind(c(0, 0, 0)), mode_y = "segments",
                               radius = 1, test_level = "unit")
  expect_identical(lone$hit_unit, NA)
  expect_false(lone$collide)
})


test_that("vcg_detect_collision test_level = 'whole' answers yes or no", {
  origin <- rbind(c(0, 0, 0))
  chain <- rbind(c(-2, 0.5, 0), c(2, 0.5, 0), c(2, 0.1, 0), c(-2, 0.1, 0))

  hit <- vcg_detect_collision(origin, chain, mode_y = "segments", radius = 0.6,
                              test_level = "whole")

  # No per-unit vector at this level -- just the answer and one witness.
  expect_named(hit, c("collide", "representation", "summary"))
  expect_type(hit$collide, "logical")
  expect_length(hit$collide, 1L)
  expect_identical(nrow(hit$representation), 1L)

  # Which unit wins is timing-dependent across threads, so assert membership in
  # the full scan's hit set rather than its identity.
  full <- vcg_detect_collision(origin, chain, mode_y = "segments", radius = 0.6)
  expect_true(hit$representation$unit %in% which(full$hit_unit))

  miss <- vcg_detect_collision(origin, chain, mode_y = "segments", radius = 0.01,
                               test_level = "whole")
  expect_false(miss$collide)
  expect_identical(nrow(miss$representation), 0L)

  # `collide` agrees with a full scan across all three modes.
  a <- vcg_sphere()
  near <- vcg_sphere()
  near$vb[1L, ] <- near$vb[1L, ] + 1
  far <- vcg_sphere()
  far$vb[1L, ] <- far$vb[1L, ] + 10

  for (b in list(near, far)) {
    expect_identical(
      vcg_detect_collision(a, b, radius = 0, test_level = "whole")$collide,
      vcg_detect_collision(a, b, radius = 0)$collide
    )
  }
})


test_that("vcg_detect_collision accepted set grows monotonically with radius", {
  set.seed(13)
  roi <- matrix(rnorm(60), ncol = 3L)
  query <- matrix(rnorm(300), ncol = 3L)

  previous <- integer(0)
  for (radius in c(0.1, 0.3, 0.6, 1, 1.5, 3)) {
    hit <- which(vcg_detect_collision(roi, query, radius = radius)$hit_unit)
    expect_true(all(previous %in% hit))
    previous <- hit
  }
  expect_gt(length(previous), 0L)

  # A radius that spans the whole bounding box must accept everything.
  span <- vcg_detect_collision(roi, query, radius = 1e6)
  expect_true(all(span$hit_unit))
})


test_that("vcg_detect_collision handles mesh-to-mesh proximity", {
  a <- vcg_sphere()
  b <- vcg_sphere()
  b$vb[1L, ] <- b$vb[1L, ] + 3

  # Two unit spheres 3 apart leave a surface gap of about 1.
  expect_false(vcg_detect_collision(a, b, radius = 0.9)$collide)
  expect_true(vcg_detect_collision(a, b, radius = 1.1)$collide)

  # Overlapping spheres genuinely intersect, so some faces are at distance 0.
  c <- vcg_sphere()
  c$vb[1L, ] <- c$vb[1L, ] + 1
  overlap <- vcg_detect_collision(a, c, radius = 0)
  expect_true(overlap$collide)
  expect_equal(min(overlap$representation$distance), 0, tolerance = 1e-6)

  # The reported index must point at a real face of `a`.
  expect_true(all(overlap$representation$x_index >= 1L &
                    overlap$representation$x_index <= ncol(a$it)))
})


test_that("vcg_detect_collision can count the interior of a closed surface", {
  roi <- vcg_sphere()
  centre <- rbind(c(0, 0, 0))

  expect_false(vcg_detect_collision(roi, centre, radius = 0.1)$collide)

  inside <- vcg_detect_collision(roi, centre, radius = 0.1,
                                 include_interior = TRUE)
  expect_true(inside$collide)
  expect_equal(inside$representation$distance, 0)
  # An interior hit touches no element of `x`, so there is none to name.
  expect_identical(inside$representation$x_index, NA_integer_)

  # The interior test applies at every level.
  for (level in c("element", "unit", "whole")) {
    expect_true(vcg_detect_collision(roi, centre, radius = 0.1,
                                     test_level = level,
                                     include_interior = TRUE)$collide)
  }

  # A point well outside stays a miss either way.
  outside <- rbind(c(10, 10, 10))
  expect_false(vcg_detect_collision(roi, outside, radius = 0.1,
                                    include_interior = TRUE)$collide)

  # A segment that lies wholly inside is only caught with the interior test.
  seg <- rbind(c(-0.2, 0, 0), c(0.2, 0, 0))
  expect_false(vcg_detect_collision(roi, seg, mode_y = "segments",
                                    radius = 0.1)$collide)
  expect_true(vcg_detect_collision(roi, seg, mode_y = "segments", radius = 0.1,
                                   include_interior = TRUE)$collide)

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

  expect_error(vcg_detect_collision(roi, pts, radius = 1, test_level = "group"),
               "'arg' should be one of")

  # 2-column input is padded to 3D, as elsewhere in the package.
  flat <- vcg_detect_collision(matrix(c(0, 0), ncol = 2L),
                               matrix(c(3, 4), ncol = 2L), radius = 10)
  expect_equal(flat$representation$distance, 5)
})


test_that("vcg_detect_collision handles degenerate input", {
  roi <- matrix(rnorm(30), ncol = 3L)

  # Empty target.
  empty <- vcg_detect_collision(roi, matrix(numeric(0), ncol = 3L), radius = 1)
  expect_length(empty$hit_unit, 0L)
  expect_identical(nrow(empty$representation), 0L)
  expect_false(empty$collide)

  # Empty ROI: nothing can collide.
  none <- vcg_detect_collision(matrix(numeric(0), ncol = 3L), roi, radius = 1e6)
  expect_false(none$collide)

  # A single row in `segments` mode carries no segment at all.
  lone <- vcg_detect_collision(roi, rbind(c(0, 0, 0)), mode_y = "segments",
                               radius = 10)
  expect_identical(lone$hit_unit, NA)

  # A zero-length segment behaves like the point it degenerates to.
  degenerate <- vcg_detect_collision(rbind(c(0, 0, 0)),
                                     rbind(c(1, 0, 0), c(1, 0, 0)),
                                     mode_y = "segments", radius = 10)
  expect_equal(degenerate$representation$distance, 1, tolerance = 1e-6)

  # Coincident points are at distance 0, and radius 0 still counts as contact.
  touching <- vcg_detect_collision(rbind(c(1, 2, 3)), rbind(c(1, 2, 3)),
                                   radius = 0)
  expect_true(touching$collide)
  expect_equal(touching$representation$distance, 0)
})
