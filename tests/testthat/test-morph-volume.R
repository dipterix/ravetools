# Foreground cube spanning `from:to` along each axis of a `size`-cube volume
morph_test_cube <- function(size, from, to) {
  x <- array(0L, dim = rep(size, 3))
  x[from:to, from:to, from:to] <- 1L
  x
}

test_that("morph_volume dilates with a ball whose radius is honored", {
  # a voxel grows to itself, 6 faces and 12 edges (the corners are too far)
  x <- array(0L, dim = c(11, 11, 11))
  x[6, 6, 6] <- 1L
  dilated <- morph_volume(x, "dilate", radius = 1)
  testthat::expect_equal(sum(dilated), 19)
  testthat::expect_equal(dilated[7, 6, 6], 1L)
  testthat::expect_equal(dilated[7, 7, 6], 1L)
  testthat::expect_equal(dilated[7, 7, 7], 0L)

  # integer points with x^2 + y^2 + z^2 <= 6: 1 + 6 + 12 + 8 + 6 + 24 + 24
  testthat::expect_equal(sum(morph_volume(x, "dilate", radius = 2)), 81)
})

test_that("morph_volume takes a radius per axis", {
  # a voxel grows to the ellipsoid with semi-axes 2.5, 1.5 and 0.5 voxels:
  # 5 voxels along the first axis and 3 along the second, in its own slice;
  # beside the 5 voxels of the long axis, the rows next to it hold 3 each
  x <- array(0L, dim = c(11, 11, 11))
  x[6, 6, 6] <- 1L
  dilated <- morph_volume(x, "dilate", radius = c(2, 1, 0))
  testthat::expect_equal(which(apply(dilated, 1, any)), 4:8)
  testthat::expect_equal(which(apply(dilated, 2, any)), 5:7)
  testthat::expect_equal(which(apply(dilated, 3, any)), 6L)
  testthat::expect_equal(sum(dilated), 11)

  # the same radius on every axis is the ball
  testthat::expect_equal(morph_volume(x, "dilate", radius = c(2, 2, 2)),
                         morph_volume(x, "dilate", radius = 2))
  cube <- morph_test_cube(9, 2, 8)
  testthat::expect_equal(morph_volume(cube, "erode", radius = c(1, 1, 1)),
                         morph_volume(cube, "erode", radius = 1))

  # an erosion along the third axis only: the 5-voxel cube loses its first
  # and last slices
  eroded <- morph_volume(morph_test_cube(9, 3, 7), "erode", radius = c(0, 0, 1))
  testthat::expect_equal(sum(eroded), 75)
  testthat::expect_equal(which(apply(eroded, 3, any)), 4:6)
  testthat::expect_equal(which(apply(eroded, 1, any)), 3:7)
})

test_that("morph_volume erodes the voxels next to the background", {
  # a 5-voxel cube loses its outer layer and keeps the 3-voxel core
  x <- morph_test_cube(9, 3, 7)
  eroded <- morph_volume(x, "erode", radius = 1)
  testthat::expect_equal(sum(eroded), 27)
  testthat::expect_equal(eroded[5, 5, 5], 1L)
  testthat::expect_equal(eroded[3, 5, 5], 0L)
})

test_that("morph_volume does not erode at the edge of the image", {
  x <- array(1L, dim = c(5, 5, 5))
  testthat::expect_equal(sum(morph_volume(x, "erode", radius = 2)), 125)
})

test_that("morph_volume follows the axes of a volume that is not a cube", {
  x <- array(0L, dim = c(4, 5, 6))
  x[1, 2, 3] <- 1L
  dilated <- morph_volume(x, "dilate", radius = 1)
  testthat::expect_equal(dim(dilated), c(4L, 5L, 6L))
  testthat::expect_equal(dilated[2, 2, 3], 1L)
  testthat::expect_equal(dilated[1, 3, 3], 1L)
  testthat::expect_equal(dilated[1, 2, 4], 1L)
  # the 5 neighbors on the side of `x - 1` fall out of the volume
  testthat::expect_equal(sum(dilated), 14)
})

test_that("morph_volume keeps the largest component and leaves its holes", {
  x <- morph_test_cube(12, 2, 6)
  x[4, 4, 4] <- 0L
  x[9:10, 9:10, 9:10] <- 1L
  largest <- morph_volume(x, "largest_component")
  testthat::expect_equal(sum(largest), 124)
  testthat::expect_equal(largest[4, 4, 4], 0L)
  testthat::expect_equal(largest[9, 9, 9], 0L)
})

test_that("morph_volume fills enclosed holes and leaves open cavities", {
  x <- morph_test_cube(9, 3, 7)
  x[5, 5, 5] <- 0L
  testthat::expect_equal(morph_volume(x, "fill_holes"), morph_test_cube(9, 3, 7))

  # a cavity that reaches the background through a one-voxel channel
  x[5, 5, 6:7] <- 0L
  testthat::expect_equal(sum(morph_volume(x, "fill_holes")), 125 - 3)
})

test_that("morph_volume takes voxels above zero as the foreground", {
  x <- array(0, dim = c(5, 5, 5))
  x[1, 1, 1] <- 0.5
  x[2, 1, 1] <- 2
  x[3, 1, 1] <- -1
  x[4, 1, 1] <- NA
  out <- morph_volume(x, "dilate", radius = 0)
  testthat::expect_equal(out[1:4, 1, 1], c(1L, 1L, 0L, 0L))
  testthat::expect_equal(sum(out), 2)
})

test_that("morph_volume refuses inputs it cannot process", {
  testthat::expect_error(morph_volume(matrix(1, 3, 3), "erode"), "3-dimensional")
  x <- array(1, dim = c(3, 3, 3))
  testthat::expect_error(morph_volume(x, "erode", radius = -1), "radius")
  testthat::expect_error(morph_volume(x, "erode", radius = NA), "radius")
  testthat::expect_error(morph_volume(x, "erode", radius = c(1, 1)), "radius")
  testthat::expect_error(morph_volume(x, "erode", radius = c(1, -1, 1)), "radius")
  testthat::expect_error(morph_volume(x, "erode", radius = c(1, NA, 1)), "radius")
  testthat::expect_error(morph_volume(x, "shrink"))
})
