test_that("resample3D", {
  arrayDim <- rev(c(256, 128, 16))
  fromArray <- array(rnorm(prod(arrayDim)), rev(arrayDim) + 1)
  oldVoxToWorld <- matrix(nrow = 4, byrow = TRUE, c(
    1.5, 0, 0.8, -200,
    0.01, 0, -1, 10,
    0, -0.8, 0.1, 60,
    0, 0, 0, 1
  ))
  newVoxToWorld <- matrix(nrow = 4, byrow = TRUE, c(
    15, 0, 0.08, -140,
    0.1, -0.5, 0, 30,
    0, 0.1, 0.1, -13.6,
    0, 0, 0, 1
  ))

  # Make sure the resulting data has valid data so this test can be valid
  expect_true(all(solve(oldVoxToWorld) %*% newVoxToWorld %*% c(arrayDim / 2,1) > 0))
  expect_true(all((solve(oldVoxToWorld) %*% newVoxToWorld %*% c(arrayDim / 2,1))[1:3] < dim(fromArray)))

  reList <- resample3D(arrayDim, fromArray, t(newVoxToWorld), t(oldVoxToWorld), na = NA_real_)

  # reList[[2]] is the point transform of indexing space from new array to old array
  dim(reList[[2]]) <- c(4,4)
  mat_actual <- reList[[2]]
  mat_expected <- solve(oldVoxToWorld) %*% newVoxToWorld

  mat <- mat_actual %*% solve(mat_expected)
  diag(mat) <- diag(mat) - 1

  # should be about 1e-14, but 1e-5 should be fine
  expect_lt(max(abs(mat)), 1e-5)

  get_expected_array <- function(na_fill = NA_real_) {
    # compare with manual
    expected_array <- array(na_fill, arrayDim)

    # construct index set in ijk1 format
    idx <- t(cbind(arrayInd(seq_len(prod(arrayDim)), arrayDim) - 1, 1))

    # calculate index in old array
    idx <- round(solve(oldVoxToWorld) %*% newVoxToWorld %*% idx)[1:3, , drop = FALSE]

    # remove outofbound indexings
    idx[idx < 0 | idx >= dim(fromArray)] <- NA

    # get 1-index
    idx <- colSums(idx * c(1, cumprod(dim(fromArray)))[1:3]) + 1

    expected_array[!is.na(idx)] <- fromArray[idx[!is.na(idx)]]
    expected_array
  }

  expected_array <- get_expected_array(NA_real_)
  expect_equal(storage.mode(reList[[1]]), storage.mode(fromArray))
  expect_equal(reList[[1]], expected_array)



  # check if other formats work too

  # integer
  fromArray <- fromArray * 10
  storage.mode(fromArray) <- "integer"
  reList <- resample3D(arrayDim, fromArray, t(newVoxToWorld), t(oldVoxToWorld), na = NA_integer_)
  dim(reList[[2]]) <- c(4,4)
  expect_equal(reList[[2]], expected = mat_actual)
  expect_equal(storage.mode(reList[[1]]), "integer")
  expect_equal(reList[[1]], get_expected_array(NA_integer_))

  # string
  storage.mode(fromArray) <- "character"
  reList <- resample3D(arrayDim, fromArray, t(newVoxToWorld), t(oldVoxToWorld), na = NA_character_)
  dim(reList[[2]]) <- c(4,4)
  expect_equal(reList[[2]], expected = mat_actual)
  expect_equal(storage.mode(reList[[1]]), "character")
  expect_equal(reList[[1]], get_expected_array(NA_character_))

  # complex
  storage.mode(fromArray) <- "complex"
  reList <- resample3D(arrayDim, fromArray, t(newVoxToWorld), t(oldVoxToWorld), na = NA_complex_)
  dim(reList[[2]]) <- c(4,4)
  expect_equal(reList[[2]], expected = mat_actual)
  expect_equal(storage.mode(reList[[1]]), "complex")
  expect_equal(reList[[1]], get_expected_array(NA_complex_))

  # raw
  fromArray <- array(sample(255, size = length(fromArray), replace = TRUE), dim(fromArray))
  storage.mode(fromArray) <- "raw"
  reList <- resample3D(arrayDim, fromArray, t(newVoxToWorld), t(oldVoxToWorld), na = as.raw(0))
  dim(reList[[2]]) <- c(4,4)
  expect_equal(reList[[2]], expected = mat_actual)
  expect_equal(storage.mode(reList[[1]]), "raw")
  expect_equal(reList[[1]], get_expected_array(as.raw(0)))

  # logical
  fromArray <- array(sample(c(TRUE, FALSE, NA_integer_), size = length(fromArray), replace = TRUE), dim(fromArray))
  storage.mode(fromArray) <- "logical"
  reList <- resample3D(arrayDim, fromArray, t(newVoxToWorld), t(oldVoxToWorld), na = as.logical(NA))
  dim(reList[[2]]) <- c(4,4)
  expect_equal(reList[[2]], expected = mat_actual)
  expect_equal(storage.mode(reList[[1]]), "logical")
  expect_equal(reList[[1]], get_expected_array(as.logical(NA)))

})
