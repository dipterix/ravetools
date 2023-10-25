test_that("Quaternion", {

  library(testthat)
  library(ravetools)
  q <- new_quaternion()

  expect_equal(q$to_array(), c(0, 0, 0, 1))

})
