test_that("VCG functions", {

  n <- 20
  A <- matrix(rnorm(n * 3), nrow = n)
  B <- matrix(rnorm(n * 6), nrow = n * 2)
  result <- vcg_kdtree_nearest(
    target = B, query = A,
     k = 1
  )
  nearest_index <- apply(A, 1, function(pt) {
    which.min(colSums((t(B) - pt) ^ 2))
  })

  expect_true(all(result$index == nearest_index))

})
