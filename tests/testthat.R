library(testthat)
library(ravetools)

# Make sure use 2 cores to comply with the CRAN policy
ravetools_threads(n_threads = 2L)
options(rgl.useNULL=TRUE)
test_check("ravetools")
