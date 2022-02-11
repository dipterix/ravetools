library(testthat)
library(ravetools)

# Make sure use 2 cores to comply with the CRAN policy
RcppParallel::setThreadOptions(numThreads = 2L)
test_check("ravetools")
