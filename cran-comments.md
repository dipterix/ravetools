## R CMD check results

0 errors | 0 warnings | 0 note


### Notes, fixes

* CRAN reports gcc-UBSAN errors

The error was because of two factors:
1. The Intel TBB library used by `RcppParallel` package failed the check
2. There was a bug in my code: `col2Idx = (*col2_ptr) - 1;` where `*col2_ptr` is an `int` type, this statement overflows when `*col2_ptr` is `NA_INTEGER`

To solve these issues:
1. I extracted non-TBB part from `RcppParallel` and put them into `inst/include` folder so the package is independent from `RcppParallel`. (The authorship and credit has been given properly in both `README` and `DESCRIPTION` files)
2. The integer overflow has been resolved by checking `*col2_ptr` with conditions `R_finite(*col2_ptr) && *col2_ptr >= 1`, such that `NA_INTEGER` enters the `else` clause.

