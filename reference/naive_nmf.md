# A naive implementation of non-negative matrix factorization

A pure-R vanilla implementation assuming inputs are non-negative
matrices without `NA`.

## Usage

``` r
naive_nmf(x, k, tol = c(1e-04, 1e-08), max_iters = 10000, verbose = TRUE)
```

## Arguments

- x:

  a matrix, or can be converted into a matrix; all negative or missing
  values will be treated as zero

- k:

  decomposition rank

- tol:

  stop criteria, a numeric of two; the first number is the tolerance for
  root-mean-squared residuals, relative to the largest number in `x`;
  the second number is the tolerance for weight differences; any
  stopping criteria met will result in the stop of iteration

- max_iters:

  maximum iterations

- verbose:

  whether to report the progress; logical or a positive integer (of step
  intervals)

## Value

A list of weights (non-negative template matrix `W` and non-negative
`H`) and errors (root mean squared error of fitted, matrix `W`, and `W`
versus their previous iteration, respectively).

## Examples

``` r



x <- stats::toeplitz(.9 ^ (0:31))

nmf <- naive_nmf(x, k = 7, verbose = FALSE)

fitted <- nmf$W %*% nmf$H

oldpar <- par(mfrow = c(1, 2))
on.exit({ par(oldpar )})

image(x, zlim = c(0, 1), main = "Input")

image(fitted, zlim = c(0, 1),
      main = sprintf("Fitted with rank=%d", nmf$k))




```
