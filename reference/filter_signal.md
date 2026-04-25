# Filter one-dimensional signal

The function is written from the scratch. The result has been compared
against the 'Matlab' `filter` function with one-dimensional real inputs.
Other situations such as matrix `b` or multi-dimensional `x` are not
implemented. For double filters (forward-backward), see
[`filtfilt`](https://dipterix.org/ravetools/reference/filtfilt.md).

## Usage

``` r
filter_signal(b, a, x, z)
```

## Arguments

- b:

  one-dimensional real numerical vector, the moving-average coefficients
  of an `ARMA` filter

- a:

  the auto-regressive (recursive) coefficients of an `ARMA` filter

- x:

  numerical vector input (real value)

- z:

  initial condition, must have length of `n-1`, where `n` is the maximum
  of lengths of `a` and `b`; default is all zeros

## Value

A list of two vectors: the first vector is the filtered signal; the
second vector is the final state of `z`

## Examples

``` r

t <- seq(0, 1, by = 0.01)
x <- sin(2 * pi * t * 2.3)
bf <- gsignal::butter(2, c(0.15, 0.3))

res <- filter_signal(bf$b, bf$a, x)
y <- res[[1]]
z <- res[[2]]

## Matlab (2022a) equivalent:
# t = [0:0.01:1];
# x = sin(2 * pi * t * 2.3);
# [b,a] = butter(2,[.15,.3]);
# [y,z] = filter(b, a, x)

```
