# Computer reciprocal condition number of an 'Arma' filter

Test whether the filter is numerically stable for
[`filtfilt`](https://dipterix.org/ravetools/reference/filtfilt.md).

## Usage

``` r
rcond_filter_ar(a)
```

## Arguments

- a:

  auto-regression coefficient, numerical vector; the first element must
  not be zero

## Value

Reciprocal condition number of matrix `z1`, used in
[`filtfilt`](https://dipterix.org/ravetools/reference/filtfilt.md). If
the number is less than `.Machine$double.eps`, then
[`filtfilt`](https://dipterix.org/ravetools/reference/filtfilt.md) will
fail.

## See also

[`check_filter`](https://dipterix.org/ravetools/reference/check_filter.md)

## Examples

``` r

# Butterworth filter with low-pass at 0.1 Hz (order = 4)
filter <- butter(4, 0.1, "low")

# TRUE
rcond_filter_ar(filter$a) > .Machine$double.eps
#> [1] TRUE

diagnose_filter(filter$b, filter$a, 500)


# Bad filter (order is too high)
filter <- butter(50, 0.1, "low")

rcond_filter_ar(filter$a) > .Machine$double.eps
#> [1] FALSE

# filtfilt needs to inverse a singular matrix
diagnose_filter(filter$b, filter$a, 500)

```
