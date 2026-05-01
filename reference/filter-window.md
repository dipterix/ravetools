# Filter window functions

Filter window functions

## Usage

``` r
hanning(n)

hamming(n)

blackman(n)

blackmannuttall(n)

blackmanharris(n)

flattopwin(n)

bohmanwin(n)
```

## Arguments

- n:

  number of time-points in window

## Value

A numeric vector of window with length `n`

## Examples

``` r

hanning(10)
#>  [1] 0.0000000 0.1169778 0.4131759 0.7500000 0.9698463 0.9698463 0.7500000
#>  [8] 0.4131759 0.1169778 0.0000000
hamming(11)
#>  [1] 0.0800000 0.1678522 0.3978522 0.6821478 0.9121478 1.0000000 0.9121478
#>  [8] 0.6821478 0.3978522 0.1678522 0.0800000
blackmanharris(21)
#>  [1] 0.000060000 0.001791203 0.010982331 0.039190758 0.103011489 0.217470000
#>  [7] 0.385892669 0.590993400 0.793833511 0.944304639 1.000000000 0.944304639
#> [13] 0.793833511 0.590993400 0.385892669 0.217470000 0.103011489 0.039190758
#> [19] 0.010982331 0.001791203 0.000060000
```
