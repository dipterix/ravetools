# Frequency response of digital filter

Compute the z-plane frequency response of an `ARMA` model.

## Usage

``` r
freqz2(b, a = 1, fs = 2 * pi, n = 512, whole = FALSE, ...)
```

## Arguments

- b:

  the moving-average coefficients of an `ARMA` model

- a:

  the auto-regressive coefficients of an `ARMA` filter; default is `1`

- fs:

  sampling frequency in `Hz`

- n:

  number of points at which to evaluate the frequency response; default
  is `512`

- whole:

  whether to evaluate beyond `Nyquist` frequency; default is false

- ...:

  ignored

## Value

A list of frequencies and corresponding responses in complex vector
