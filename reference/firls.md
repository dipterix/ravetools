# Least-squares linear-phase `FIR` filter design

Produce a linear phase filter from the weighted mean squared such that
error in the specified bands is minimized.

## Usage

``` r
firls(N, freq, A, W = NULL, ftype = "", legacy = FALSE)
```

## Arguments

- N:

  filter order, must be even (if odd, then will be increased by one)

- freq:

  vector of frequency points in the range from 0 to 1, where 1
  corresponds to the `Nyquist` frequency.

- A:

  vector of the same length as `freq` containing the desired amplitude
  at each of the points specified in `freq`.

- W:

  weighting function that contains one value for each band that weights
  the mean squared error in that band. `W` must be half the length of
  `freq`.

- ftype:

  transformer type; default is `""`; alternatively, `'h'` or `'hilbert'`
  for 'Hilbert' transformer.

- legacy:

  whether to use the legacy implementations, which uses
  [`qr.solve`](https://rdrr.io/r/base/qr.html) instead of faster
  [`chol`](https://rdrr.io/r/base/chol.html)

## Value

The `FIR` filter coefficients with class `'Arma'`. The moving average
coefficient is a vector of length `n+1`.
