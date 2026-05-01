# Low-level FFTW3 wrappers

Thin R bindings around the FFTW3 library. These are **low-level**
routines exposed primarily for advanced users and other packages that
need maximum throughput. They perform minimal input checking and follow
FFTW conventions (e.g. unnormalized inverse transforms, one-sided
real-to-complex spectra). For most user code prefer
[`fft`](https://rdrr.io/r/stats/fft.html),
[`mvfft`](https://rdrr.io/r/stats/fft.html), or higher-level helpers in
this package such as
[`convolve`](https://dipterix.org/ravetools/reference/convolve.md),
[`pwelch`](https://dipterix.org/ravetools/reference/pwelch.md),
[`multitaper`](https://dipterix.org/ravetools/reference/multitaper.md),
and the filtering utilities.

**Warning:** the API is intentionally close to FFTW's C interface and
may change between releases. Outputs match the corresponding base R
transforms up to floating-point round-off.

## Arguments

- data:

  Numeric (real) or complex input. For 2D/3D variants, a matrix or
  3-dimensional array. For `mvfftw_*` functions, a matrix whose columns
  are transformed independently.

- HermConj:

  Integer `0` or `1`. When `1`, return the full Hermitian-symmetric
  spectrum of length `N` (matching
  [`stats::fft`](https://rdrr.io/r/stats/fft.html)); when `0`, return
  only the non-redundant one-sided half of length `floor(N/2) + 1`.

- inverse:

  Integer `0` (forward) or `1` (inverse). The inverse transform is
  **unnormalized**; divide by the number of elements to invert a forward
  transform, mirroring `stats::fft(., inverse = TRUE)`.

- fftwplanopt:

  Integer planner effort: `0` = `FFTW_ESTIMATE` (default, fast
  planning), `1` = `FFTW_MEASURE`, `2` = `FFTW_PATIENT`, `3` =
  `FFTW_EXHAUSTIVE`.

- retrows:

  Integer; expected number of rows of the time-domain signal for
  `mvfftw_c2r` when reconstructing from a one-sided spectrum.

- ret:

  Optional pre-allocated output buffer of the correct type and length;
  pass `NULL` (default) to let the function allocate one.

## Value

A complex (or real, for `*_c2r`) vector / matrix / array matching the
corresponding base R transform up to floating-point error.

## Details

All functions preserve their `data` argument: the input buffer is copied
internally before planning when needed, so callers may safely reuse
`data` after the call. For multi-dimensional variants, axis ordering
follows R (column-major) conventions.

## Examples

``` r
set.seed(1)

## --- 1D real-to-complex --------------------------------------------------
x <- rnorm(16)
a <- ravetools::fftw_r2c(x, HermConj = 1)
b <- stats::fft(x)
all.equal(a, b)                                  # TRUE (within tol)
#> [1] TRUE

# one-sided spectrum (length floor(N/2)+1)
a_half <- ravetools::fftw_r2c(x, HermConj = 0)
all.equal(a_half, b[seq_len(length(x) %/% 2 + 1)])
#> [1] TRUE

## --- 1D complex-to-complex ----------------------------------------------
z <- complex(real = rnorm(16), imaginary = rnorm(16))
all.equal(ravetools::fftw_c2c(z, inverse = 0), stats::fft(z))
#> [1] TRUE
all.equal(ravetools::fftw_c2c(z, inverse = 1), stats::fft(z, inverse = TRUE))
#> [1] TRUE

## --- 1D complex-to-real (inverse of fftw_r2c) ---------------------------
# Using the full Hermitian spectrum:
xr <- ravetools::fftw_c2r(a, HermConj = 1) / length(x)
all.equal(xr, x)
#> [1] TRUE

## --- Multivariate (column-wise) ----------------------------------------
M <- matrix(rnorm(32), nrow = 8, ncol = 4)
all.equal(ravetools::mvfftw_r2c(M, HermConj = 1), stats::mvfft(M + 0i))
#> [1] TRUE

Mz <- matrix(complex(real = rnorm(32), imaginary = rnorm(32)),
             nrow = 8, ncol = 4)
all.equal(ravetools::mvfftw_c2c(Mz, inverse = 0), stats::mvfft(Mz))
#> [1] TRUE
all.equal(ravetools::mvfftw_c2c(Mz, inverse = 1),
          stats::mvfft(Mz, inverse = TRUE))
#> [1] TRUE

# one-sided -> back to real signal
Mh <- ravetools::mvfftw_r2c(M, HermConj = 0)
Mr <- ravetools::mvfftw_c2r(Mh, retrows = nrow(M)) / nrow(M)
all.equal(Mr, M)
#> [1] TRUE

## --- 2D ----------------------------------------------------------------
X2 <- matrix(rnorm(20), nrow = 5, ncol = 4)
all.equal(ravetools::fftw_r2c_2d(X2, HermConj = 1), stats::fft(X2 + 0i))
#> [1] TRUE

Z2 <- matrix(complex(real = rnorm(20), imaginary = rnorm(20)),
             nrow = 5, ncol = 4)
all.equal(ravetools::fftw_c2c_2d(Z2, inverse = 0), stats::fft(Z2))
#> [1] TRUE
all.equal(ravetools::fftw_c2c_2d(Z2, inverse = 1),
          stats::fft(Z2, inverse = TRUE))
#> [1] TRUE

## --- 3D ----------------------------------------------------------------
X3 <- array(rnorm(60), dim = c(5, 4, 3))
all.equal(ravetools::fftw_r2c_3d(X3, HermConj = 1), stats::fft(X3 + 0i))
#> [1] TRUE

Z3 <- array(complex(real = rnorm(60), imaginary = rnorm(60)),
            dim = c(5, 4, 3))
all.equal(ravetools::fftw_c2c_3d(Z3, inverse = 0), stats::fft(Z3))
#> [1] TRUE
all.equal(ravetools::fftw_c2c_3d(Z3, inverse = 1),
          stats::fft(Z3, inverse = TRUE))
#> [1] TRUE
```
