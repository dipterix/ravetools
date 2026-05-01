# Check 'Arma' filter

Check 'Arma' filter

## Usage

``` r
check_filter(b, a, w = NULL, r_expected = NULL, fs = NULL)
```

## Arguments

- b:

  moving average (`MA`) polynomial coefficients.

- a:

  auto-regressive (`AR`) polynomial coefficients.

- w:

  normalized frequency, ranging from 0 to 1, where 1 is 'Nyquist'

- r_expected:

  attenuation in decibel of each `w`

- fs:

  sample rate, used to infer the frequencies and formatting print
  message, not used in calculation; leave it blank by default

## Value

A list of power estimation and the reciprocal condition number of the
`AR` coefficients.

## Examples

``` r


# create a butterworth filter with -3dB (half-power) at [1, 5] Hz
# and -60dB stop-band attenuation at [0.5, 6] Hz

sample_rate <- 20
nyquist <- sample_rate / 2

specs <- buttord(
  Wp = c(1, 5) / nyquist,
  Ws = c(0.5, 6) / nyquist,
  Rp = 3,
  Rs = 60
)
filter <- butter(specs)

# filter quality is poor because the AR-coefficients
# creates singular matrix with unstable inverse,
# this will cause `filtfilt` to fail
check_filter(
  b = filter$b, a = filter$a,

  # frequencies (normalized) where power is evaluated
  w = c(1, 5, 0.5, 6) / nyquist,

  # expected power
  r_expected = c(3, 3, 60, 60)

)
#> <RAVE filter quality test>
#> Attenuation: 
#>   Freq=0.1 xPi rad/s, mag=-4.113 dB (expected=3 dB)
#>   Freq=0.5 xPi rad/s, mag=-3 dB (expected=3 dB)
#>   Freq=0.05 xPi rad/s, mag=-129.7 dB (expected=60 dB)
#>   Freq=0.6 xPi rad/s, mag=-63.23 dB (expected=60 dB)
#> Reciprocal condition number: 4.3e-21 < .Machine$double.eps
#> 
#> WARNING: 
#>  * Unstable autoregressive (AR) polynomial coefficients

```
