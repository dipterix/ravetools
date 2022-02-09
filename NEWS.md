# ravetools 0.0.1

* Added `README` file to demonstrate basic usage
* Added a `NEWS.md` file to track changes to the package.
* Re-implemented `decimate` with `FIR` filters creating the same results as in `Matlab`
* Added `detrend` function to 
* Added `diagnose_channel` to visually inspect channel signals
* Added `morlet_wavelet` to enable fast and memory efficient wavelet decomposition; the result agrees with existing `Matlab` code with floating errors (`10^-7`)
* Added `multitaper`
* Added `pwelch` (`Welch` periodogram)
* Added `notch_filter` to remove line noise

