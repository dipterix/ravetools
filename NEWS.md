# ravetools 0.0.2

This is an initial version of `ravetools`. Although a bare minimal set of signal processing functions are provided, it is sufficient to perform preprocess pipelines on most `iEEG` signals. Some functions are added from the `dipsaus` package, with considerable performance improvement. The `C++` functions have been tested on all major platforms, with different architectures (`ARM`, `i386`, `x64`).

### Documentation

* Added `README` file to demonstrate basic usage
* Added a `NEWS.md` file to track changes to the package.

### Signal processing functions
* Re-implemented `decimate` with `FIR` filters creating the same results as in `Matlab`
* Added `detrend` function to 
* Added `diagnose_channel` to visually inspect channel signals
* Added `morlet_wavelet` to enable fast and memory efficient wavelet decomposition; the result agrees with existing `Matlab` code with floating errors (`10^-7`)
* Added `multitaper`
* Added `pwelch` (`Welch` periodogram)
* Added `notch_filter` to remove line noise

### High-performance functions

The following functions are implemented in `C++` parallel. They tend to be faster than normal base-R implementations, depending on the number of `CPU` cores used.

* Added `collapse` to collapse arrays
* Added `shift_array` to shift array along certain indices
* Added `fast_cov` to calculate `pearson` covariance matrix in parallel
* Added `baseline_array` to calculate baseline arrays with multiple margins
