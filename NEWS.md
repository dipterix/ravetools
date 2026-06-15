# ravetools 0.2.7

* Added `color_ramp_continuous` helper to create a vectorized color-mapping function from a palette and a numeric domain
* Fixed `make_error` to properly protect the error class string from garbage collection, preventing class names from being silently lost in nested error calls
* Guarded `design_filter` against `bandstop` configurations where automatically derived transition bandwidths overlap each other, producing unordered frequency breakpoints
* Added `mris_curvature` to estimate per-vertex mean, Gaussian, and principal curvatures by fitting an osculating quadratic surface to each vertex's two-ring neighborhood
* Added `mris_smooth` to smooth a triangular surface mesh via iterative `Laplacian` or `Taubin` smoothing
* Added `mris_inflate` to inflate a cortical surface mesh by iteratively smoothing and expanding toward a sphere while preserving total surface area
* Added `mris_sphere` to project an inflated cortical surface onto a sphere and relax metric distortion via area-preserving gradient descent
* Added `mris_remesh` to re-mesh a triangular surface to uniform edge length using the `Botsch` and `Kobbelt` (2004) isotropic re-meshing algorithm
* Added `mris_make_surfaces` to localize white-matter and `pial` surfaces from an intensity volume
* Added `vcg_average_edge_length` and `vcg_max_edge_length` to query edge-length statistics of a mesh
* Added `vcg_count_edge_defects` to count non-manifold edge defects in a mesh
* Added `vcg_fix_defects` to repair mesh topology defects via `VCG`
* Added `vcg_subdivide_max_edge_length` to subdivide only edges that exceed a specified length, avoiding global super-sampling and improving `Dijkstra` path quality
* Added `vcg_mesh_patch` to extract a geodesic patch of a mesh bounded by a closed curve of user-supplied `waypoint` vertices
* All `mesh3d` objects returned by `ravetools` now carry the `ravetools_mesh3d` class, enabling a `plot()` generic that dispatches automatically to `plot_mesh_polygon` or `plot_mesh_dotcloud`
* `register_volume3d` (`SyN` mode) re-implemented from the `ANTs` paper in pure `Rcpp`/`RcppEigen`; no longer depends on `RNiftyReg` or any external registration library
* `register_volume3d` gains a `points` argument to additionally warp a matrix of anatomical coordinates (e.g. cortical surface vertices) through the estimated transform, enabling direct cortical mapping without a separate resampling step
* Exported `apply_transform3d` to apply a previously estimated rigid, `affine`, or `SyN` transform to a new 3D volume


# ravetools 0.2.6

* Added `plot_mesh_dotcloud` for rendering one or more `mesh3d` objects as an orthographic rim-lit dot cloud in base R (no `rgl` dependency); supports per-vertex colors, depth-gradient palettes, per-mesh alpha, `side` filtering, and painter's-algorithm depth sorting across multiple meshes
* Added `plot_mesh_polygon` for rendering one or more `mesh3d` objects as flat-shaded Lambert-lit triangles in base R; point-cloud meshes (no face matrix) are automatically substituted by small sphere instances; supports camera-facing clipping (`mesh_clipping`), `side` filtering, per-mesh alpha blending, and configurable `shadow_color`, `light_intensity`, and `ambient_intensity`
* `plot_mesh_dotcloud` and `plot_mesh_polygon` gain a `clipping_plane` argument (a 4-element numeric normal-and-offset vector, or a list of such vectors) to discard geometry on one side of arbitrary planes before rendering, and a per-mesh `clipping_plane_enabled` logical toggle to opt individual meshes in or out of clipping
* Exported `ensure_mesh3d` for coercing various surface formats (`mesh3d`, `ieegio_surface`, `fs.surface`, `surf.asc`) to a canonical `mesh3d` object; previously internal-only
* Fixed `freqz2` plot legend: the cutoff legend used the raw `cutoffs` argument (and could index past the color palette, yielding `NA` colors) when a requested cutoff did not actually intersect the magnitude curve; the legend now only lists cutoffs that are drawn, and labels the unit as `dB` instead of the filter's frequency unit


# ravetools 0.2.5

* `morlet_wavelet` gains a `segment_length` argument (default `NULL`) that processes long signals (e.g. multi-hour recordings) in overlapping segments using batched `mvfftw_c2c` convolutions, dramatically reducing peak memory and FFT cost while preserving the legacy result on the signal interior
* `baseline_array` now supports seven baseline methods: `"percentage"`, `"sqrt_percentage"`, `"decibel"`, `"zscore"`, `"sqrt_zscore"`, `"db_zscore"`, and `"subtract_mean"`; newly added methods include decibel z-score normalization
* Added spike sorting helpers for normality test, spike detection, per-channel waveform extraction, and `Haar` wavelet feature extraction
* Replaced `wavelets` package dependency with `waveslim` for discrete wavelet transform used in spike sorting utilities
* `project_plane` gains an `n_iters` argument to control the number of projection iterations (default is `5`)
* Applied lint fixes and improved code style consistency across all source files
* Added `carla` implementing the Common Average Re-referencing by Least Anti-Correlation (`CARLA`) algorithm (see `CITATION`) for selecting an optimal subset of channels as the common average reference in `CCEP` data
* Added `crp` implementing the Canonical Response Parameterization (`CRP`) method for characterizing single-trial evoked responses (e.g. `CCEP`s): estimates response duration, extracts the canonical response shape, and reports per-trial weights, `SNR`, and explained variance
* Exposed low-level `FFTW3` wrappers (`fftw_r2c`, `fftw_c2r`, `fftw_c2c`, `mvfftw_r2c`, `mvfftw_c2c`, `mvfftw_c2r`, and their `2D`/`3D` variants) with memory bugs fixed; these are thin bindings for advanced users requiring maximum throughput
* Fixed `design_filter_fir` band-pass scaling reference frequency to use the average of the transition-band midpoints (matching `MATLAB` `fir1`/`scale_filter`), and fixed band-stop scaling to always normalize at `DC`; previously both could produce incorrect gain for asymmetric transition bands
* `filtfilt` now accepts a `Sos` (second-order sections) object from `gsignal` and delegates to `gsignal::filtfilt` in that case
* Added `catmull_rom_3d` for smooth `Catmull-Rom` spline interpolation through 3D point sequences, including closest-point projection from an arbitrary 3D point onto the curve
* Fixed `qfac` handling in the matrix to `quaternion` conversion so that right-handed transforms (negative determinant) are correctly represented
* `fir1` filter computation is significantly faster via an optimized internal implementation

# ravetools 0.2.4

* Added a naive implementation of non-negative matrix factorization in pure R
* Added finding and aligning stimulation pulses (for continuous stimulation)
* Added `find_peaks` to provide finding peaks along a trace of signal

# ravetools 0.2.3

* Added `gammatone_fast` filters to obtain the audio envelope at different frequencies
* Added `vcg_subset_certex` to subset mesh
* Added `vcg_subdivision` to up-sample mesh
* Added plane-generating function to create mesh for plane

# ravetools 0.2.2

* Fixed `clang20` warning and removed problematic `vcglib` code that use pointers after free.

# ravetools 0.2.1

* Using `std::nearbyint` instead `std::round` to round numbers to comply to `IEC-60559` standard that half numbers round to nearest even integers
* Implemented low-level `3D` volume sampling in `C++`: `resample_3d_volume`
* Exported `gsignal::resample` (however, some edits might be needed in the future to produce the same results as `Matlab`)
* Supported converting `ieegio` surface geometries to `mesh3d` to be used in `VCG` related functions
* Added `K-D` tree search to find closest points among two point clouds

# ravetools 0.2.0

* Added `vcg_raycaster` to find intersection of rays and given mesh object.

# ravetools 0.1.9

* Hot fix by adding `vctrs` to `Suggests` to fix the issue that fails the unit test due to an update in `testthat` package. 

# ravetools 0.1.8

* `design_filter` ensures frequency window cuts off within 0 to `Nyquist`
* `fftfilt` allows matrix input
* `filtfilt` with `a=1` (`FIR` filter) calls `fftfilt` to speed up
* Fixed `pwelch` incorrect power calculation. The results agrees mostly with `Matlab` function
* Using `hamming` window as the default in `pwelch`; exported additional window options such as `blackman` families, `bohmanwin`, `flattopwin`, and `hanning`

# ravetools 0.1.7

* Fixed a `c++` template issue caused by `vcglib`, which fails to compile under `clang19`
* `decimate` agrees with `Matlab` now (both `fir` and `iir` filters)
* Implemented `freqz2` to obtain the frequency response of a digital filter (similar to `gsignal::freqz` but with more accurate cutoff frequency calculation and more customize plots) 
* Added filter diagnostic plot `diagnose_filter`
* Added `check_filter` to obtain expected magnitude at given frequency and reciprocal condition number
* Removed dependency `signal` and use `gsignal` instead
* Added `design_filter` for both `fir` and `iir` filters, allowing both entry and intermediate users to design band-pass/stop, low/high pass filters easily. 
* The `iir` filter order generated from `design_filter` will be checked against `rcond_filter_ar` (reciprocal condition number) to make improve the numeric stability
* Fixed (hopefully) memory issues caused by `dijkstra`; the package passed `asan`, `valgrind` test provided by `rhub2`


# ravetools 0.1.6

* Fixed `dijkstra` method occasionally causing memory error. New method is much faster now.
* Fixed `plot.pwelch` not displaying the signal names correctly.
* Removed `c++17` requirement and supports 11, 14, and 17 standards

# ravetools 0.1.5

* Fixed `FIR1` filter
* Fixed `pwelch` throwing warnings when signal is zero (or zero power)
* Updated authorship

# ravetools 0.1.4

* Implemented `dijkstra` to find shortest paths in mesh
* Migrated and incorporated `vcglib`
* Fixed a `C++` template issue via type explicit calls
* Fixed `fir1` filter when band-passing signals with incorrect `n`
* `pwelch` plot works with zero power now; `mv_pwelch` plot error fixed

# ravetools 0.1.3

* Rewrote `band_pass2` to avoid `NA` generated when upper band frequency is `Nyquist`
* Added `Vector3`, `Matrix4`, `Quaternion` for in-place calculation
* Added support for `WASM`
* Fixed issues reported by `CRAN`: "format string is not a string literal (potentially insecure)"

# ravetools 0.1.2

* Compatible with the latest `filearray`
* Exported `grow_volume`
* `mesh_from_volume` no longer throw errors if the mesh does not form a manifold

# ravetools 0.1.1

* Fixed a precision issue that caused test failure on some machine

# ravetools 0.1.0

* Added `fill_surface` to fill in volume based on given surface mesh
* Added `mesh_from_volume` to generate mesh from volume. This function can be used together with `fill_surface` to generate surface envelope
* Added `register_volume` to align two imaging data using linear or non-linear registration
* Added `fftw` on `2D` image and `3D` volume data
* Added convolution for `1D`, `2D`, `3D` data using `FFT` 

# ravetools 0.0.9

* Fixed `pwelch` frequency not starting from zero issue
* Upgraded `TinyThread` using the latest pull-request to `RcppParallel`

# ravetools 0.0.8

* Added `interpolate_stimulation` to detect stimulation signals within the response and interpolate with smooth signals
* The package now imports `splines`
* Added `fast_quantile` and `fast_mvquantile` to improve the quantile/median calculation speed
* Fixed the `plot_signals` plotting range too large when signals have large values (such as stimulation)
* Fixed `TinyThreads` library memory leak issues
* Simplified `diagnose_channel`, avoid duplicated `pwelch` calculation

# ravetools 0.0.7

* Added signal `filter`, `filtfilt` that produce the same results as `Matlab` (with precision error)
* Added two ways to perform band-pass filters
* Allows multiple channels through `pwelch` as a row-major matrix to speed up calculation
* Added `wavelet_cycles_suggest` to provide default calculation of wavelet cycles
* Added internal argument `signature` to wavelet to resolve potential cache conflicts when running in multiple processes. (This allows `RAVE` to run wavelet on multiple subjects at the same time)

# ravetools 0.0.6

* Added decibel average in `pwelch`
* Allowed `pwelch` sampling frequency to be greater than the signal length
* Adjusted parameters diagnostic plot and `pwelch` plot to properly handle graph text, margin, axis
* Added `plot_signals` to plot multiple functional data within the same canvas

# ravetools 0.0.5

* Exposed `C++` code to `inst/includes` so other users can dynamically link to the functions (https://github.com/dipterix/ravetools/issues/5)
* Removed confusing in-place arguments in the `fftw` related code
* Corrected `fftw` plans to respect the flags
* Added `C++` to convert raw binary bytes to `uint`, `int`, `float`, and `string`

# ravetools 0.0.4

Parallel processes might use different temporary directory paths. To improve the performance, it is recommended to set a shared temporary directory, hence this version

* Allows temporary directories to be set via environment variable `RAVETOOLS_TEMPDIR` or option `ravetools.tempdir`. 

# ravetools 0.0.3

This version fixes a memory issue reported by `CRAN` check (`gcc-UBSAN`). 

* There is a potential integer overflow where `NA_INTEGER` is subtracted by one before being converted to `R_xlen_t` type. This update fixes this issue
* Removed `RcppParallel` and copied part of it into `inst/include` folder, with `TBB` removed under the `GPL-3` license framework.


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
