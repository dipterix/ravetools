# Native 3D volume registration (`'rigid'`, `'affine'`, or `'SyN'`)

Self-contained image registration for 3D volumes, implemented purely in
RcppEigen (no other external registration library). It mirrors the core
behavior of 'ANTs' `antsRegistration`: a multi-resolution,
physical-shift scaled gradient-descent optimizer driving a similarity
metric, working entirely in the anatomical `RAS`
(right-anterior-superior) space. Each volume carries its own `vox2ras`
(0-indexed voxel index to anatomical `RAS`) \\4\times 4\\ transform, so
volumes with different sampling, orientation, or field of view are
aligned correctly.

## Usage

``` r
register_volume3d(
  source,
  target,
  source_vox2ras = NULL,
  target_vox2ras = NULL,
  source_mask = NULL,
  target_mask = NULL,
  source_points = NULL,
  target_points = NULL,
  points_weight = 0.5,
  weights = NULL,
  type = c("rigid", "affine", "syn", "syn_only"),
  metric = "mattes",
  shrink_factors = c(4, 2, 1),
  smoothing_sigmas = c(2, 1, 0),
  iterations = c(1000, 500, 250),
  sampling_rate = 0.2,
  interpolation = "trilinear",
  number_of_bins = 32L,
  seed = 1L,
  init_transform = NULL,
  syn_iterations = c(40, 20, 0),
  syn_sigma = 3,
  verbose = TRUE
)
```

## Arguments

- source:

  the moving volume to be aligned, a 3D array (for example a `'CT'`);
  the result transform maps `target` into this image's space. May also
  be a `list` of already co-registered 3D arrays (multiple modalities
  such as `'T1'`, `'T2'`, an atlas constraint); all source channels must
  share the same dimensions and `vox2ras`

- target:

  the fixed/reference volume to align to, a 3D array (for example a
  `'MRI'`), or a `list` of co-registered arrays matching `source`
  channel-for-channel

- source_vox2ras, target_vox2ras:

  4x4 (or 3x4) matrices mapping the 0-indexed voxel coordinate
  (column-row-slice, `'C'`-style starting from 0, complying with
  `'NIfTI'`) to the `RAS` coordinate system; if `NULL`, the function
  looks for a `"vox2ras"` attribute on the array

- source_mask, target_mask:

  optional 3D mask arrays restricting where the metric is evaluated;
  default `NULL` (no mask, evaluate everywhere). A `target_mask` (on the
  target/fixed grid, same dimensions as `target`) limits which target
  voxels drive the registration; a `source_mask` (on the source/moving
  grid) drops samples that map outside it. Non-zero voxels are included.
  Besides focusing the alignment on a region of interest, masks speed
  things up by skipping background: the linear stage samples only inside
  the mask, and the deformable stage skips warping voxels outside the
  (dilated) mask. One mask per grid, shared across channels

- source_points, target_points:

  optional `N x 3` matrices of corresponding landmark coordinates that
  add a surface/landmark term to the **deformable**
  (`"syn"`/`"syn_only"`) stage; default `NULL` (no term). Row `i` of
  `target_points` (in the target/fixed `RAS`) and row `i` of
  `source_points` (in the source/moving `RAS`) must be the *same*
  anatomical location, for example corresponding cortical-surface
  vertices from a `'FreeSurfer'` spherical registration. The term pulls
  the warp of each target point onto its source correspondent,
  recovering cortical folding (`gyrification`) that the intensity metric
  blurs over, while the image metric still drives deep-brain
  `subcortical` structures. Must be supplied together with equal row
  counts. **Points must be in the same `RAS` frame as
  `source_vox2ras`/`target_vox2ras`** (note `'FreeSurfer'` surfaces use
  `surface/tkr` `RAS`, which differs from scanner `RAS` by `c_ras`)

- points_weight:

  relative weight of the landmark term against the image metric in the
  deformable stage; default `0.5`. Larger values follow the landmarks
  more closely. The sparse landmark force is attenuated by `syn_sigma`
  smoothing, so this typically needs tuning for a given point count and
  spacing

- weights:

  optional numeric weights, one per source/target pair, controlling each
  channel's contribution to the deformable cost; default is equal
  weighting. Weights are normalized internally to sum to 1. Only the
  deformable (`"syn"`/`"syn_only"`) stage is multivariate; the linear
  stage always uses the first (primary) pair

- type:

  type of transform to estimate; one of `'rigid'` (6 degrees of
  freedom), `'affine'` (12), `"syn"` (`affine` followed by a `SDR` -
  symmetric `diffeomorphic` deformation), or `'syn_only'` (`deformable`
  stage only — no `affine` is estimated; `init_transform` is used
  directly as the starting `affine`, useful when you already have a good
  linear alignment)

- metric:

  similarity metric: `"mattes"` (Mattes mutual information, the default,
  best for cross-modal such as `'CT'`-`'MRI'`), `"cc"` (normalized
  cross-correlation, for same-modality), or `"meansquares"` (mean
  squared intensity difference). With multiple channels, supply either a
  single metric (used for every pair) or a vector with one metric per
  source/target pair

- shrink_factors:

  integer down-sampling factors, one per resolution level (coarsest
  first); default `c(4, 2, 1)`

- smoothing_sigmas:

  Gaussian smoothing applied at each level, in voxels, same length as
  `shrink_factors`; default `c(2, 1, 0)`

- iterations:

  maximum optimizer iterations per level; default `c(1000, 500, 250)`

- sampling_rate:

  fraction of fixed voxels sampled to evaluate the metric (speeds up
  large volumes); default `0.2`

- interpolation:

  output interpolation used when warping each modality onto the target
  grid: `'trilinear'` (default), `'nearest'` (keeps label or
  segmentation values intact), or `"bspline"` (cubic `Catmull-Rom`,
  higher quality, slower). Like `metric`, supply a single value (applied
  to every channel) or one per source/target pair. This affects only the
  returned warped images, never the optimization's internal sampling
  (always `'trilinear'`, so every channel still produces smooth
  gradients)

- number_of_bins:

  number of histogram bins for the `"mattes"` metric; default `32`

- seed:

  random seed for the voxel sampler, for reproducibility

- init_transform:

  optional 4x4 initial `RAS`-to-`RAS` transform (fixed to moving) to
  start from

- syn_iterations, syn_sigma:

  deformable stage controls (only used when `type = "syn"` or
  `"syn_only"`): per-level iteration counts and the Gaussian
  regularization sigma (in voxels) applied to the update field. Both are
  recycled to the number of levels, so `syn_sigma` may be a vector to
  vary the regularization per stage (e.g. `c(3, 3, 1)` to relax it at
  the finest level for sharper detail)

- verbose:

  logical; if `TRUE` (default) print per-level and per-iteration
  progress to the console, including the current stage, shrink factor,
  smoothing sigma, cost metric, and step size (linear) or maximum
  displacement (deformable); useful to monitor convergence on large
  volumes

## Value

A list with:

- `transform`:

  the estimated 4x4 `RAS`-to-`RAS` linear transform mapping `target`
  (fixed) coordinates to `source` (moving) coordinates

- `image`:

  the (primary) `source` resampled onto the `target` grid

- `images`:

  a list with every source channel resampled onto the `target` grid
  using its own `interpolation`; `image` is the first element. For
  single-image input this is a length-1 list

- `forward_field`,`inverse_field`:

  (only for `"syn"`) the deformation fields

- `metric_trace`:

  the metric value across optimizer iterations

- `type`,`metric`:

  echoes of the inputs

## See also

[`apply_transform3d`](https://dipterix.org/ravetools/reference/apply_transform3d.md),
[`resample_3d_volume`](https://dipterix.org/ravetools/reference/resample_3d_volume.md)

## Examples

``` r
# \donttest{

# --- synthetic same-modality example -------------------------------------
nd <- c(50, 50, 50)
vox2ras <- diag(4); vox2ras[1:3, 4] <- -25
blob <- function(cx, cy, cz, s = 6) {
  g <- expand.grid(x = 0:(nd[1]-1), y = 0:(nd[2]-1), z = 0:(nd[3]-1))
  array(exp(-((g$x-cx)^2 + (g$y-cy)^2 + (g$z-cz)^2) / (2*s^2)), nd)
}
target <- blob(25, 25, 25)
source <- blob(28, 23, 26)            # shifted by a known (3, -2, 1) mm

res <- register_volume3d(
  source, target,
  source_vox2ras = vox2ras, target_vox2ras = vox2ras,
  type = "rigid", metric = "cc"
)
#> [rigid] level 1/3 (shrink=4, sigma=2.0): max 1000 iterations
#>   it     1  step = 4.000000: cost = -0.915936
#>   it     2  step = 2.000000: cost = -0.979488
#>   it     3  step = 2.000000: cost = -0.999069
#>   it     4  step = 1.000000: cost = -0.996996
#>   it     5  step = 0.500000: cost = -0.999831
#>   it     6  step = 0.500000: cost = -0.999220
#>   it     7  step = 0.250000: cost = -0.999940
#>   it     8  step = 0.250000: cost = -0.999846
#>   it     9  step = 0.125000: cost = -0.999994
#>   it    10  step = 0.125000: cost = -0.999943
#>   it    11  step = 0.062500: cost = -0.999993
#>   it    12  step = 0.062500: cost = -0.999994
#>   it    13  step = 0.031250: cost = -1.000000
#>   it    14  step = 0.015625: cost = -0.999999
#>   it    15  step = 0.007812: cost = -1.000000
#>   it    16  step = 0.007812: cost = -1.000000
#>   => converged (step < minStep)  best cost=-1.000000
#> [rigid] level 2/3 (shrink=2, sigma=1.0): max 500 iterations
#>   it     1  step = 2.000000: cost = -0.966638
#>   it     2  step = 1.000000: cost = -0.991557
#>   it     3  step = 1.000000: cost = -0.999982
#>   it     4  step = 0.500000: cost = -0.998462
#>   it     5  step = 0.250000: cost = -0.999695
#>   it     6  step = 0.250000: cost = -0.999981
#>   it     7  step = 0.125000: cost = -0.999955
#>   it     8  step = 0.062500: cost = -0.999998
#>   it     9  step = 0.062500: cost = -0.999982
#>   it    10  step = 0.031250: cost = -0.999998
#>   it    11  step = 0.031250: cost = -0.999998
#>   it    12  step = 0.015625: cost = -1.000000
#>   it    13  step = 0.007812: cost = -1.000000
#>   it    14  step = 0.003906: cost = -1.000000
#>   it    15  step = 0.003906: cost = -1.000000
#>   => converged (step < minStep)  best cost=-1.000000
#> [rigid] level 3/3 (shrink=1, sigma=0.0): max 250 iterations
#>   it     1  step = 1.000000: cost = -0.992627
#>   it     2  step = 0.500000: cost = -0.998147
#>   it     3  step = 0.500000: cost = -1.000000
#>   it     4  step = 0.250000: cost = -0.999536
#>   it     5  step = 0.125000: cost = -0.999890
#>   it     6  step = 0.125000: cost = -1.000000
#>   it     7  step = 0.062500: cost = -0.999975
#>   it     8  step = 0.031250: cost = -0.999995
#>   it     9  step = 0.031250: cost = -1.000000
#>   it    10  step = 0.015625: cost = -0.999999
#>   it    11  step = 0.007812: cost = -1.000000
#>   it    12  step = 0.007812: cost = -1.000000
#>   it    13  step = 0.003906: cost = -1.000000
#>   it    14  step = 0.003906: cost = -1.000000
#>   it    15  step = 0.001953: cost = -1.000000
#>   => converged (step < minStep)  best cost=-1.000000
res$transform[1:3, 4]                 # ~ c(3, -2, 1)
#> [1]  3.0000832 -2.0002276  0.9998873


# --- multimodal registration (several co-registered channels) ------------
# Two aligned modalities (e.g. a T1 and a T2) jointly drive the deformable
# stage. Each pair may use its own metric, and `weights` set their relative
# influence (normalized internally to sum to 1). All source channels must
# share a grid; likewise all target channels.
t1_target <- blob(25, 25, 25, s = 6)
t2_target <- blob(25, 25, 25, s = 9)  # same anatomy, different contrast
t1_source <- blob(27, 24, 26, s = 6)  # both channels share the same warp
t2_source <- blob(27, 24, 26, s = 9)

res_mm <- register_volume3d(
  source = list(t1_source, t2_source),
  target = list(t1_target, t2_target),
  source_vox2ras = vox2ras, target_vox2ras = vox2ras,
  weights = c(2, 1),                  # T1 counts twice as much as T2
  metric = c("mattes", "cc"),         # one metric per channel
  type = "syn"
)
#> [rigid] level 1/3 (shrink=4, sigma=2.0): max 1000 iterations
#>   it     1  step = 4.000000: cost = -0.473960
#>   it     2  step = 2.000000: cost = -0.529658
#>   it     3  step = 2.000000: cost = -0.589115
#>   it     4  step = 2.000000: cost = -0.559668
#>   it     5  step = 1.000000: cost = -0.564823
#>   it     6  step = 1.000000: cost = -0.587755
#>   it     7  step = 0.500000: cost = -0.584957
#>   it     8  step = 0.250000: cost = -0.590449
#>   it     9  step = 0.250000: cost = -0.587515
#>   it    10  step = 0.125000: cost = -0.590149
#>   it    11  step = 0.125000: cost = -0.590498
#>   it    12  step = 0.062500: cost = -0.590606
#>   it    13  step = 0.031250: cost = -0.590649
#>   it    14  step = 0.015625: cost = -0.590657
#>   it    15  step = 0.007812: cost = -0.590661
#>   => converged (step < minStep)  best cost=-0.590661
#> [rigid] level 2/3 (shrink=2, sigma=1.0): max 500 iterations
#>   it     1  step = 2.000000: cost = -0.441226
#>   it     2  step = 1.000000: cost = -0.508360
#>   it     3  step = 1.000000: cost = -0.570802
#>   it     4  step = 0.500000: cost = -0.555848
#>   it     5  step = 0.250000: cost = -0.570353
#>   it     6  step = 0.250000: cost = -0.570802
#>   it     7  step = 0.125000: cost = -0.572873
#>   it     8  step = 0.062500: cost = -0.572469
#>   it     9  step = 0.031250: cost = -0.572822
#>   it    10  step = 0.031250: cost = -0.572873
#>   it    11  step = 0.015625: cost = -0.572886
#>   it    12  step = 0.007812: cost = -0.572889
#>   it    13  step = 0.003906: cost = -0.572890
#>   => converged (step < minStep)  best cost=-0.572890
#> [rigid] level 3/3 (shrink=1, sigma=0.0): max 250 iterations
#>   it     1  step = 1.000000: cost = -0.470544
#>   it     2  step = 0.500000: cost = -0.508023
#>   it     3  step = 0.500000: cost = -0.538966
#>   it     4  step = 0.250000: cost = -0.523130
#>   it     5  step = 0.125000: cost = -0.529009
#>   it     6  step = 0.125000: cost = -0.538154
#>   it     7  step = 0.062500: cost = -0.539603
#>   it     8  step = 0.031250: cost = -0.540068
#>   it     9  step = 0.031250: cost = -0.537787
#>   it    10  step = 0.015625: cost = -0.539982
#>   it    11  step = 0.015625: cost = -0.535920
#>   it    12  step = 0.007812: cost = -0.537438
#>   it    13  step = 0.007812: cost = -0.540149
#>   it    14  step = 0.003906: cost = -0.538544
#>   it    15  step = 0.001953: cost = -0.539447
#>   => converged (step < minStep)  best cost=-0.540149
#> [affine] level 1/3 (shrink=4, sigma=2.0): max 1000 iterations
#>   it     1  step = 4.000000: cost = -0.496476
#>   it     2  step = 2.000000: cost = -0.539590
#>   it     3  step = 2.000000: cost = -0.522749
#>   it     4  step = 1.000000: cost = -0.569244
#>   it     5  step = 1.000000: cost = -0.558248
#>   it     6  step = 0.500000: cost = -0.577268
#>   it     7  step = 0.500000: cost = -0.575533
#>   it     8  step = 0.250000: cost = -0.583549
#>   it     9  step = 0.250000: cost = -0.583164
#>   it    10  step = 0.125000: cost = -0.585193
#>   it    11  step = 0.125000: cost = -0.585787
#>   it    12  step = 0.125000: cost = -0.586112
#>   it    13  step = 0.062500: cost = -0.586543
#>   it    14  step = 0.062500: cost = -0.586806
#>   it    15  step = 0.062500: cost = -0.586986
#>   it    16  step = 0.031250: cost = -0.587155
#>   it    17  step = 0.031250: cost = -0.587299
#>   it    18  step = 0.031250: cost = -0.587438
#>   it    19  step = 0.031250: cost = -0.587574
#>   it    20  step = 0.031250: cost = -0.587707
#>   it    21  step = 0.031250: cost = -0.587837
#>   it    22  step = 0.031250: cost = -0.587965
#>   it    23  step = 0.031250: cost = -0.588092
#>   it    24  step = 0.031250: cost = -0.588219
#>   it    25  step = 0.031250: cost = -0.588346
#>   it    26  step = 0.031250: cost = -0.588474
#>   it    27  step = 0.031250: cost = -0.588605
#>   it    28  step = 0.031250: cost = -0.588739
#>   it    29  step = 0.031250: cost = -0.588878
#>   it    30  step = 0.031250: cost = -0.589022
#>   it    31  step = 0.031250: cost = -0.589170
#>   it    32  step = 0.031250: cost = -0.589324
#>   it    33  step = 0.031250: cost = -0.589482
#>   it    34  step = 0.031250: cost = -0.589644
#>   it    35  step = 0.031250: cost = -0.589809
#>   it    36  step = 0.031250: cost = -0.589974
#>   it    37  step = 0.031250: cost = -0.590143
#>   it    38  step = 0.031250: cost = -0.590313
#>   it    39  step = 0.031250: cost = -0.590483
#>   it    40  step = 0.031250: cost = -0.590650
#>   it    41  step = 0.031250: cost = -0.590823
#>   it    42  step = 0.031250: cost = -0.590993
#>   it    43  step = 0.031250: cost = -0.591163
#>   it    44  step = 0.031250: cost = -0.591332
#>   it    45  step = 0.031250: cost = -0.591497
#>   it    46  step = 0.031250: cost = -0.591655
#>   it    47  step = 0.031250: cost = -0.591814
#>   it    48  step = 0.031250: cost = -0.591967
#>   it    49  step = 0.031250: cost = -0.592111
#>   it    50  step = 0.031250: cost = -0.592243
#>   it    51  step = 0.031250: cost = -0.592364
#>   it    52  step = 0.031250: cost = -0.592465
#>   it    53  step = 0.031250: cost = -0.592550
#>   it    54  step = 0.031250: cost = -0.592621
#>   it    55  step = 0.031250: cost = -0.592677
#>   it    56  step = 0.015625: cost = -0.592727
#>   it    57  step = 0.015625: cost = -0.592767
#>   it    58  step = 0.015625: cost = -0.592805
#>   it    59  step = 0.015625: cost = -0.592833
#>   it    60  step = 0.015625: cost = -0.592856
#>   it    61  step = 0.015625: cost = -0.592873
#>   it    62  step = 0.015625: cost = -0.592884
#>   it    63  step = 0.007812: cost = -0.592897
#>   it    64  step = 0.007812: cost = -0.592907
#>   it    65  step = 0.007812: cost = -0.592917
#>   it    66  step = 0.007812: cost = -0.592927
#>   it    67  step = 0.007812: cost = -0.592937
#>   it    68  step = 0.007812: cost = -0.592946
#>   it    69  step = 0.007812: cost = -0.592954
#>   it    70  step = 0.007812: cost = -0.592963
#>   it    71  step = 0.007812: cost = -0.592971
#>   it    72  step = 0.007812: cost = -0.592979
#>   it    73  step = 0.007812: cost = -0.592987
#>   it    74  step = 0.007812: cost = -0.592995
#>   it    75  step = 0.007812: cost = -0.593003
#>   it    76  step = 0.007812: cost = -0.593011
#>   it    77  step = 0.007812: cost = -0.593019
#>   it    78  step = 0.007812: cost = -0.593027
#>   it    79  step = 0.007812: cost = -0.593034
#>   it    80  step = 0.007812: cost = -0.593041
#>   it    81  step = 0.007812: cost = -0.593049
#>   it    82  step = 0.007812: cost = -0.593057
#>   it    83  step = 0.007812: cost = -0.593066
#>   it    84  step = 0.007812: cost = -0.593075
#>   it    85  step = 0.007812: cost = -0.593085
#>   it    86  step = 0.007812: cost = -0.593096
#>   it    87  step = 0.007812: cost = -0.593107
#>   it    88  step = 0.007812: cost = -0.593118
#>   it    89  step = 0.007812: cost = -0.593129
#>   it    90  step = 0.007812: cost = -0.593141
#>   it    91  step = 0.007812: cost = -0.593153
#>   it    92  step = 0.007812: cost = -0.593166
#>   it    93  step = 0.007812: cost = -0.593183
#>   it    94  step = 0.007812: cost = -0.593202
#>   it    95  step = 0.007812: cost = -0.593220
#>   it    96  step = 0.007812: cost = -0.593239
#>   it    97  step = 0.007812: cost = -0.593258
#>   it    98  step = 0.007812: cost = -0.593277
#>   it    99  step = 0.007812: cost = -0.593296
#>   it   100  step = 0.007812: cost = -0.593316
#>   it   101  step = 0.007812: cost = -0.593337
#>   it   102  step = 0.007812: cost = -0.593358
#>   it   103  step = 0.007812: cost = -0.593379
#>   it   104  step = 0.007812: cost = -0.593401
#>   it   105  step = 0.007812: cost = -0.593423
#>   it   106  step = 0.007812: cost = -0.593446
#>   it   107  step = 0.007812: cost = -0.593469
#>   it   108  step = 0.007812: cost = -0.593493
#>   it   109  step = 0.007812: cost = -0.593517
#>   it   110  step = 0.007812: cost = -0.593541
#>   it   111  step = 0.007812: cost = -0.593566
#>   it   112  step = 0.007812: cost = -0.593591
#>   it   113  step = 0.007812: cost = -0.593616
#>   it   114  step = 0.007812: cost = -0.593642
#>   it   115  step = 0.007812: cost = -0.593668
#>   it   116  step = 0.007812: cost = -0.593694
#>   it   117  step = 0.007812: cost = -0.593721
#>   it   118  step = 0.007812: cost = -0.593748
#>   it   119  step = 0.007812: cost = -0.593775
#>   it   120  step = 0.007812: cost = -0.593802
#>   it   121  step = 0.007812: cost = -0.593829
#>   it   122  step = 0.007812: cost = -0.593857
#>   it   123  step = 0.007812: cost = -0.593884
#>   it   124  step = 0.007812: cost = -0.593912
#>   it   125  step = 0.007812: cost = -0.593939
#>   it   126  step = 0.007812: cost = -0.593967
#>   it   127  step = 0.007812: cost = -0.593994
#>   it   128  step = 0.007812: cost = -0.594022
#>   it   129  step = 0.007812: cost = -0.594049
#>   it   130  step = 0.007812: cost = -0.594076
#>   it   131  step = 0.007812: cost = -0.594103
#>   it   132  step = 0.007812: cost = -0.594130
#>   it   133  step = 0.007812: cost = -0.594157
#>   it   134  step = 0.007812: cost = -0.594184
#>   it   135  step = 0.007812: cost = -0.594210
#>   it   136  step = 0.007812: cost = -0.594237
#>   it   137  step = 0.007812: cost = -0.594263
#>   it   138  step = 0.007812: cost = -0.594289
#>   it   139  step = 0.007812: cost = -0.594314
#>   it   140  step = 0.007812: cost = -0.594340
#>   it   141  step = 0.007812: cost = -0.594365
#>   it   142  step = 0.007812: cost = -0.594390
#>   it   143  step = 0.007812: cost = -0.594414
#>   it   144  step = 0.007812: cost = -0.594438
#>   it   145  step = 0.007812: cost = -0.594462
#>   it   146  step = 0.007812: cost = -0.594486
#>   it   147  step = 0.007812: cost = -0.594509
#>   it   148  step = 0.007812: cost = -0.594531
#>   it   149  step = 0.007812: cost = -0.594554
#>   it   150  step = 0.007812: cost = -0.594576
#>   it   151  step = 0.007812: cost = -0.594597
#>   it   152  step = 0.007812: cost = -0.594618
#>   it   153  step = 0.007812: cost = -0.594639
#>   it   154  step = 0.007812: cost = -0.594659
#>   it   155  step = 0.007812: cost = -0.594679
#>   it   156  step = 0.007812: cost = -0.594698
#>   it   157  step = 0.007812: cost = -0.594717
#>   it   158  step = 0.007812: cost = -0.594735
#>   it   159  step = 0.007812: cost = -0.594753
#>   it   160  step = 0.007812: cost = -0.594771
#>   it   161  step = 0.007812: cost = -0.594788
#>   it   162  step = 0.007812: cost = -0.594804
#>   it   163  step = 0.007812: cost = -0.594820
#>   it   164  step = 0.007812: cost = -0.594836
#>   it   165  step = 0.007812: cost = -0.594851
#>   it   166  step = 0.007812: cost = -0.594866
#>   it   167  step = 0.007812: cost = -0.594880
#>   it   168  step = 0.007812: cost = -0.594893
#>   it   169  step = 0.007812: cost = -0.594907
#>   it   170  step = 0.007812: cost = -0.594920
#>   it   171  step = 0.007812: cost = -0.594932
#>   it   172  step = 0.007812: cost = -0.594944
#>   it   173  step = 0.007812: cost = -0.594956
#>   it   174  step = 0.007812: cost = -0.594967
#>   it   175  step = 0.007812: cost = -0.594978
#>   it   176  step = 0.007812: cost = -0.594989
#>   it   177  step = 0.007812: cost = -0.594999
#>   it   178  step = 0.007812: cost = -0.595008
#>   it   179  step = 0.007812: cost = -0.595018
#>   it   180  step = 0.007812: cost = -0.595027
#>   it   181  step = 0.007812: cost = -0.595035
#>   it   182  step = 0.007812: cost = -0.595044
#>   it   183  step = 0.007812: cost = -0.595052
#>   it   184  step = 0.007812: cost = -0.595060
#>   it   185  step = 0.007812: cost = -0.595067
#>   it   186  step = 0.007812: cost = -0.595075
#>   it   187  step = 0.007812: cost = -0.595082
#>   it   188  step = 0.007812: cost = -0.595089
#>   it   189  step = 0.007812: cost = -0.595096
#>   it   190  step = 0.007812: cost = -0.595103
#>   it   191  step = 0.007812: cost = -0.595109
#>   it   192  step = 0.007812: cost = -0.595116
#>   it   193  step = 0.007812: cost = -0.595122
#>   it   194  step = 0.007812: cost = -0.595128
#>   it   195  step = 0.007812: cost = -0.595134
#>   it   196  step = 0.007812: cost = -0.595140
#>   it   197  step = 0.007812: cost = -0.595144
#>   it   198  step = 0.007812: cost = -0.595148
#>   => converged (step < minStep)  best cost=-0.595148
#> [affine] level 2/3 (shrink=2, sigma=1.0): max 500 iterations
#>   it     1  step = 2.000000: cost = -0.436550
#>   it     2  step = 1.000000: cost = -0.520698
#>   it     3  step = 1.000000: cost = -0.582069
#>   it     4  step = 0.500000: cost = -0.582235
#>   it     5  step = 0.250000: cost = -0.595145
#>   it     6  step = 0.250000: cost = -0.593523
#>   it     7  step = 0.125000: cost = -0.595179
#>   it     8  step = 0.062500: cost = -0.595617
#>   it     9  step = 0.031250: cost = -0.595766
#>   it    10  step = 0.015625: cost = -0.595816
#>   it    11  step = 0.015625: cost = -0.595842
#>   it    12  step = 0.015625: cost = -0.595843
#>   it    13  step = 0.007812: cost = -0.595863
#>   it    14  step = 0.007812: cost = -0.595872
#>   it    15  step = 0.007812: cost = -0.595879
#>   it    16  step = 0.007812: cost = -0.595884
#>   it    17  step = 0.003906: cost = -0.595890
#>   it    18  step = 0.003906: cost = -0.595895
#>   it    19  step = 0.003906: cost = -0.595899
#>   it    20  step = 0.003906: cost = -0.595904
#>   it    21  step = 0.003906: cost = -0.595909
#>   it    22  step = 0.003906: cost = -0.595914
#>   it    23  step = 0.003906: cost = -0.595919
#>   it    24  step = 0.003906: cost = -0.595924
#>   it    25  step = 0.003906: cost = -0.595929
#>   it    26  step = 0.003906: cost = -0.595934
#>   it    27  step = 0.003906: cost = -0.595940
#>   it    28  step = 0.003906: cost = -0.595945
#>   it    29  step = 0.003906: cost = -0.595951
#>   it    30  step = 0.003906: cost = -0.595957
#>   it    31  step = 0.003906: cost = -0.595962
#>   it    32  step = 0.003906: cost = -0.595968
#>   it    33  step = 0.003906: cost = -0.595973
#>   it    34  step = 0.003906: cost = -0.595978
#>   it    35  step = 0.003906: cost = -0.595984
#>   it    36  step = 0.003906: cost = -0.595989
#>   it    37  step = 0.003906: cost = -0.595995
#>   it    38  step = 0.003906: cost = -0.596000
#>   it    39  step = 0.003906: cost = -0.596005
#>   it    40  step = 0.003906: cost = -0.596010
#>   it    41  step = 0.003906: cost = -0.596016
#>   it    42  step = 0.003906: cost = -0.596021
#>   it    43  step = 0.003906: cost = -0.596026
#>   it    44  step = 0.003906: cost = -0.596031
#>   it    45  step = 0.003906: cost = -0.596036
#>   it    46  step = 0.003906: cost = -0.596041
#>   it    47  step = 0.003906: cost = -0.596046
#>   it    48  step = 0.003906: cost = -0.596051
#>   it    49  step = 0.003906: cost = -0.596056
#>   it    50  step = 0.003906: cost = -0.596061
#>   it    51  step = 0.003906: cost = -0.596066
#>   it    52  step = 0.003906: cost = -0.596071
#>   it    53  step = 0.003906: cost = -0.596076
#>   it    54  step = 0.003906: cost = -0.596080
#>   it    55  step = 0.003906: cost = -0.596085
#>   it    56  step = 0.003906: cost = -0.596090
#>   it    57  step = 0.003906: cost = -0.596095
#>   it    58  step = 0.003906: cost = -0.596099
#>   it    59  step = 0.003906: cost = -0.596104
#>   it    60  step = 0.003906: cost = -0.596108
#>   it    61  step = 0.003906: cost = -0.596113
#>   it    62  step = 0.003906: cost = -0.596117
#>   it    63  step = 0.003906: cost = -0.596122
#>   it    64  step = 0.003906: cost = -0.596126
#>   it    65  step = 0.003906: cost = -0.596131
#>   it    66  step = 0.003906: cost = -0.596135
#>   it    67  step = 0.003906: cost = -0.596140
#>   it    68  step = 0.003906: cost = -0.596144
#>   it    69  step = 0.003906: cost = -0.596149
#>   it    70  step = 0.003906: cost = -0.596153
#>   it    71  step = 0.003906: cost = -0.596157
#>   it    72  step = 0.003906: cost = -0.596162
#>   it    73  step = 0.003906: cost = -0.596166
#>   it    74  step = 0.003906: cost = -0.596171
#>   it    75  step = 0.003906: cost = -0.596175
#>   it    76  step = 0.003906: cost = -0.596179
#>   it    77  step = 0.003906: cost = -0.596184
#>   it    78  step = 0.003906: cost = -0.596188
#>   it    79  step = 0.003906: cost = -0.596192
#>   it    80  step = 0.003906: cost = -0.596196
#>   it    81  step = 0.003906: cost = -0.596201
#>   it    82  step = 0.003906: cost = -0.596205
#>   it    83  step = 0.003906: cost = -0.596209
#>   it    84  step = 0.003906: cost = -0.596213
#>   it    85  step = 0.003906: cost = -0.596217
#>   it    86  step = 0.003906: cost = -0.596221
#>   it    87  step = 0.003906: cost = -0.596226
#>   it    88  step = 0.003906: cost = -0.596230
#>   it    89  step = 0.003906: cost = -0.596234
#>   it    90  step = 0.003906: cost = -0.596238
#>   it    91  step = 0.003906: cost = -0.596242
#>   it    92  step = 0.003906: cost = -0.596245
#>   it    93  step = 0.003906: cost = -0.596249
#>   it    94  step = 0.003906: cost = -0.596253
#>   it    95  step = 0.003906: cost = -0.596257
#>   it    96  step = 0.003906: cost = -0.596261
#>   it    97  step = 0.003906: cost = -0.596264
#>   it    98  step = 0.003906: cost = -0.596268
#>   it    99  step = 0.003906: cost = -0.596271
#>   it   100  step = 0.003906: cost = -0.596275
#>   it   101  step = 0.003906: cost = -0.596278
#>   it   102  step = 0.003906: cost = -0.596281
#>   it   103  step = 0.003906: cost = -0.596285
#>   it   104  step = 0.003906: cost = -0.596288
#>   it   105  step = 0.003906: cost = -0.596291
#>   it   106  step = 0.003906: cost = -0.596294
#>   it   107  step = 0.003906: cost = -0.596296
#>   it   108  step = 0.003906: cost = -0.596299
#>   it   109  step = 0.003906: cost = -0.596301
#>   => converged (step < minStep)  best cost=-0.596301
#> [affine] level 3/3 (shrink=1, sigma=0.0): max 250 iterations
#>   it     1  step = 1.000000: cost = -0.524449
#>   it     2  step = 0.500000: cost = -0.548692
#>   it     3  step = 0.250000: cost = -0.559561
#>   it     4  step = 0.250000: cost = -0.565430
#>   it     5  step = 0.250000: cost = -0.562228
#>   it     6  step = 0.125000: cost = -0.567481
#>   it     7  step = 0.125000: cost = -0.565607
#>   it     8  step = 0.062500: cost = -0.567467
#>   it     9  step = 0.062500: cost = -0.567487
#>   it    10  step = 0.031250: cost = -0.567719
#>   it    11  step = 0.015625: cost = -0.567724
#>   it    12  step = 0.007812: cost = -0.567732
#>   it    13  step = 0.003906: cost = -0.567735
#>   it    14  step = 0.003906: cost = -0.567737
#>   it    15  step = 0.003906: cost = -0.567740
#>   it    16  step = 0.003906: cost = -0.567743
#>   it    17  step = 0.003906: cost = -0.567746
#>   it    18  step = 0.003906: cost = -0.567749
#>   it    19  step = 0.003906: cost = -0.567752
#>   it    20  step = 0.003906: cost = -0.567754
#>   it    21  step = 0.003906: cost = -0.567757
#>   it    22  step = 0.003906: cost = -0.567760
#>   it    23  step = 0.003906: cost = -0.567763
#>   it    24  step = 0.003906: cost = -0.567765
#>   it    25  step = 0.003906: cost = -0.567768
#>   it    26  step = 0.003906: cost = -0.567771
#>   it    27  step = 0.003906: cost = -0.567773
#>   it    28  step = 0.003906: cost = -0.567776
#>   it    29  step = 0.003906: cost = -0.567779
#>   it    30  step = 0.003906: cost = -0.567781
#>   it    31  step = 0.003906: cost = -0.567784
#>   it    32  step = 0.003906: cost = -0.567786
#>   it    33  step = 0.003906: cost = -0.567789
#>   it    34  step = 0.003906: cost = -0.567791
#>   it    35  step = 0.003906: cost = -0.567794
#>   it    36  step = 0.003906: cost = -0.567796
#>   it    37  step = 0.003906: cost = -0.567798
#>   it    38  step = 0.003906: cost = -0.567801
#>   it    39  step = 0.003906: cost = -0.567803
#>   it    40  step = 0.003906: cost = -0.567805
#>   it    41  step = 0.003906: cost = -0.567808
#>   it    42  step = 0.003906: cost = -0.567810
#>   it    43  step = 0.003906: cost = -0.567812
#>   it    44  step = 0.003906: cost = -0.567815
#>   it    45  step = 0.003906: cost = -0.567817
#>   it    46  step = 0.003906: cost = -0.567819
#>   it    47  step = 0.003906: cost = -0.567822
#>   it    48  step = 0.003906: cost = -0.567824
#>   it    49  step = 0.003906: cost = -0.567826
#>   it    50  step = 0.003906: cost = -0.567828
#>   it    51  step = 0.003906: cost = -0.567830
#>   it    52  step = 0.003906: cost = -0.567833
#>   it    53  step = 0.003906: cost = -0.567835
#>   it    54  step = 0.003906: cost = -0.567837
#>   it    55  step = 0.003906: cost = -0.567839
#>   it    56  step = 0.003906: cost = -0.567841
#>   it    57  step = 0.003906: cost = -0.567843
#>   it    58  step = 0.003906: cost = -0.567845
#>   it    59  step = 0.001953: cost = -0.567846
#>   it    60  step = 0.001953: cost = -0.567848
#>   it    61  step = 0.001953: cost = -0.567849
#>   it    62  step = 0.001953: cost = -0.567850
#>   it    63  step = 0.001953: cost = -0.567852
#>   it    64  step = 0.001953: cost = -0.567853
#>   it    65  step = 0.001953: cost = -0.567854
#>   it    66  step = 0.001953: cost = -0.567856
#>   it    67  step = 0.001953: cost = -0.567857
#>   it    68  step = 0.001953: cost = -0.567858
#>   it    69  step = 0.001953: cost = -0.567860
#>   it    70  step = 0.001953: cost = -0.567861
#>   it    71  step = 0.001953: cost = -0.567862
#>   it    72  step = 0.001953: cost = -0.567864
#>   it    73  step = 0.001953: cost = -0.567865
#>   it    74  step = 0.001953: cost = -0.567866
#>   it    75  step = 0.001953: cost = -0.567867
#>   it    76  step = 0.001953: cost = -0.567869
#>   it    77  step = 0.001953: cost = -0.567870
#>   it    78  step = 0.001953: cost = -0.567871
#>   it    79  step = 0.001953: cost = -0.567872
#>   it    80  step = 0.001953: cost = -0.567873
#>   it    81  step = 0.001953: cost = -0.567874
#>   it    82  step = 0.001953: cost = -0.567876
#>   it    83  step = 0.001953: cost = -0.567877
#>   it    84  step = 0.001953: cost = -0.567878
#>   it    85  step = 0.001953: cost = -0.567879
#>   it    86  step = 0.001953: cost = -0.567880
#>   it    87  step = 0.001953: cost = -0.567881
#>   it    88  step = 0.001953: cost = -0.567882
#>   it    89  step = 0.001953: cost = -0.567883
#>   it    90  step = 0.001953: cost = -0.567884
#>   it    91  step = 0.001953: cost = -0.567885
#>   it    92  step = 0.001953: cost = -0.567886
#>   it    93  step = 0.001953: cost = -0.567887
#>   it    94  step = 0.001953: cost = -0.567888
#>   it    95  step = 0.001953: cost = -0.567889
#>   it    96  step = 0.001953: cost = -0.567890
#>   it    97  step = 0.001953: cost = -0.567891
#>   it    98  step = 0.001953: cost = -0.567892
#>   it    99  step = 0.001953: cost = -0.567893
#>   it   100  step = 0.001953: cost = -0.567894
#>   it   101  step = 0.001953: cost = -0.567895
#>   it   102  step = 0.001953: cost = -0.567896
#>   it   103  step = 0.001953: cost = -0.567897
#>   it   104  step = 0.001953: cost = -0.567898
#>   it   105  step = 0.001953: cost = -0.567898
#>   it   106  step = 0.001953: cost = -0.567899
#>   it   107  step = 0.001953: cost = -0.567900
#>   it   108  step = 0.001953: cost = -0.567901
#>   it   109  step = 0.001953: cost = -0.567902
#>   it   110  step = 0.001953: cost = -0.567903
#>   it   111  step = 0.001953: cost = -0.567904
#>   it   112  step = 0.001953: cost = -0.567905
#>   it   113  step = 0.001953: cost = -0.567906
#>   it   114  step = 0.001953: cost = -0.567906
#>   it   115  step = 0.001953: cost = -0.567907
#>   it   116  step = 0.001953: cost = -0.567908
#>   it   117  step = 0.001953: cost = -0.567909
#>   it   118  step = 0.001953: cost = -0.567909
#>   it   119  step = 0.001953: cost = -0.567910
#>   it   120  step = 0.001953: cost = -0.567911
#>   it   121  step = 0.001953: cost = -0.567911
#>   it   122  step = 0.001953: cost = -0.567912
#>   it   123  step = 0.001953: cost = -0.567913
#>   it   124  step = 0.001953: cost = -0.567914
#>   it   125  step = 0.001953: cost = -0.567914
#>   it   126  step = 0.001953: cost = -0.567915
#>   it   127  step = 0.001953: cost = -0.567916
#>   it   128  step = 0.001953: cost = -0.567916
#>   it   129  step = 0.001953: cost = -0.567917
#>   it   130  step = 0.001953: cost = -0.567918
#>   it   131  step = 0.001953: cost = -0.567918
#>   it   132  step = 0.001953: cost = -0.567919
#>   it   133  step = 0.001953: cost = -0.567919
#>   it   134  step = 0.001953: cost = -0.567920
#>   it   135  step = 0.001953: cost = -0.567921
#>   it   136  step = 0.001953: cost = -0.567921
#>   it   137  step = 0.001953: cost = -0.567922
#>   it   138  step = 0.001953: cost = -0.567922
#>   it   139  step = 0.001953: cost = -0.567923
#>   it   140  step = 0.001953: cost = -0.567923
#>   it   141  step = 0.001953: cost = -0.567924
#>   it   142  step = 0.001953: cost = -0.567925
#>   it   143  step = 0.001953: cost = -0.567925
#>   it   144  step = 0.001953: cost = -0.567926
#>   it   145  step = 0.001953: cost = -0.567926
#>   it   146  step = 0.001953: cost = -0.567927
#>   it   147  step = 0.001953: cost = -0.567927
#>   it   148  step = 0.001953: cost = -0.567928
#>   it   149  step = 0.001953: cost = -0.567928
#>   it   150  step = 0.001953: cost = -0.567929
#>   => converged (step < minStep)  best cost=-0.567929
#> [SyN] level 1/3 (shrink=4, sigma=2.0, flow=3.0, channels=2): max 40 iterations
#>   it     1: cost = -0.892844  field_max = 0.80 mm
#>   it     2: cost = -0.904248  field_max = 1.60 mm
#>   it     3: cost = -0.914074  field_max = 2.40 mm
#>   it     4: cost = -0.924056  field_max = 3.19 mm
#>   it     5: cost = -0.931777  field_max = 3.99 mm
#>   it     6: cost = -0.938930  field_max = 4.79 mm
#>   it     7: cost = -0.944740  field_max = 5.59 mm
#>   it     8: cost = -0.950270  field_max = 6.36 mm
#>   it     9: cost = -0.955244  field_max = 7.06 mm
#>   it    10: cost = -0.960497  field_max = 7.67 mm
#>   it    11: cost = -0.963909  field_max = 8.39 mm
#>   it    12: cost = -0.966920  field_max = 9.12 mm
#>   it    13: cost = -0.970232  field_max = 9.83 mm
#>   it    14: cost = -0.973304  field_max = 10.40 mm
#>   it    15: cost = -0.975462  field_max = 10.85 mm
#>   it    16: cost = -0.978081  field_max = 11.08 mm
#>   it    17: cost = -0.979111  field_max = 11.27 mm
#>   it    18: cost = -0.979603  field_max = 11.43 mm
#>   it    19: cost = -0.980324  field_max = 11.56 mm
#>   it    20: cost = -0.981028  field_max = 11.67 mm
#>   it    21: cost = -0.981958  field_max = 11.81 mm
#>   it    22: cost = -0.982852  field_max = 11.94 mm
#>   it    23: cost = -0.984064  field_max = 12.05 mm
#>   it    24: cost = -0.984800  field_max = 12.20 mm
#>   it    25: cost = -0.985892  field_max = 12.33 mm
#>   it    26: cost = -0.986956  field_max = 12.50 mm
#>   it    27: cost = -0.987708  field_max = 12.57 mm
#>   it    28: cost = -0.988188  field_max = 12.64 mm
#>   it    29: cost = -0.988802  field_max = 12.70 mm
#>   it    30: cost = -0.989299  field_max = 12.75 mm
#>   it    31: cost = -0.989959  field_max = 12.77 mm
#>   it    32: cost = -0.990740  field_max = 12.79 mm
#>   it    33: cost = -0.991660  field_max = 12.98 mm
#>   it    34: cost = -0.992293  field_max = 13.15 mm
#>   it    35: cost = -0.992826  field_max = 13.32 mm
#>   it    36: cost = -0.993394  field_max = 13.49 mm
#>   it    37: cost = -0.993905  field_max = 13.70 mm
#>   it    38: cost = -0.994207  field_max = 14.00 mm
#>   it    39: cost = -0.994230  field_max = 14.16 mm
#>   it    40: cost = -0.993726  field_max = 14.30 mm
#> [SyN] level 2/3 (shrink=2, sigma=1.0, flow=3.0, channels=2): max 20 iterations
#>   it     1: cost = -0.896596  field_max = 13.98 mm
#>   it     2: cost = -0.896885  field_max = 14.10 mm
#>   it     3: cost = -0.896447  field_max = 14.22 mm
#>   it     4: cost = -0.895603  field_max = 14.33 mm
#>   it     5: cost = -0.894657  field_max = 14.45 mm
#>   it     6: cost = -0.893685  field_max = 14.55 mm
#>   it     7: cost = -0.892752  field_max = 14.61 mm
#>   it     8: cost = -0.891746  field_max = 14.68 mm
#>   it     9: cost = -0.890929  field_max = 14.75 mm
#>   it    10: cost = -0.890349  field_max = 14.80 mm
#>   it    11: cost = -0.889983  field_max = 14.84 mm
#>   it    12: cost = -0.889733  field_max = 14.90 mm
#>   it    13: cost = -0.889552  field_max = 14.94 mm
#>   it    14: cost = -0.889323  field_max = 15.00 mm
#>   it    15: cost = -0.889113  field_max = 15.05 mm
#>   it    16: cost = -0.888916  field_max = 15.10 mm
#>   it    17: cost = -0.888720  field_max = 15.15 mm
#>   it    18: cost = -0.888473  field_max = 15.20 mm
#>   it    19: cost = -0.888188  field_max = 15.25 mm
#>   it    20: cost = -0.887909  field_max = 15.30 mm
#> [SyN] level 3/3 (shrink=1, sigma=0.0, flow=3.0, channels=2): max 0 iterations
res_mm$transform[1:3, 4]
#> [1]  2.0017331 -1.0002336  0.9990778

# }
```
