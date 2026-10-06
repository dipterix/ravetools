# Prior-guided tissue segmentation with a `Gaussian` mixture and an `MRF`

Segments a 3D intensity volume into `K` tissue classes with a finite
`Gaussian` mixture model, optional spatial prior probability maps, and a
mean-field `Markov` random field (`MRF`) that favors spatially coherent
labels. This is a self-contained re-implementation of the `Atropos`
algorithm of `'ANTs'` (`Avants` and colleagues, 2011) in Rcpp; the
defaults reproduce the `'ANTsPy'` `atropos` call `i = 'Kmeans[3]'`,
`m = '[0.2,1x1x1]'`, `c = '[5,0]'`, `priorweight = 0.25`.

## Usage

``` r
segment_volume_tissue_gmm(
  volume,
  mask = NULL,
  priors = NULL,
  n_classes = 3L,
  prior_weight = 0.25,
  mrf_beta = 0.2,
  mrf_radius = 1L,
  iterations = 5L,
  tolerance = 0,
  vox2ras = NULL,
  verbose = FALSE
)
```

## Arguments

- volume:

  a 3D numeric (or integer, logical) array of intensities; a `"vox2ras"`
  attribute, if present, is used when `vox2ras` is `NULL`. Trailing
  singleton dimensions (for example a `NIfTI` volume read as
  `nx x ny x nz x 1`) are dropped from `volume`, `mask` and `priors`,
  and the outputs are always 3D. Extreme outlier intensities inside the
  mask (hot voxels, raw `CT` values) should be truncated or clipped
  beforehand, since the `k-means` seeds are spread over the intensity
  range (see 'Initialization'). The result also depends on the absolute
  intensity scale, because the `Gaussian` density is floored at `1e-10`
  as in `'ANTs'` (see 'Details'): intensities in very large units (class
  standard deviations beyond roughly `1e7`) should be rescaled first

- mask:

  optional 3D array of the same dimensions as `volume`; non-zero (and
  non-`NA`) voxels are segmented, everything else receives label `0` and
  zero posteriors. The default `NULL` uses every voxel with a finite
  intensity. Voxels with non-finite intensities are always excluded from
  the mask

- priors:

  optional list of `K` non-negative 3D arrays on the grid of `volume`,
  the spatial prior probability of each class (for example tissue
  probability maps warped from a template with
  `apply_transform3d_volume`). They are normalized per voxel; at an
  isolated voxel where all priors are zero a uniform prior is used, but
  priors that are zero at every voxel inside the mask (for example
  integer-truncated probability maps, or maps warped onto the wrong
  grid) are an error, as is a prior that is zero everywhere inside the
  mask. The number of classes is the number of priors, and the output
  classes follow the order of the list. Values outside the mask are
  ignored

- n_classes:

  number of classes for the `k-means` initialization (default `3`, the
  classic `CSF`, gray matter, white matter); ignored when `priors` are
  given, unless it is passed explicitly and conflicts with the number of
  priors, which is an error

- prior_weight:

  weight \\w \in \[0, 1\]\\ of the spatial priors in the posterior
  (default `0.25`), with the meaning of the `priorweight` of `'ANTs'`;
  `0` uses the priors only to initialize the model, `1` ignores the
  image and segments by the priors alone (see 'Details'). Ignored
  without `priors`

- mrf_beta:

  `MRF` smoothing factor \\\beta \ge 0\\ (default `0.2`); larger values
  produce smoother label maps, `0` disables the `MRF`

- mrf_radius:

  `MRF` neighborhood radius in voxels, either one non-negative integer
  for all axes or one per axis (default `1`, the 26-connected
  neighborhood); non-integer values are an error and `0` disables the
  `MRF`. A radius beyond the grid extent is equivalent to the extent
  along that axis and is clamped to it

- iterations:

  maximum number of `EM` iterations (default `5`); the iterations
  usually stop earlier, see `tolerance`

- tolerance:

  non-negative convergence threshold on the mean maximum posterior
  probability over the mask (default `0`), with the meaning of the
  convergence threshold of `'ANTs'` `Atropos`: after each iteration
  (except the first) the change of this value from the previous
  iteration is compared with `tolerance`, and the iterations stop as
  soon as the change is below it. With the default `0` the iterations
  stop as soon as the value decreases (as `Atropos` does with
  `c = [n,0]`), otherwise they run up to `iterations`; the posteriors
  and labels of the last iteration performed are returned (see
  'Details'). Negative values are an error

- vox2ras:

  optional 4x4 (or 3x4) matrix mapping 0-indexed voxel coordinates to
  `RAS`; only its \\3\times 3\\ part matters, to express the `MRF`
  neighbor weights in physical distance. Default `NULL` looks for a
  `"vox2ras"` attribute on `volume` and otherwise uses voxel units (unit
  spacing)

- verbose:

  logical; print per-iteration progress (default `FALSE`)

## Value

A list with

- `segmentation`:

  integer 3D array, `0` outside the mask and `1..K` inside, the `argmax`
  of the posteriors

- `posteriors`:

  list of `K` double 3D arrays with the posterior probability of each
  class (summing to one inside the mask, `0` outside)

- `means`, `sds`, `proportions`:

  the final class intensity means, standard deviations and mixture
  proportions, estimated from the returned posteriors

- `trace`:

  the mean maximum posterior probability over the mask after each
  iteration performed (one entry per iteration, so its length is the
  number of iterations run; see `tolerance`)

The `segmentation` and each posterior carry the `"vox2ras"` attribute
when one is known.

## Details

**Model.** For a voxel \\i\\ inside the mask with intensity \\y_i\\, the
posterior of class \\k\\ follows the default (`Socrates`) posterior
formulation of `Atropos`, \$\$P_i(k) \propto s_k(i)^{w}\\\left(N(y_i
\mid \mu_k, \sigma_k^2)\\ M_k(i)\right)^{1-w},\$\$ where \\w\\ is
`prior_weight`, \$\$s_k(i) = \frac{\pi_k\\ p_k(i)}{\sum_c \pi_c\\
p_c(i)}\$\$ is the spatial prior (the per-voxel normalized prior
\\p_k(i)\\ re-weighted by the mixture proportions \\\pi_k\\) and
\$\$M_k(i) = \frac{\exp\left(\beta \sum_j \omega\_{ij} P_j(k)\right)}{%
\sum_c \exp\left(\beta \sum_j \omega\_{ij} P_j(c)\right)}\$\$ is the
`MRF` term, with \\\beta\\ equal to `mrf_beta`, the sums over the
\\(2r+1)^3 - 1\\ neighbors \\j\\ inside the mask, their posteriors
\\P_j(k)\\ and the inverse physical distance \\\omega\_{ij} = 1 /
d\_{ij}\\ between the two voxels. As in `'ANTs'`, \\s_k\\, \\M_k\\ and
the `Gaussian` density are each floored at `1e-10` before the powers are
taken, so a zero prior is a very strong but not absolute constraint.
Consequently, with `prior_weight = 0` the priors only initialize the
model; with `prior_weight = 1` the image and the `MRF` are ignored and
the segmentation is \\\arg\max_k \pi_k p_k(i)\\; and without `priors`
the spatial prior is the same constant for every class (`Atropos`
applies no prior weight then), so the posterior is proportional to
\\N(y_i \mid \mu_k, \sigma_k^2)\\ M_k(i)\\ and the mixture proportions
do not enter it. Each iteration performs one synchronous E-step (every
voxel reads the posteriors of the previous iteration, initially the
one-hot initial labels) followed by an M-step that re-estimates
\\\mu_k\\ and \\\sigma_k\\ (posterior-weighted mean and unbiased
weighted variance, the `ITK` estimators used by `Atropos`, with a small
variance floor relative to the intensity variance inside the mask) and
\\\pi_k\\ (the mean posterior over the mask).

**Intensity scale.** Because the floor of `1e-10` is applied to the
(unnormalized) `Gaussian` density, whose peak is \\1 /
(\sqrt{2\pi}\\\sigma_k)\\, the result is not invariant to the intensity
scale, exactly as in `Atropos`: once a class standard deviation exceeds
roughly `1e7` the floor starts to flatten its likelihood, and beyond
roughly `1e9` the `EM` inflates the variances until every likelihood is
floored and the segmentation collapses to a single class. Intensities in
such units should be rescaled (or truncated) before segmenting; the
usual `MRI` and `CT` ranges are unaffected.

**Convergence.** After each iteration the mean maximum posterior
probability over the mask (the returned `trace`) is compared with the
previous iteration's value, and the iterations stop as soon as the
signed change is below `tolerance`; the first iteration never stops.
This is the stopping rule of `Atropos`, so with the default
`tolerance = 0` the iterations stop after the first iteration in which
this value decreased, and the posteriors, labels and parameters of that
iteration are returned (there is no rollback to the previous iteration).
The value itself equals the posterior probability that `Atropos` reports
only without `priors` or with `prior_weight = 0`, and with
`mrf_beta = 0`: `Atropos` computes its measure during its
iterated-conditional-modes sweep, with the neighbors' hard labels and
the normalized prior \\p_k(i)\\ itself (not re-weighted by the
proportions) instead of \\s_k(i)\\, so with weighted priors or an active
`MRF` the two measures, and hence the iteration at which the two
implementations stop, can differ.

**Initialization.** Without `priors`, a deterministic `k-means` (Lloyd
iterations seeded at equally spaced intensities between the minimum and
the maximum inside the mask, as in `'ANTs'`) clusters the masked
intensities and the classes are ordered by ascending mean, so with three
classes label `1` is the darkest tissue. A handful of extreme
intensities stretches that range so much that a seed can end up without
any voxel, which is reported as an error naming the intensity range:
truncate or clip the intensities first (the intensity truncation of
`bias_correction_n4` does this), tighten the mask, or reduce
`n_classes`. With `priors`, each voxel starts at the class with the
largest prior (voxels where every prior is zero start unlabeled), the
initial class statistics are weighted by that prior value and the
initial proportions are the label fractions; the class order is the
order of the list.

**Relation to `'ANTs'` `Atropos`.** The posterior formula above, the
`Gaussian` likelihood, the `k-means` initialization, the
inverse-distance neighborhood weights, the probability floor, the
parameter estimators and the stopping rule are those of `Atropos` with
its defaults (the `Socrates` formulation with mixture proportions and no
annealing), so `prior_weight` is interchangeable with its `priorweight`.
The remaining difference is the `MRF` term: `Atropos` plugs the
neighbors' hard labels (updated by an iterated-conditional-modes sweep)
into \\M_k(i)\\, whereas here the neighbors' posteriors are used (a
mean-field update), so the two can differ along tissue boundaries when
`mrf_beta > 0`; with `mrf_beta = 0` the same model is computed, although
the convergence measure (see above) and therefore the number of
iterations can still differ when `prior_weight > 0`. The posteriors
returned here are the ones the labels were taken from, whereas `Atropos`
re-evaluates the probability images it writes out with the parameters
updated in its last iteration. Two safeguards also differ: a class that
receives no voxel at initialization falls back to prior-weighted
statistics instead of a zero likelihood, and priors that are zero
everywhere inside the mask are an error.

**Memory.** All computation is restricted to the mask; the working
memory is about \\2 \times 8 K N\_{mask}\\ bytes for the two posterior
buffers (plus \\8 K N\_{mask}\\ when `prior_weight > 0`) and one integer
per voxel of the full grid.

## References

`Avants` BB, `Tustison` NJ, Wu J, Cook PA, Gee `JC` (2011). An open
source multivariate framework for n-tissue segmentation with evaluation
on public data. *`Neuroinformatics`*, 9(4), 381-400.
[doi:10.1007/s12021-011-9109-y](https://doi.org/10.1007/s12021-011-9109-y)

`Zhang` Y, Brady M, Smith S (2001). Segmentation of brain MR images
through a hidden Markov random field model and the
expectation-maximization algorithm. *IEEE Transactions on Medical
Imaging*, 20(1), 45-57.
[doi:10.1109/42.906424](https://doi.org/10.1109/42.906424)

`Dempster` AP, Laird NM, Rubin DB (1977). Maximum likelihood from
incomplete data via the EM algorithm. *Journal of the Royal Statistical
Society, Series B*, 39(1), 1-38.
[doi:10.1111/j.2517-6161.1977.tb01600.x](https://doi.org/10.1111/j.2517-6161.1977.tb01600.x)

## See also

[`register_volume3d`](https://dipterix.org/ravetools/reference/register_volume3d.md)
to bring template priors into the subject space

## Examples

``` r

# A toy phantom: nested spheres with three intensity classes
nd <- 32
ctr <- (nd - 1) / 2
idx <- arrayInd(seq_len(nd^3), rep(nd, 3)) - 1
r <- sqrt(rowSums((idx - ctr)^2))
truth <- integer(nd^3)
truth[r <= 14] <- 1L
truth[r <= 10] <- 2L
truth[r <= 6] <- 3L
dim(truth) <- rep(nd, 3)
mask <- truth > 0
set.seed(1)
volume <- array(0, dim(truth))
volume[mask] <- c(30, 80, 130)[truth[mask]] + rnorm(sum(mask), sd = 12)

# k-means initialization, classes ordered by intensity
res <- segment_volume_tissue_gmm(volume, mask = mask)
res$means
#> [1]  29.90402  79.58744 129.79314
table(truth = truth[mask], segmented = res$segmentation[mask])
#>      segmented
#> truth    1    2    3
#>     1 7290   22    0
#>     2   21 3286    5
#>     3    0    7  905

# the same with spatial priors (here a softened version of the truth, in a
# shuffled order) that fix the class order
priors <- lapply(c(3, 1, 2), function(k) {
  p <- (truth == k) * 1
  p <- (p + 0.1) / 1.3   # a soft, not quite informative prior
  p
})
res2 <- segment_volume_tissue_gmm(volume, mask = mask, priors = priors,
                                  prior_weight = 0.25)
res2$means       # follows the prior order: bright, dark, medium
#> [1] 129.93406  29.96174  79.68156
```
