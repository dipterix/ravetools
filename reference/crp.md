# Canonical Response Parameterization (`CRP`)

Parameterizes single-trial evoked responses (e.g. cortico-cortical
evoked potentials, `CCEPs`) using the Canonical Response
Parameterization method (see 'Citation'). The function estimates the
response duration \\\tau_R\\, the time after stimulus at which the
evoked response has its most consistent, shared structure across trials.
The estimator is obtained from the time course of cross-trial projection
magnitudes, extracts the canonical response shape \\C(t)\\ via a linear
kernel-trick PCA on the trial matrix truncated at \\\tau_R\\, and
reports per-trial weights, residuals, signal-to-noise, explained
variance and extraction-significance statistics.

This is an R translation of `CRP_method.m` (and the surrounding
artifact-rejection / duration-uncertainty logic in `CRP_illustration.m`)
from the upstream MATLAB reference implementation; see ‘References’.

## Usage

``` r
crp(
  x,
  time,
  t_start = 0.015,
  t_end = 1,
  remove_artifacts = TRUE,
  artifact_interval = c("full", "tR"),
  artifact_p_threshold = 1e-05,
  threshold_quantile = 0.98,
  time_step = 5L
)
```

## Arguments

- x:

  numeric matrix of single-trial evoked voltages with shape
  `time x trials` (i.e. rows are timepoints, columns are trials); the
  matrix orientation matches the variable `V` / `data` in the MATLAB
  reference. At least two trials are required.

- time:

  numeric vector of length `nrow(x)` giving the stimulus-aligned time
  (in seconds) of each row of `x`; must be monotonically increasing and
  span `[t_start, t_end]`.

- t_start, t_end:

  numeric scalars, post-stimulation start and end times (in seconds)
  defining the analysis window. Defaults match the MATLAB illustration
  (`0.015 s` to `1 s`).

- remove_artifacts:

  logical; if `TRUE` (the default), an initial `CRP` pass is run to
  identify outlier/artifact trials, which are then dropped before the
  final pass. See ‘Details’.

- artifact_interval:

  character, one of `"full"` (the default, matching the active option in
  the MATLAB illustration) or `"tR"`; selects whether per-trial outlier
  statistics are computed on the projection magnitudes for the full
  window or only at the response duration \\\tau_R\\.

- artifact_p_threshold:

  numeric, p-value threshold below which a trial is flagged as artifact
  (provided its mean projection is also below the cohort mean); defaults
  to `1e-5`.

- threshold_quantile:

  numeric in `(0, 1)`; the fraction of the peak mean projection
  magnitude used to derive the duration-uncertainty bounds `tau_R_lower`
  and `tau_R_upper`. Defaults to `0.98` as in the manuscript.

- time_step:

  integer, sampling step (in samples) used when sweeping candidate
  response duration; defaults to `5L`, matching the MATLAB `t_step`.
  Larger values are faster but smooth the projection profile.

## Value

A named list with the following elements:

- `parameters`:

  a list of single-trial parameterizations (`crp_params` in MATLAB),
  with fields: `V_tR` (voltage matrix truncated to \\\tau_R\\), `al`
  (alpha coefficient weights of \\C\\ into each trial), `C` (canonical
  shape, the first kernel-PCA component), `ep` (residual epsilon after
  removing \\\alpha C\\ from each trial), `tR` (response duration in
  seconds), `params_times` (times for parameterized data),
  `avg_trace_tR` (simple average trace, truncated to \\\tau_R\\), `al_p`
  (alpha-prime: \\\alpha\\ scaled by \\\sqrt{\textrm{length}(C)}\\),
  `epep_root` (per-trial residual norm), `Vsnr` (per-trial
  signal-to-noise), `expl_var` (per-trial explained variance).

- `projections`:

  a list of projection-stage outputs (`crp_projs` in MATLAB), with
  fields: `proj_tpts` (projection times in seconds), `S_all` (matrix of
  projection magnitudes; rows are non-redundant trial pairs, columns are
  durations), `mean_proj_profile`, `var_proj_profile`, `tR_index`
  (column of `S_all` corresponding to \\\tau_R\\), `avg_trace_input`
  (mean trace over the full window), `stat_indices` (row indices of
  `S_all` used for the significance t-tests; non-overlapping comparison
  pairs), `t_value_tR`, `p_value_tR`, `t_value_full`, `p_value_full`.

- `bad_trials`:

  integer vector of trial indices (into the original `x`) flagged and
  removed as artifacts; `integer(0)` when none are removed or when
  `remove_artifacts = FALSE`.

- `tau_R_lower`, `tau_R_upper`:

  numeric scalars, lower and upper threshold-crossing times (in seconds)
  bracketing \\\tau_R\\ at the `threshold_quantile` fraction of the peak
  mean projection magnitude.

- `t_start`, `t_end`:

  the analysis window used.

- `sample_rate`:

  numeric, sampling rate inferred from `time`.

## Details

Briefly, the algorithm proceeds in three stages:

1.  For a sweep of candidate durations \\k\\, compute pairwise
    L2-normalized cross-projection magnitudes between trials truncated
    to \\\[0, k\]\\. The duration that maximizes the mean projection
    magnitude is taken as the response duration \\\tau_R\\.

2.  Apply linear kernel-trick PCA to the trial matrix truncated to
    \\\tau_R\\; the first principal component is the canonical response
    shape \\C(t)\\.

3.  Project \\C(t)\\ into each trial to obtain per-trial weights
    \\\alpha_k\\; the residual \\\epsilon_k = V_k - \alpha_k C\\
    summarizes trial-by-trial deviation from the canonical shape.

Significance is assessed by a one-sided t-test on the off-diagonal
projection magnitudes against zero, restricted to a non-overlapping
subset of comparison pairs to avoid double-counting.

When `remove_artifacts = TRUE`, the function performs an initial `CRP`
pass and runs an unpaired t-test for each trial comparing the
projections it participates in against all other off-diagonal
projections. Trials with `p < artifact_p_threshold` *and* mean
projection below the cohort mean are dropped, and `CRP` is re-run.

## References

The `CRP` algorithm is described in
[doi:10.1371/journal.pcbi.1011105](https://doi.org/10.1371/journal.pcbi.1011105)
, with a reference MATLAB implementation at
<https://github.com/kaijmiller/crp_scripts>. See `citation("ravetools")`
for the full bibliographic entry.

## Examples

``` r


set.seed(42)

# Synthetic CCEP-like data: shared canonical shape with per-trial scaling
n_time <- 500L
n_trials <- 20L
tt <- seq(-0.05, 1, length.out = n_time)
canonical <- exp(-((tt - 0.10) / 0.05)^2) -
             0.5 * exp(-((tt - 0.30) / 0.10)^2)
V <- outer(canonical, runif(n_trials, 0.5, 1.5)) +
     matrix(rnorm(n_time * n_trials, sd = 0.3), n_time, n_trials)
V <- V * 100

res <- crp(V, tt)
res$parameters$tR
#> [1] 0.3918838
res$projections$p_value_tR
#> [1] 4.183469e-90


matplot(tt, V, type = 'l', lty = 1, col = "#80808080", las = 1,
        xlab = "Time (s)", ylab = bquote(mu * V),
        main = bquote(tau[R] ~ ": estimated ERP duration"))

# Mean response excluding bad trials
V_mean <- rowMeans(V[, -res$bad_trials])
lines(tt, V_mean, col = "black", lwd = 2)

# Plot CRP response stopped at tR
sel_lower <- tt <= res$tau_R_lower & tt > 0
sel_tR <- tt <= res$parameters$tR & tt > 0
sel_upper <- tt <= res$tau_R_upper & tt > 0

lines(tt[sel_upper], V_mean[sel_upper], col = "red", lwd = 2)
lines(tt[sel_tR], V_mean[sel_tR], col = "orange2", lwd = 2)
lines(tt[sel_lower], V_mean[sel_lower], col = "cyan", lwd = 2)

idx_lower <- which.min(abs(tt - res$tau_R_lower))
idx_tR <- which.min(abs(tt - res$parameters$tR))
idx_upper <- which.min(abs(tt - res$tau_R_upper))
idx <- c(idx_lower, idx_tR, idx_upper)
points(tt[idx], V_mean[idx], col = c("cyan", "orange2", "red"))
text(tt[idx], V_mean[idx] + 20,
     labels = expression(tau[lb], tau[R], tau[up]),
     adj = c(0.5, 0))

# The underlying duration
underlying_duration <- max(tt[which(abs(canonical) > 0.1)])
abline(v = underlying_duration, lty = 3)
text(x = underlying_duration, y = 100,
     labels = " Underlying duration", adj = 0, cex = 0.6)


```
