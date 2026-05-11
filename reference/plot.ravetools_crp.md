# Plot CRP results

S3 plot method for objects of class `ravetools_crp` returned by
[`crp`](https://dipterix.org/ravetools/reference/crp.md). Produces a
three-panel figure:

1.  Single-trial traces over the full analysis window with the mean and
    the scaled canonical shape \\C(t)\\ overlaid (the shape is drawn
    only up to \\\tau_R\\, so the cut-off is itself informative).

2.  Per-trial \\\alpha'\\ weights sorted in ascending order.

3.  Mean cross-trial projection profile with vertical lines marking
    \\\tau\_{lb}\\, \\\tau_R\\ and \\\tau\_{ub}\\.

## Usage

``` r
# S3 method for class 'ravetools_crp'
plot(x, ...)
```

## Arguments

- x:

  an object of class `ravetools_crp` as returned by
  [`crp`](https://dipterix.org/ravetools/reference/crp.md).

- ...:

  additional graphical parameters passed to
  [`par`](https://rdrr.io/r/graphics/par.html) (e.g. `mar`, `cex.axis`);
  currently unused beyond restoring the previous `par` state on exit.

## Value

Invisibly returns `x`.
