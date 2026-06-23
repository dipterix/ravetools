# Plot electrode clustering results

S3 plot method for `ravetools_crp_cluster` objects returned by
[`crp_cluster`](https://dipterix.org/ravetools/reference/crp_cluster.md).
Draws a two-panel figure: the per-cluster basis profile curves, and the
electrode-by-electrode similarity matrix ordered by cluster.

## Usage

``` r
# S3 method for class 'ravetools_crp_cluster'
plot(x, ...)
```

## Arguments

- x:

  a `ravetools_crp_cluster` object.

- ...:

  ignored.

## Value

Invisibly returns `x`.
