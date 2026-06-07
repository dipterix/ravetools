# Project a surface onto a sphere and relax metric distortion

Projects a closed triangular surface mesh radially onto a sphere and
then iteratively relaxes the metric distortion this introduces
("unfolding"), so that geodesic distances and face orientations stay
close to those of the input surface - a step typically used to prepare
an inflated cortical surface for spherical registration.

## Usage

``` r
mris_sphere(
  mesh,
  target_radius = 100,
  n_averages = 64L,
  niterations = 25L,
  l_dist = 1,
  l_area = 1,
  momentum = 0.9,
  dt = 0.05,
  verbose = FALSE
)
```

## Arguments

- mesh:

  triangular mesh of class `'mesh3d'`, typically an inflated surface
  (see
  [`mris_inflate`](https://dipterix.org/ravetools/reference/mris_inflate.md)).
  Must be closed, manifold, and genus-0 - the same hard precondition as
  [`mris_inflate`](https://dipterix.org/ravetools/reference/mris_inflate.md).
  `mris_sphere` will raise an error if such defects are detected.

- target_radius:

  radius of the target sphere. Default `100`.

- n_averages:

  number of gradient-averaging passes applied to the distance-term
  gradient each iteration (the same neighborhood-averaging mechanism
  [`mris_inflate`](https://dipterix.org/ravetools/reference/mris_inflate.md)
  uses). Default `64`; reduced from the much larger averaging size used
  in production pipelines (tuned for cortical surfaces with hundreds of
  thousands of vertices) to a size more practical for typical `'mesh3d'`
  inputs.

- niterations:

  number of unfolding iterations. Default `25`.

- l_dist:

  distance-preservation coefficient. Default `1.0`.

- l_area:

  folded-face repulsion coefficient. Default `1.0`.

- momentum:

  momentum coefficient. Default `0.9`.

- dt:

  time step. Default `0.05`.

- verbose:

  logical; print per-iteration progress. Default `FALSE`.

## Value

The projected and relaxed surface as a `'mesh3d'` object with `vb`,
`it`, and `normals`.

## Details

**This is a reduced procedure, not a complete reproduction of any
particular reference implementation.** The full unfolding procedure
described in the literature integrates seven weighted energy terms
against a separate reference surface through a multi-resolution,
multi-thousand line optimization pipeline, which is out of scope for
this package. Instead, this implementation keeps the two dominant terms
and integrates them with the same momentum-based machinery
[`mris_inflate`](https://dipterix.org/ravetools/reference/mris_inflate.md)
uses:

- Distance term (`l_dist`):

  Restoring force pulling each vertex back towards the input mesh's own
  original distances to its neighbors.

- Area term (`l_area`):

  Repulsive force that acts only on folded/negative-area faces, pushing
  their vertices apart - this is the actual "unfolding" mechanism.

Both terms are integrated via momentum integration with gradient
averaging and a 1 mm per-step displacement cap, exactly as
[`mris_inflate`](https://dipterix.org/ravetools/reference/mris_inflate.md)
does.

A further simplification: the reference procedure relaxes a freshly
spherical-projected surface against a separate, previously-loaded
white-matter reference surface for the distance term; since this
function takes a single mesh as input, the input mesh's own metric
(captured before projection) is used as the reference instead - the same
convention
[`mris_inflate`](https://dipterix.org/ravetools/reference/mris_inflate.md)
already uses, and the practical analogue (the white-matter surface is
normally what gets inflated to produce a typical input to this kind of
unfolding step).

## Examples

``` r
if (is_not_cran()) {

sphere <- vcg_sphere(sub_division = 3L)

# deform so there is metric distortion to relax
sphere$vb[1, ] <- sphere$vb[1, ] * (1 + 0.2 * rnorm(ncol(sphere$vb)))

# Fix defects
sphere <- vcg_fix_defects(sphere)

result <- mris_sphere(
  sphere,
  n_averages = 8L,
  niterations = 10L,
  target_radius = 1,
  verbose = TRUE
)

plot_mesh_polygon(list(sphere, result),
                  col = list("gray", "red"),
                  alpha = c(0.3, 0.9))



}
#> mris_sphere: building adjacency for 642 vertices, 1280 faces
#> mris_sphere: projecting onto sphere of radius 1.00
#> mris_sphere: starting unfolding, initial neg_area=0.4999 (total_area=12.9821)
#>   iter   5: neg_area=0.53725 (total_area=13.0182)
#>   iter  10: neg_area=0.61820 (total_area=13.0965)
#> mris_sphere: done (10 iterations)

```
