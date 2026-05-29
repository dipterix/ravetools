# Render one or more meshes as flat-shaded triangles in base R

Projects each triangular face onto a 2D plane using an orthographic
camera, shades it with a single color proportional to how directly it
faces the camera (Lambert), depth-sorts all faces across all meshes, and
draws them in a single
[`polygon`](https://rdrr.io/r/graphics/polygon.html) call.

Meshes without faces (point clouds) are substituted by a small sphere
([`vcg_sphere`](https://dipterix.org/ravetools/reference/vcg_sphere.md))
centered at each point and scaled by `cex`; they then participate in the
same rasterizer as ordinary faced meshes.

A camera-facing clipping pass discards triangles whose outward normal
points along the camera ray (signed \\n \cdot z\_{cam}\\ \\\> 1 -
\mathrm{mesh\\clipping}\\), peeling the front cap off the surface so the
back wall (and any interior meshes) become visible. Set
`mesh_clipping = 0` to disable clipping. Point-cloud meshes (those
rendered as substitute
[`vcg_sphere`](https://dipterix.org/ravetools/reference/vcg_sphere.md)
instances) are exempt from this clip so they remain solid even when the
enclosing surface is peeled.

Multiple meshes share a single depth space: all faces are projected,
sorted, and drawn together so the painter's algorithm works correctly
across meshes.

## Usage

``` r
plot_mesh_polygon(
  mesh,
  eye = c(0, 0, 1000),
  lookat = c(0, 0, 0),
  up = c(0, 1, 0),
  col = c("white", "gray30"),
  cex = 1,
  add = FALSE,
  axes = FALSE,
  asp = 1,
  xlim = NULL,
  ylim = NULL,
  xlab = "",
  ylab = "",
  side = c("both", "front", "back"),
  mesh_clipping = 0,
  sphere_subdivision = 1L,
  alpha = 1,
  shadow_color = NULL,
  light_intensity = 1,
  ambient_intensity = 0.2,
  ...
)
```

## Arguments

- mesh:

  a `'mesh3d'` object, or a list of `'mesh3d'` objects. Meshes without a
  face matrix (`$it`) are rendered as sphere instances (one
  [`vcg_sphere`](https://dipterix.org/ravetools/reference/vcg_sphere.md)
  per vertex, radius `cex`).

- eye:

  numeric vector of length 3 - camera position in world space.

- lookat:

  numeric vector of length 3 - the world-space point the camera is
  looking at.

- up:

  numeric vector of length 3 - world-space "up" direction; defaults to
  `c(0, 1, 0)`.

- col:

  base color(s) per mesh. Same forms as
  [`plot_mesh_dotcloud`](https://dipterix.org/ravetools/reference/plot_mesh_dotcloud.md):
  a single color, a depth-gradient vector, a per-vertex character
  vector, or a list of any of these (one element per mesh). Default
  `c("white", "gray30")`.

- cex:

  radius of the substitute sphere used for point-cloud meshes (world
  units). Has no effect on meshes that already have faces. Default `1`.

- add:

  logical; if `TRUE` the faces are added to an existing plot instead of
  opening a new one. Default `FALSE`.

- axes, asp, xlim, ylim, xlab, ylab:

  passed to
  [`plot.default`](https://rdrr.io/r/graphics/plot.default.html) when
  `add = FALSE`.

- side:

  which side of each triangle to render. One of `"both"` (default, shows
  all triangles), `"front"`, or `"back"`.

- mesh_clipping:

  numeric in `[0, 1]` controlling camera-facing clipping: any triangle
  whose signed Lambert dot product (\\n \cdot z\_{cam}\\, i.e. positive
  for front-facing) exceeds `1 - mesh_clipping` is discarded. A value of
  `0` keeps the full surface; `0.5` peels a 60-degree cap off the front
  and reveals the back wall (whose absolute Lambert shade is still
  large, so it renders brightly); `1` clips every front-facing triangle
  and shows the interior only. Back-facing triangles are never clipped
  by this rule. Default `0`.

- sphere_subdivision:

  integer subdivision level forwarded to
  [`vcg_sphere`](https://dipterix.org/ravetools/reference/vcg_sphere.md)
  when substituting point clouds. Higher values give smoother spheres at
  the cost of more triangles. Default `1` (80 triangles per sphere).

- alpha:

  numeric in `[0, 1]`; one value per mesh (recycled), sets the
  transparency of each mesh (`1` = fully opaque, `0` = fully
  transparent). Faces belonging to a mesh with `alpha = 0` are dropped
  entirely. Because faces are drawn in back-to-front (painter's) order,
  alpha blending across multiple meshes follows that order. Default `1`.

- shadow_color:

  color used for fully unlit (grazing/back) faces. The Lambert shade
  linearly interpolates from `shadow_color` (at `shade = 0`) toward the
  face color lit by a white light (at `shade = 1`), so unlit faces fade
  to this color instead of an implicit black background. Default
  `par("fg")`, so the shadow contrasts with the current device's canvas.

- light_intensity:

  non-negative scalar controlling the brightness of the (white) light
  source: at `shade = 1` the face color is multiplied by
  `light_intensity` (clamped to `[0, 1]`). `1` (default) reproduces the
  face color exactly when fully lit; smaller values darken the whole
  mesh, larger values saturate.

- ambient_intensity:

  scalar in `[0, 1]` acting as a lower bound on the Lambert shade so
  grazing/back faces still receive some light and never collapse fully
  to `shadow_color`. An effective shade of
  `max(shade, ambient_intensity)` is used in the bg/light lerp. Default
  `0.2`.

- ...:

  additional graphical parameters forwarded to
  [`plot.default`](https://rdrr.io/r/graphics/plot.default.html) (new
  plot only).

## Value

Invisibly returns a list with components `xlim` and `ylim` (the plot
limits used).

## Details

Limitations of the base-R polygon path (no `rgl`):

- Flat shading only - one color per triangle. No per-vertex color
  interpolation.

- No depth buffer - faces are depth-sorted by centroid (painter's
  algorithm). Interpenetrating triangles can render in the wrong order.

- Anti-aliasing seams can appear between adjacent triangles on raster
  devices; Cairo-based devices (`png(type = "cairo")`,
  [`svg()`](https://rdrr.io/r/grDevices/cairo.html),
  [`pdf()`](https://rdrr.io/r/grDevices/pdf.html)) produce cleaner
  output than the default quartz/X11 path.

## See also

[`plot_mesh_dotcloud`](https://dipterix.org/ravetools/reference/plot_mesh_dotcloud.md),
[`vcg_isosurface`](https://dipterix.org/ravetools/reference/vcg_isosurface.md)

## Examples

``` r

mesh <- vcg_isosurface(left_hippocampus_mask)

# Surface alone
plot_mesh_polygon(
  mesh,
  eye    = c(150, 0, 0),
  lookat = c(0, 0, 0),
  up     = c(0, 0, 1),
  col    = "steelblue"
)


# Surface + electrode point cloud (rendered as small icospheres)
n_elec <- 20
electrodes <- structure(
  list(vb = matrix(rnorm(3 * n_elec, sd = 5), 3, n_elec)),
  class = "mesh3d"
)
plot_mesh_polygon(
  mesh = list(mesh, electrodes),
  eye    = c(150, 0, 0),
  lookat = c(0, 0, 0),
  up     = c(0, 0, 1),
  col    = list("steelblue", "red"),
  cex    = 1.5
)

```
