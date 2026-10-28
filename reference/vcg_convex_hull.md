# Convex hull of a point cloud or surface

Computes the smallest convex closed surface that contains a set of
points (for a surface, its vertices), with the `Quickhull`
implementation of the bundled `VCG` library or with a port of the
`quickhull3d` library to 'C++'.

## Usage

``` r
vcg_convex_hull(x, method = c("vcg", "quickhull3d"), tolerance = NA)
```

## Arguments

- x:

  a matrix of point coordinates with 3 columns and one row per point, or
  a surface (see 'Coercing Surface Inputs' below), in which case its
  vertices are used and its faces ignored

- method:

  algorithm: `"vcg"` (default) for the `Quickhull` of `VCG` (`Barber`
  and colleagues, 1996), or `"quickhull3d"` for a port of the
  `quickhull3d` library, which follows `Lloyd`'s `QuickHull3D` and
  merges neighboring faces that lie in nearly the same plane; see
  'Details'

- tolerance:

  non-negative distance, in the units of `x`: a point closer than this
  to the plane of a hull face counts as lying on the face rather than
  outside it. The default `NA` uses the method's own value: `1e-10`
  times the diagonal of the bounding box for `"vcg"`, and for
  `"quickhull3d"` the rounding error of the coordinates, three times the
  machine precision times the sum, over the three axes, of the largest
  absolute coordinate. Values below that rounding error are raised to
  it, because below it the visibility tests contradict each other

## Value

A closed triangular mesh of class `'mesh3d'` with:

- `vb`:

  coordinates of the hull vertices, each one exactly an input point, in
  the input order

- `it`:

  triangular faces, ordered counter-clockwise when seen from outside
  (face normals point outward)

- `normals`:

  area-weighted vertex normal vectors

and attribute `"vertex_index"`: for each hull vertex, the row of `x` (or
the vertex index of the surface) it comes from, increasing;
`attr(hull, "vertex_index")[hull$it]` indexes the faces by input point.

## Details

Both methods start from a tetrahedron of extreme points and repeatedly
add the point farthest outside a face, replacing the faces that point
can see. Points within `tolerance` of a face plane count as on the face.

The `VCG` implementation keeps every face a triangle, so a flat region
of the hull (common for points on a voxel grid) is split into arbitrary
triangles, and it relies on the tolerance alone to stay consistent under
rounding. The `quickhull3d` method also merges neighboring faces that
lie in one plane or are not convex within the tolerance into polygons,
which it splits into triangles fanning out from one corner; such a fan
can contain triangles of zero area when a polygon has three corners on a
line.

Points that are identical, or all lie on one line or in one plane,
enclose no volume and raise an error. So does a result that fails the
check for a closed, consistently oriented surface, a numerical failure
that a larger `tolerance` or the other method can avoid.

## Coercing Surface Inputs

The surface objects are converted to `'mesh3d'` object before applying
further calculations.

When `surface` is a surface ieegio object, the returned `mesh3d$vb`
contains vertices that have been left-multiplied by
`surface$geometry$transforms[[1]]` (the first transform stored in the
geometry, typically the `ScannerAnat` or voxel-to-world transform).

**Breaking change:** Earlier versions (before 0.2.6) of ravetools
returned the raw `surface$geometry$vertices` without applying any
transform, so downstream code often multiplied by
`surface$geometry$transforms[[1]]` (or an equivalent) manually before
working in world space. Such code will now *double* apply the transform
and produce incorrect coordinates. If you previously applied a transform
from `surface$geometry$transforms` by hand after calling a ravetools
mesh function on an `'ieegio_surface'`, remove that manual step.

Surfaces with an empty or missing `geometry$transforms` list (for
example, surfaces produced by ieegio's `volume_to_surface`, which stores
an identity transform) are unaffected.

If `geometry$transforms` contains multiple transforms targeting
different coordinate spaces, only the first one is used. Callers that
need a specific target space should select and apply that transform
themselves before calling ravetools mesh functions.

## References

`Barber` CB, `Dobkin` DP, `Huhdanpaa` H (1996). The `quickhull`
algorithm for convex hulls. *`ACM` Transactions on Mathematical
Software*, 22(4), 469-483.

`Lloyd JE`. `QuickHull3D`: a robust 3D convex hull algorithm in Java.
<https://www.cs.ubc.ca/~lloyd/java/quickhull3d.html>

`Poppe` M. `quickhull3d`: a robust `quickhull` implementation to find
the convex hull of a set of 3D points.
<https://github.com/mauriciopoppe/quickhull3d>

## Examples

``` r

# the 8 corners of a unit cube and 100 points inside it
set.seed(1)
x <- rbind(
  as.matrix(expand.grid(0:1, 0:1, 0:1)),
  matrix(runif(300), ncol = 3)
)
hull <- vcg_convex_hull(x)

# the rows of `x` on the hull are the corners
attr(hull, "vertex_index")
#> [1] 1 2 3 4 5 6 7 8

vcg_mesh_volume(hull)
#> [1] 1

if(is_not_cran()) {

data("left_hippocampus_mask")
mesh <- vcg_isosurface(left_hippocampus_mask)
hull <- vcg_convex_hull(mesh, method = "quickhull3d")

rgl_view({
  rgl_call("shade3d", mesh, col = 2)
  rgl_call("wire3d", hull, col = 3)
})

}
#> Package `rgl` is not installed. Please install `rgl` to use this function.
#> Error in loadNamespace(name): there is no package called ‘rgl’
```
