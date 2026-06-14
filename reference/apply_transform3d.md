# Apply a linear `RAS` transform to resample a 3D volume

Warps a moving volume onto a reference grid given a 4x4 `RAS`-to-`RAS`
transform (such as the `transform` returned by
[`register_volume3d`](https://dipterix.org/ravetools/reference/register_volume3d.md)).
The transform maps reference (fixed) `RAS` coordinates to moving `RAS`
coordinates.

## Usage

``` r
apply_transform3d(
  volume,
  vox2ras,
  transform,
  reference_dim = dim(volume),
  reference_vox2ras = vox2ras,
  interpolation = c("trilinear", "nearest", "bspline"),
  na_fill = 0
)
```

## Arguments

- volume:

  moving 3D array to resample

- vox2ras:

  the moving volume's voxel-to-`RAS` 4x4 transform

- transform:

  4x4 `RAS`-to-`RAS` transform (fixed to moving)

- reference_dim:

  output dimension (the fixed grid); defaults to `dim(volume)`

- reference_vox2ras:

  the fixed grid's voxel-to-`RAS` transform; defaults to `vox2ras`

- interpolation:

  `'trilinear'` (default), `'nearest'`, or `'bspline'` (cubic
  `Catmull-Rom`)

- na_fill:

  value for out-of-bounds voxels; default `0`

## Value

The resampled volume on the reference grid, with a `'vox2ras'` attribute
equal to `reference_vox2ras`.

## See also

[`register_volume3d`](https://dipterix.org/ravetools/reference/register_volume3d.md)
