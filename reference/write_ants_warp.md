# Read and write an `'ANTs'` deformation (warp) field

`write_ants_warp` stores a dense displacement field (such as the
`forward_field` / `inverse_field` from
[`register_volume3d`](https://dipterix.org/ravetools/reference/register_volume3d.md))
as an 'ANTs' 5-D warp `'NIfTI'` (`intent_code` `1007`, `LPS` vectors)
readable by `antsApplyTransforms`; `read_ants_warp` reads it back.
Requires the suggested freesurferformats package.

The field is stored on disk with displacement vectors in `LPS` (the x
and y components are negated relative to the `RAS` field) and an `sform`
equal to `vox2ras`, so the file describes the fixed grid's physical
space exactly as 'ANTs' expects. Optionally the `RAS` `affine` can be
hidden in the header text fields (`descrip`/`aux_file`); 'ANTs' ignores
those, and recovery on read is opt-in.

## Usage

``` r
write_ants_warp(
  field,
  file,
  vox2ras = attr(field, "vox2ras"),
  affine = NULL,
  direction = c("forward", "inverse")
)

read_ants_warp(file, recover_affine = FALSE)
```

## Arguments

- field:

  a `(nx, ny, nz, 3)` array of `RAS` displacements; the `vox2ras`
  attribute is used when `vox2ras` is not supplied

- file:

  path to the warp `.nii` / `.nii.gz` file

- vox2ras:

  the fixed-grid \\4\times 4\\ voxel-to-`RAS` matrix

- affine:

  optional \\4\times 4\\ `RAS` `affine` to embed (with reduced
  precision) in the header; `NULL` (default) embeds nothing

- direction:

  `"forward"` or `"inverse"`; recorded as the header direction flag when
  `affine` is embedded

- recover_affine:

  logical; if `TRUE`, attempt to decode a hidden `affine` from the
  header and attach it as the `"transform"` attribute (with a
  `"direction"` attribute). Default `FALSE`

## Value

`write_ants_warp` returns `file` invisibly; `read_ants_warp` returns the
`(nx, ny, nz, 3)` `RAS` field with a `"vox2ras"` attribute (and, if
recovered, `"transform"` / `"direction"` attributes).

## See also

[`register_volume3d`](https://dipterix.org/ravetools/reference/register_volume3d.md),
[`save_registration`](https://dipterix.org/ravetools/reference/save_registration.md)
