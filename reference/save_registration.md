# Save or load a registration result in `'ANTs'`-compatible files

`save_registration` writes a
[`register_volume3d`](https://dipterix.org/ravetools/reference/register_volume3d.md)
result to disk as 'ANTs'-style files: an `'ITK'` `affine` `.mat`, the
forward/inverse warp `'NIfTI'`(s) when present, and a small `'DCF'`
manifest (the same key-value format as `'DESCRIPTION'`) recording what
each file is plus the registration parameters. `load_registration` reads
any of those back: a `.mat` returns a \\4\times 4\\ `RAS` matrix, a warp
`.nii`/`.nii.gz` returns the field, and a `.dcf` manifest reassembles
the full object.

## Usage

``` r
save_registration(x, path, ...)

# S3 method for class 'ravetools_register_volume3d'
save_registration(x, path, prefix = "registration", compress = TRUE, ...)

# Default S3 method
save_registration(x, path, prefix = "registration", ...)

load_registration(file, recover_affine_from_header = FALSE)
```

## Arguments

- x:

  a `ravetools_register_volume3d` object (or, for the default method, a
  \\4\times 4\\ matrix)

- path:

  output directory, or the manifest (`.dcf`) path

- ...:

  passed to methods

- prefix:

  file-name prefix; defaults to `"registration"`

- compress:

  logical; write warp fields as `.nii.gz` (default) or `.nii`

- file:

  a `.mat`, `.nii`/`.nii.gz`, or `.dcf` path

- recover_affine_from_header:

  logical; only used as a last resort when loading. If `TRUE` and the
  `affine` `.mat` is missing (or a warp file is loaded directly), the
  `affine` is recovered from the warp header's hidden encoding. Off by
  default because that encoding is reduced-precision

## Value

`save_registration` returns the manifest path invisibly.
`load_registration` returns a \\4\times 4\\ matrix, a field array (with
a `"vox2ras"` attribute), or a `ravetools_register_volume3d` object,
depending on the file type.

## See also

[`register_volume3d`](https://dipterix.org/ravetools/reference/register_volume3d.md),
[`write_ants_transform`](https://dipterix.org/ravetools/reference/write_ants_transform.md),
[`write_ants_warp`](https://dipterix.org/ravetools/reference/write_ants_warp.md)
