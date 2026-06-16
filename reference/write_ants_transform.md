# Read and write an `'ITK'`/`'ANTs'` `affine` transform

`write_ants_transform` stores a \\4\times 4\\ `RAS`-to-`RAS` `affine`
(such as the `transform` returned by
[`register_volume3d`](https://dipterix.org/ravetools/reference/register_volume3d.md))
as an 'ITK' `.mat` file (a `'MATLAB'` level-5 binary, `LPS` convention)
compatible with `antsApplyTransforms`. `read_ants_transform` reads such
a file back into a \\4\times 4\\ `RAS` matrix. The reader folds a
non-zero `ITK` center of rotation into the translation, so transforms
written by 'ANTs' are read correctly.

## Usage

``` r
write_ants_transform(transform, file)

read_ants_transform(file)
```

## Arguments

- transform:

  a \\4\times 4\\ (or \\3\times 4\\) `RAS` `affine`

- file:

  path to the `.mat` file

## Value

`write_ants_transform` returns `file` invisibly; `read_ants_transform`
returns a \\4\times 4\\ `RAS` matrix.

## See also

[`register_volume3d`](https://dipterix.org/ravetools/reference/register_volume3d.md),
[`save_registration`](https://dipterix.org/ravetools/reference/save_registration.md)

## Examples

``` r
tf <- tempfile(fileext = ".mat")
m <- diag(4); m[1:3, 4] <- c(3, -2, 1)
write_ants_transform(m, tf)
read_ants_transform(tf)
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    0    3
#> [2,]    0    1    0   -2
#> [3,]    0    0    1    1
#> [4,]    0    0    0    1
```
