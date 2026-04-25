# Create a `Matrix4` instance for `'Affine'` transform

Create a `Matrix4` instance for `'Affine'` transform

## Usage

``` r
new_matrix4()

as_matrix4(m)
```

## Arguments

- m:

  a matrix or a vector to be converted to the `Matrix4` instance; `m`
  must be one of the followings: for matrices, the dimension must be
  `4x4`, `3x4` (the last row will be `0 0 0 1`), or `3x3` (linear
  transform); for vectors, the length must be `16`, `12` (will append
  `0 0 0 1` internally), `3` (translation), or `1` (scale).

## Value

A `Matrix4` instance

## See also

[`new_vector3`](https://dipterix.org/ravetools/reference/new_vector3.md),
[`new_quaternion`](https://dipterix.org/ravetools/reference/new_quaternion.md)
