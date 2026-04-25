# Create a `Vector3` instance to store '3D' points

Create instances that mimic the `'three.js'` syntax.

## Usage

``` r
new_vector3(x = 0, y = 0, z = 0)

as_vector3(v)
```

## Arguments

- x, y, z:

  numeric, must have the same length, `'xyz'` positions

- v:

  R object to be converted to `Vector3` instance

## Value

A `Vector3` instance

## See also

[`new_matrix4`](https://dipterix.org/ravetools/reference/new_matrix4.md),
[`new_quaternion`](https://dipterix.org/ravetools/reference/new_quaternion.md)

## Examples

``` r
vec3 <- new_vector3(
  x = 1:9,
  y = 9:1,
  z = rep(c(1,2,3), 3)
)

vec3[]
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#> [1,]    1    2    3    4    5    6    7    8    9
#> [2,]    9    8    7    6    5    4    3    2    1
#> [3,]    1    2    3    1    2    3    1    2    3

# transform
m <- new_matrix4()

# rotation xy plane by 30 degrees
m$make_rotation_z(pi / 6)
#> <Matrix4>
#>           [,1]       [,2] [,3] [,4]
#> [1,] 0.8660254 -0.5000000    0    0
#> [2,] 0.5000000  0.8660254    0    0
#> [3,] 0.0000000  0.0000000    1    0
#> [4,] 0.0000000  0.0000000    0    1

vec3$apply_matrix4(m)
#> <Vector3: len=9>
#>            x        y z
#> 1 -3.6339746 8.294229 1
#> 2 -2.2679492 7.928203 2
#> 3 -0.9019238 7.562178 3
#> 4  0.4641016 7.196152 1
#> 5  1.8301270 6.830127 2
#> ...

vec3[]
#>           [,1]      [,2]       [,3]      [,4]     [,5]     [,6]     [,7]
#> [1,] -3.633975 -2.267949 -0.9019238 0.4641016 1.830127 3.196152 4.562178
#> [2,]  8.294229  7.928203  7.5621778 7.1961524 6.830127 6.464102 6.098076
#> [3,]  1.000000  2.000000  3.0000000 1.0000000 2.000000 3.000000 1.000000
#>          [,8]     [,9]
#> [1,] 5.928203 7.294229
#> [2,] 5.732051 5.366025
#> [3,] 2.000000 3.000000

as_vector3(c(1,2,3))
#> <Vector3: len=1>
#>   x y z
#> 1 1 2 3
```
