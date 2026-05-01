# Cast rays to intersect with mesh

Cast rays to intersect with mesh

## Usage

``` r
vcg_raycaster(
  x,
  ray_origin,
  ray_direction,
  max_distance = Inf,
  both_sides = FALSE
)
```

## Arguments

- x:

  surface mesh

- ray_origin:

  a matrix with 3 rows or a vector of length 3, the positions of ray
  origin

- ray_direction:

  a matrix with 3 rows or a vector of length 3, the direction of the
  ray, will be normalized to length 1

- max_distance:

  positive maximum distance to cast the normalized ray; default is
  infinity. Any invalid distances (negative, zero, or `NA`) will be
  interpreted as unset.

- both_sides:

  whether to inverse the ray (search both positive and negative ray
  directions); default is false

## Value

A list of ray casting results: whether any intersection is found,
position and face normal of the intersection, distance of the ray, and
the index of the intersecting face (counted from 1)

## Examples

``` r

library(ravetools)
sphere <- vcg_sphere(normals = FALSE)
sphere$vb[1:3, ] <- sphere$vb[1:3, ] + c(10, 10, 10)
vcg_raycaster(
  x = sphere,
  ray_origin = array(c(0, 0, 0, 1, 0, 0), c(3, 2)),
  ray_direction = c(1, 1, 1)
)
#> $has_intersection
#> [1] TRUE TRUE
#> 
#> $intersection
#>          [,1]      [,2]
#> [1,] 9.425265 10.336223
#> [2,] 9.425265  9.336223
#> [3,] 9.425265  9.336223
#> 
#> $normals
#>            [,1]       [,2]
#> [1,] -0.5773502  0.3685302
#> [2,] -0.5773502 -0.6649041
#> [3,] -0.5773503 -0.6496830
#> 
#> $face_index
#> [1] 1217 1179
#> 
#> $distance
#> [1] 16.32504 16.17081
#> 
#> $ray_origin
#>      [,1] [,2]
#> [1,]    0    1
#> [2,]    0    0
#> [3,]    0    0
#> 
#> $ray_direction
#>           [,1]      [,2]
#> [1,] 0.5773503 0.5773503
#> [2,] 0.5773503 0.5773503
#> [3,] 0.5773503 0.5773503
#> 
```
