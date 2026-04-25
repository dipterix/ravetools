# Safe ways to call package `'rgl'` without requiring `'x11'`

Internally used for example show-cases. Please install package `'rgl'`
manually to use these functions.

## Usage

``` r
rgl_call(FUN, ...)

rgl_view(expr, quoted = FALSE, env = parent.frame())

rgl_plot_normals(x, length = 1, lwd = 1, col = 1, ...)
```

## Arguments

- FUN:

  `'rgl'` function name

- ...:

  passed to `'rgl'` function

- expr:

  expression within which `'rgl'` functions are called

- quoted:

  whether `expr` is quoted

- env:

  environment in which `expr` is evaluated

- x:

  triangular `'mesh3d'` object

- length, lwd, col:

  normal vector length, size, and color

## Examples

``` r

# Make sure the example does not run when compiling
# or check the package
if(FALSE) {

  volume <- array(0, dim = c(8,8,8))
  volume[4:5, 4:5, 4:5] <- 1
  mesh <- mesh_from_volume(volume, verbose = FALSE)

  rgl_view({

    rgl_call("shade3d", mesh, col = 3)
    rgl_plot_normals(mesh)

  })

}

```
