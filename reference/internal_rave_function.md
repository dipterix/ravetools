# Get external function from 'RAVE'

Internal function used for examples relative to 'RAVE' project and
should not be used directly.

## Usage

``` r
internal_rave_function(name, pkg, inherit = TRUE, on_missing = NULL)
```

## Arguments

- name:

  function or variable name

- pkg:

  'RAVE' package name

- inherit:

  passed to [`get0`](https://rdrr.io/r/base/exists.html)

- on_missing:

  default value to return of no function is found

## Value

Function object if found, otherwise `on_missing`.
