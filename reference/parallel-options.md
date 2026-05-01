# Set or get thread options

Set or get thread options

## Usage

``` r
detect_threads()

ravetools_threads(n_threads = "auto", stack_size = "auto")
```

## Arguments

- n_threads:

  number of threads to set

- stack_size:

  Stack size (in bytes) to use for worker threads. The default used for
  `"auto"` is 2MB on 32-bit systems and 4MB on 64-bit systems.

## Value

`detect_threads` returns an integer of default threads that is
determined by the number of `CPU` cores; `ravetools_threads` returns
nothing.

## Examples

``` r

detect_threads()
#> [1] 4

ravetools_threads(n_threads = 2)
```
