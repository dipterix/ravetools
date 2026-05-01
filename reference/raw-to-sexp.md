# Convert raw vectors to R vectors

Convert raw vectors to R vectors

## Usage

``` r
raw_to_uint8(x)

raw_to_uint16(x)

raw_to_uint32(x)

raw_to_int8(x)

raw_to_int16(x)

raw_to_int32(x)

raw_to_int64(x)

raw_to_float(x)

raw_to_string(x)
```

## Arguments

- x:

  raw vector of bytes

## Value

Numeric vectors, except for `raw_to_string`, which returns a string.

## Details

For numeric conversions, the function names are straightforward. For
example, `raw_to_uintN` converts raw vectors to unsigned integers, and
`raw_to_intN` converts raw vectors to signed integers. The number `'N'`
stands for the number of bits used to store the integer. For example
`raw_to_uint8` uses 8 bits (1 byte) to store an integer, hence the value
range is `0-255`.

The input data length must be multiple of the element size represented
by the underlying data. For example `uint16` integer uses 16 bites, and
one raw number uses 8 bits, hence two raw vectors can form one unsigned
integer-16. That is, `raw_to_uint16` requires the length of input to be
multiple of two. An easy calculation is: the length of `x` times 8, must
be divided by `'N'` (see last paragraph for definition).

The returned data uses the closest available R native data type that can
fully represent the data. For example, R does not have single `float`
type, hence `raw_to_float` returns `double` type, which can represent
all possible values in `float`. For `raw_to_uint32`, the potential value
range is `0 - (2^32-1)`. This exceeds the limit of R integer type
`(-2^31) - (2^31-1)`. Therefore, the returned values will be real
(double float) data type.

There is no native data type that can store integer-64 data in R,
package `bit64` provides `integer64` type, which will be used by
`raw_to_int64`. Currently there is no solution to convert raw to
unsigned integer-64 type.

`raw_to_string` converts raw to character string. This function respects
`null` character, hence is slightly different than the native
[`rawToChar`](https://rdrr.io/r/base/rawConversion.html), which
translates raw byte-by-byte. If each raw byte represents a valid
character, then the above two functions returns the same result.
However, when the characters represented by raw bytes are invalid,
`raw_to_string` will stop parsing and returns only the valid characters,
while [`rawToChar`](https://rdrr.io/r/base/rawConversion.html) will
still try to parse, and most likely to result in errors. Please see
Examples for comparisons.

## Examples

``` r

# 0x00, 0x7f, 0x80, 0xFF
x <- as.raw(c(0, 127, 128, 255))

raw_to_uint8(x)
#> [1]   0 127 128 255

# The first bit becomes the integer sign
# 128 -> -128, 255 -> -1
raw_to_int8(x)
#> [1]    0  127 -128   -1

## Comments based on little endian system

# 0x7f00 (32512), 0xFF80 (65408 unsigned, or -128 signed)
raw_to_uint16(x)
#> [1] 32512 65408
raw_to_int16(x)
#> [1] 32512  -128

# 0xFF807F00 (4286611200 unsigned, -8356096 signed)
raw_to_uint32(x)
#> [1] 4286611200
raw_to_int32(x)
#> [1] -8356096

# ---------------------------- String ---------------------------

# ASCII case: all valid
x <- charToRaw("This is an ASCII string")

raw_to_string(x)
#> [1] "This is an ASCII string"
rawToChar(x)
#> [1] "This is an ASCII string"

x <- c(charToRaw("This is the end."),
       as.raw(0),
       charToRaw("*** is invalid"))

# rawToChar will raise error
raw_to_string(x)
#> [1] "This is the end."

# ---------------------------- Integer64 ------------------------
# Runs on little endian system
x <- as.raw(c(0x80, 0x00, 0x7f, 0x80, 0xFF, 0x50, 0x7f, 0x00))

# Calculate bitstring, which concaternates the followings
# 10000000 (0x80), 00000000 (0x00), 01111111 (0x7f), 10000000 (0x80),
# 11111111 (0xFF), 01010000 (0x50), 01111111 (0x7f), 00000000 (0x00)

if(.Platform$endian == "little") {
  bitstring <- paste0(
    "00000000011111110101000011111111",
    "10000000011111110000000010000000"
  )
} else {
  bitstring <- paste0(
    "00000001000000001111111000000001",
    "11111111000010101111111000000000"
  )
}

# This is expected value
bit64::as.integer64(structure(
  bitstring,
  class = "bitstring"
))
#> integer64
#> [1] 35836380344942720

# This is actual value
raw_to_int64(x)
#> integer64
#> [1] 35836380344942720

```
