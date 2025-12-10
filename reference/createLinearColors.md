# Create Linear Colors Based on Numeric Input

This function generates a gradient of colors based on the provided
numeric values. The colors can be adjusted to include zero and are
configurable with a specified maximum and custom color palette.

## Usage

``` r
createLinearColors(
  numbers,
  withZero = T,
  maximum = 100,
  my_colors = c("royalblue3", "white", "red")
)
```

## Arguments

- numbers:

  A numeric vector for which the color gradient is to be generated.

- withZero:

  A logical value indicating whether zero should be included in the
  color gradient. Default is TRUE.

- maximum:

  An integer specifying the maximum number of colors to be generated in
  the gradient. Default is 100.

- my_colors:

  A character vector of length three specifying the colors to be used in
  the gradient. Default is c("royalblue3", "white", "red").

## Value

A character vector of colors representing the gradient based on the
input numeric values.

## Examples

``` r
# Generate colors for a set of numbers including zero
numbers <- c(-50, -20, 0, 20, 50)
colors <- createLinearColors(numbers, withZero = TRUE, maximum = 100)
```
