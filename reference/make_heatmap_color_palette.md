# Create Heatmap Color Palette

This function generates a color palette suitable for heatmaps based on
the values in a matrix. It uses the \`createLinearColors\` function to
generate separate color gradients for positive and negative values.

## Usage

``` r
make_heatmap_color_palette(my_matrix)
```

## Arguments

- my_matrix:

  A numeric matrix for which the heatmap color palette is to be
  generated.

## Value

A character vector of colors representing the heatmap color palette
based on the input matrix values.

## Examples

``` r
# Create a sample matrix
my_matrix <- matrix(c(-3, -1, 0, 1, 3), nrow = 1)

# Generate heatmap color palette
heatmap_palette <- make_heatmap_color_palette(my_matrix)
```
