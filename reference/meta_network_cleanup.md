# meta_network_cleanup

This function cleans up a meta network data frame by removing
self-interactions, calculating the mean interaction values for
duplicated source-target pairs, and keeping only interactions with
values of 1 or -1.

## Usage

``` r
meta_network_cleanup(meta_network)
```

## Arguments

- meta_network:

  A data frame with columns 'source', 'interaction', and 'target'.

## Value

A cleaned up meta network data frame.

## Examples

``` r
# Create a meta network data frame
example_meta_network <- data.frame(
source = c("A", "B", "C", "D", "A", "B", "C", "D", "A"),
interaction = c(1, 1, 1, -1, 1, -1, 1, -1, 1),
target = c("B", "C", "D", "A", "C", "D", "A", "B", "B")
)

# Clean up the example meta network
cleaned_meta_network <- meta_network_cleanup(example_meta_network)
print(cleaned_meta_network)
#>   source target interaction
#> 1      A      B           1
#> 2      A      C           1
#> 3      B      C           1
#> 4      B      D          -1
#> 5      C      A           1
#> 6      C      D           1
#> 7      D      A          -1
#> 8      D      B          -1
```
