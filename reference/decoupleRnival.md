# DecoupleRnival

Iteratively propagate downstream input activity through a signed
directed network using the weighted mean enrichment score from decoupleR
package

## Usage

``` r
decoupleRnival(
  upstream_input = NULL,
  downstream_input,
  meta_network,
  n_layers,
  n_perm = 1000,
  downstream_cutoff = 0,
  statistic = "norm_wmean"
)
```

## Arguments

- upstream_input:

  A named vector with up_stream nodes and their corresponding activity.

- downstream_input:

  A named vector with down_stream nodes and their corresponding
  activity.

- meta_network:

  A network data frame containing signed directed prior knowledge of
  molecular interactions.

- n_layers:

  The number of layers that will be propagated upstream.

- n_perm:

  The number of permutations to use in decoupleR's algorithm.

- downstream_cutoff:

  If downstream measurments should be included above a given threshold

- statistic:

  the decoupleR stat to consider: "wmean", "norm_wmean", or "ulm"

## Value

A data frame containing the score of the nodes upstream of the
downstream input based on the iterative propagation

## Examples

``` r
# Example input data
upstream_input <- c("A" = 1, "B" = -1, "C" = 0.5)
downstream_input <- c("D" = 2, "E" = -1.5)
meta_network <- data.frame(
  source = c("A", "A", "B", "C", "C", "D", "E"),
  target = c("B", "C", "D", "E", "D", "B", "A"),
  sign = c(1, -1, -1, 1, -1, -1, 1)
)

# Run the function with the example input data
result <- decoupleRnival(upstream_input, downstream_input, meta_network, n_layers = 2, n_perm = 100)
#> [1] "Warning, this function is deprecated and will no longer receive futur support. Please use the 'moon' function instead"
#> [1] 2

# View the results
print(result)
#>   source    score
#> 1      B -1.05659
#> D      D  2.00000
#> E      E -1.50000
```
