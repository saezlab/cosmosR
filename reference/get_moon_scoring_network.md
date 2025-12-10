# Get Moon Scoring Network

This function analyzes a given meta network based on moon scores and an
upstream node. It filters and processes the network by controlling and
observing neighbours according to specified parameters. The function
returns a list containing a filtered network and updated moon scores.

## Usage

``` r
get_moon_scoring_network(
  upstream_node,
  meta_network,
  moon_scores,
  keep_upstream_node_peers = F
)
```

## Arguments

- upstream_node:

  The node from which the network analysis starts.

- meta_network:

  The complete network data.

- moon_scores:

  Scores associated with each node in the network.

- keep_upstream_node_peers:

  Logical; whether to keep peers of the upstream node. Default is FALSE.

## Value

A list with two elements: - \`SIF\`: A data frame representing the
filtered meta network. - \`ATT\`: A data frame representing the updated
moon scores.

## Examples

``` r
# Example usage (requires appropriate data structures for meta_network and moon_scores)
# result <- get_moon_scoring_network(upstream_node, meta_network, moon_scores)
```
