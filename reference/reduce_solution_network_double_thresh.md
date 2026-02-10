# reduce_solution_network_double_thresh

Extracts a subnetwork using a two-tier absolute-score threshold, then
restricts to only paths connecting upstream_input seeds to level 0
nodes, with recursive consistency filtering to ensure edge signs match
node activities.

## Usage

``` r
reduce_solution_network_double_thresh(
  decoupleRnival_res,
  meta_network,
  primary_thresh,
  secondary_thresh,
  upstream_input,
  RNA_input = NULL
)
```

## Arguments

- decoupleRnival_res:

  A data.frame with columns \`source\`, numeric \`score\`, and integer
  \`level\`.

- meta_network:

  A data.frame with columns \`source\`, \`target\`, and \`interaction\`.

- primary_thresh:

  Numeric. Absolute score cutoff for primary node selection.

- secondary_thresh:

  Numeric. Absolute score cutoff for secondary node restriction.

- upstream_input:

  A named numeric vector of upstream seed nodes.

- RNA_input:

  Optional named numeric vector of differential expression values;
  merged into ATT.

## Value

A list with: - SIF: data.frame of filtered edges
(\`source\`,\`target\`,\`interaction\`). - ATT: data.frame of node
attributes (\`source\`,\`score\`,\`level\`,\`type\`,\`RNA_input\`).

## Examples

``` r
dec_res <- data.frame(
  source = paste0("G", 1:6),
  score  = c(2.5, 1.2, 0.8, -2.2, 1.5, -0.5),
  level  = c(0, 0, 1, 0, 1, 1)
)
meta_net <- data.frame(
  source      = c("G1","G1","G2","G3","G4","G5"),
  target      = c("G2","G3","G4","G5","G2","G6"),
  interaction = c(1, -1, 1, 1, -1, 1)
)
upstream_input <- c(G1 = 1, G4 = -1)
RNA_input <- c(G1 = 0.5, G2 = -0.2, G4 = 1.1)
dbl_net <- reduce_solution_network_double_thresh(
  decoupleRnival_res = dec_res,
  meta_network       = meta_net,
  primary_thresh     = 2,
  secondary_thresh   = 1,
  upstream_input     = upstream_input,
  RNA_input          = RNA_input
)
print(dbl_net$SIF)
#>   source target interaction
#> 1     G1     G2           1
#> 5     G4     G2          -1
print(dbl_net$ATT)
#>   source score level           type RNA_input
#> 1     G1   2.5     0 upstream_input       0.5
#> 2     G2   1.2     0         level0      -0.2
#> 3     G4  -2.2     0 upstream_input       1.1
```
