# filter_incohrent_TF_target

Filters incoherent target genes from a regulatory network based on a
decoupling analysis of upstream and downstream gene expression.

## Usage

``` r
filter_incohrent_TF_target(
  decouplRnival_res,
  TF_reg_net,
  meta_network,
  RNA_input
)
```

## Arguments

- decouplRnival_res:

  A data frame resulting from the decoupleRnival function.

- TF_reg_net:

  A data frame containing prior knowledge of transcription factor (TF)
  regulatory interactions.

- meta_network:

  A network data frame containing signed directed prior knowledge of
  molecular interactions.

- RNA_input:

  A named vector containing differential gene expression data.

## Value

A network data frame containing the genes that are not incoherently
regulated by TFs.

## Examples

``` r
# Example input data
upstream_input <- c("A" = 1, "B" = -1, "C" = 0.5)
downstream_input <- c("D" = 2, "E" = -1.5)
meta_network <- data.frame(
  source = c("A", "A", "B", "C", "C", "D", "E"),
  target = c("B", "D", "D", "E", "D", "B", "A"),
  interaction = c(-1, 1, -1, 1, -1, -1, 1)
)
RNA_input <- c("A" = 1, "B" = -1, "C" = 5, "D" = -0.7, "E" = -0.3)

TF_reg_net <- data.frame(
source = c("B"),
target = c("D"),
mor = c(-1)
)

# Run the decoupleRnival function to get the upstream influence scores
upstream_scores <- decoupleRnival(upstream_input, downstream_input, meta_network, n_layers = 2, n_perm = 100)
#> [1] "Warning, this function is deprecated and will no longer receive futur support. Please use the 'moon' function instead"

filtered_network <- filter_incohrent_TF_target(upstream_scores, TF_reg_net, meta_network, RNA_input)

print(filtered_network)
#>   source target interaction
#> 1      A      B          -1
#> 2      A      D           1
#> 4      C      E           1
#> 5      C      D          -1
#> 6      D      B          -1
#> 7      E      A           1
```
