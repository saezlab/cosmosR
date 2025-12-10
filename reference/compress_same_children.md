# Compress Network by Merging Nodes with Identical Children

This function compresses a network by merging nodes that have the same
children. The input network is represented as a data frame with three
columns: source, target, and sign of interaction. The function returns a
list containing the compressed network, node signatures, and duplicated
signatures.

## Usage

``` r
compress_same_children(df, sig_input, metab_input)
```

## Arguments

- df:

  A data frame representing the network with three columns: source,
  target, and sign of interaction.

- sig_input:

  A list of input node signatures to be considered for the merging
  process.

- metab_input:

  A list of input metabolic signatures to be considered for the merging
  process.

## Value

A list containing the following elements:

- compressed_network:

  A data frame representing the compressed network.

- node_signatures:

  A list of signatures of nodes in the network after the merging
  process.

- duplicated_signatures:

  A list of duplicated signatures in the network after the merging
  process.

## Examples

``` r
# Create a sample network
df <- data.frame(source = c("A", "A", "B", "B"),
                 target = c("C", "D", "C", "D"),
                 sign_of_interaction = c(1, 1, 1, 1))

# Define input node and metabolic signatures
sig_input <- list()
metab_input <- list()

# Compress the network
result <- compress_same_children(df, sig_input, metab_input)
compressed_network <- result$compressed_network
```
