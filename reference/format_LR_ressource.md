# Format Ligand-Receptor Resource

This function formats a ligand-receptor resource by creating a gene set
with source-target pairs, converting it to a long format, and adding
default values for 'mor' and 'likelihood'.

## Usage

``` r
format_LR_ressource(ligrec_ressource)
```

## Arguments

- ligrec_ressource:

  A data frame representing the ligand-receptor resource with columns
  for source and target gene symbols.

## Value

A data frame containing the formatted ligand-receptor gene set with
columns:

- gene:

  The gene symbol from the ligand-receptor pairs.

- set:

  The set identifier combining source and target gene symbols.

- mor:

  Default value set to 1 for all entries.

- likelihood:

  Default value set to 1 for all entries.

## Examples

``` r
# Create a sample ligand-receptor resource
ligrec_ressource <- data.frame(source_genesymbol = c("L1", "L2"),
                               target_genesymbol = c("R1", "R2"))

# Format the ligand-receptor resource
formatted_geneset <- format_LR_ressource(ligrec_ressource)
```
