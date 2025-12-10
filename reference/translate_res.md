# translate_res

formats the network with readable names

## Usage

``` r
translate_res(SIF, ATT, HMDB_mapper_vec = NULL)
```

## Arguments

- SIF:

  result SIF of decoupleRnival pipeline

- ATT:

  result ATT of decoupleRnival pipeline

- HMDB_mapper_vec:

  a named vector with HMDB Ids as names and desired metabolite names as
  values.

## Value

list with network and attribute tables.

## Examples

``` r
# Create a meta network data frame
example_SIF <- data.frame(
source = c("GPX1", "Gene863__GPX1"),
target = c("Gene863__GPX1", "Metab__HMDB0003337_c"),
sign = c(1, 1)
)

example_ATT <- data.frame(
Nodes = c("GPX1", "Gene863__GPX1","Metab__HMDB0003337_c"),
sign = c(1, 1, 1)
)

example_SIF
#>          source               target sign
#> 1          GPX1        Gene863__GPX1    1
#> 2 Gene863__GPX1 Metab__HMDB0003337_c    1

data("HMDB_mapper_vec")

translated_res <- translate_res(example_SIF,example_ATT,HMDB_mapper_vec)

translated_res$SIF
#>            source                 target sign
#> 1            GPX1        Enzyme863__GPX1    1
#> 2 Enzyme863__GPX1 Metab__Oxiglutatione_c    1
```
