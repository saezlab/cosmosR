# add metabolic compartment and metab\_\_ prefix to metabolite IDs

This function adds metabolic compartments to the metabolic identifiers
provided by the user.

## Usage

``` r
prepare_metab_inputs(metab_input, compartment_codes)
```

## Arguments

- metab_input:

  a named vector with matebolic statistics as inputs and metabolite
  identifiers as names

- compartment_codes:

  a character vector, the desired compartment codes to be added.
  Possible values are "r", "c", "e", "x", "m", "l", "n" and "g"

## Value

a named vector with the compartment code and prefixed added to the names
