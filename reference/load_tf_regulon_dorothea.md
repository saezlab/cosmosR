# load transcription factor regulon

load the transcription factors from `DOROTHEA` package and converts gene
symbols to EntrezID using org.Hs.eg.db

## Usage

``` r
load_tf_regulon_dorothea(confidence = c("A", "B", "C"))
```

## Arguments

- confidence:

  strong vector (by default: c("A","B","C")). Subset of {A, B, C, D, E}.
  See the \`dorothea\` for the meaning of confidence levels. package for
  further details.

## Value

returns a PKN of a form of a data table. Each row is an interaction.
Columns names are:

\- \`tf\` transcription factor - \`confidence\` class of confidence -
\`target\` target gene - \`sign\` indicates if interaction is up (1) or
down-regulation (-1).

## Examples

``` r
load_tf_regulon_dorothea()
#> # A tibble: 13,223 × 3
#>    tf    target    sign
#>    <chr> <chr>    <dbl>
#>  1 AHR   CYP1A1       1
#>  2 AHR   CYP1A2       1
#>  3 AHR   CYP1B1       1
#>  4 AHR   FOS          1
#>  5 AHR   MYC          1
#>  6 AHR   UGT1A6       1
#>  7 AHR   ASAP1        1
#>  8 AHR   ERG          1
#>  9 AHR   VGLL4        1
#> 10 AHR   ARHGAP15     1
#> # ℹ 13,213 more rows
```
