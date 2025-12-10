# Convert gmt file to data frame

This function is designed to convert a gmt file (gene set file from
MSigDB) into a two column data frame where the first column corresponds
to omic features (genes) and the second column to associated terms
(pathway the gene belongs to). One gene can belong to several pathways.

## Usage

``` r
gmt_to_dataframe(gmtfile)
```

## Arguments

- gmtfile:

  a full path name of the gmt file to be converted

## Value

a two column data frame where the first column corresponds to omic
features and the second column to associated terms (pathways).
