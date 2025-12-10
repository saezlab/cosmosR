# Decompress Moon Result

This function decompresses the results obtained from moon analysis by
incorporating node signatures and handling duplicated parents. It merges
these details with the provided meta network data and returns a
comprehensive data frame.

## Usage

``` r
decompress_moon_result(moon_res, meta_network_compressed_list, meta_network)
```

## Arguments

- moon_res:

  A data frame containing the results of a moon analysis.

- meta_network_compressed_list:

  A list containing compressed meta network details, including node
  signatures and duplicated parents.

- meta_network:

  A data frame representing the original meta network.

## Value

A data frame which merges the moon analysis results with the meta
network data, including additional details about node signatures and
handling of duplicated parents.

## Examples

``` r
# Example usage (requires appropriate data structures for moon_res, 
# meta_network_compressed_list, and meta_network)
# decompressed_result <- decompress_moon_result(moon_res, meta_network_compressed_list, meta_network)
```
