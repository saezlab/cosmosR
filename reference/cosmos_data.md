# Create Cosmos Data

An S3 class that combines the required data into a comprehensive list.
Use the
[`preprocess_COSMOS_signaling_to_metabolism`](preprocess_COSMOS_signaling_to_metabolism.md)
or
[`preprocess_COSMOS_metabolism_to_signaling`](preprocess_COSMOS_metabolism_to_signaling.md)
to create an instance.

## Usage

``` r
cosmos_data(
  meta_network,
  tf_regulon = NULL,
  signaling_data,
  metabolic_data,
  expression_data,
  verbose = TRUE
)
```

## Arguments

- meta_network:

  Prior knowledge network (PKN). By default COSMOS use a PKN derived
  from Omnipath, STITCHdb and Recon3D. See details on the data
  [`meta_network`](meta_network.md).

- tf_regulon:

  Collection of transcription factor - target interactions. A default
  collection from dorothea can be obtained by the
  [`load_tf_regulon_dorothea`](load_tf_regulon_dorothea.md) function.

- signaling_data:

  Numerical vector, where names are signaling nodes in the PKN and
  values are from {1, 0, -1}. Continuous data will be discretized using
  the [`sign`](https://rdrr.io/r/base/sign.html) function.

- metabolic_data:

  Numerical vector, where names are metabolic nodes in the PKN and
  values are continuous values that represents log2 fold change or
  t-values from a differential analysis. These values are compared to
  the simulation results (simulated nodes can take value -1, 0 or 1).

- expression_data:

  Numerical vector that represents the results of a differential gene
  expression analysis. Names are gene names using EntrezID starting with
  an X and values are log fold change or t-values. Genes with NA values
  are considered none expressed and they will be removed from the
  TF-gene expression interactions.

- verbose:

  (default: TRUE) Reports details about the `cosmos_data` object.

## Value

`cosmos data` class instance.
