# Preprocess COSMOS Inputs For Metabolism to Signaling

Runs checks on the input data and simplifies the prior knowledge
network. Simplification includes the removal of (1) nodes that are not
reachable from signaling nodes and (2) interactions between
transcription factors and target genes if the target gene does not
respond or the response is contradictory with the change in the
transcription factor activity. Optionally, further TF activities are
estimated via network optimization via CARNIVAL and the interactions
between TF and genes are filtered again.

## Usage

``` r
preprocess_COSMOS_metabolism_to_signaling(
  meta_network = meta_network,
  tf_regulon = load_tf_regulon_dorothea(),
  signaling_data,
  metabolic_data,
  diff_expression_data = NULL,
  diff_exp_threshold = 1,
  maximum_network_depth = 8,
  expressed_genes = NULL,
  remove_unexpressed_nodes = TRUE,
  filter_tf_gene_interaction_by_optimization = TRUE,
  CARNIVAL_options = default_CARNIVAL_options("lpSolve")
)
```

## Arguments

- meta_network:

  prior knowledge network (PKN). A PKN released with COSMOS and derived
  from Omnipath, STITCHdb and Recon3D can be used. See details on the
  data [`meta_network`](meta_network.md).

- tf_regulon:

  collection of transcription factor - target interactions. A default
  collection from dorothea can be obtained by the
  [`load_tf_regulon_dorothea`](load_tf_regulon_dorothea.md) function.

- signaling_data:

  numerical vector, where names are signaling nodes in the PKN and
  values are from {1, 0, -1}. Continuous data will be discretized using
  the [`sign`](https://rdrr.io/r/base/sign.html) function.

- metabolic_data:

  numerical vector, where names are metabolic nodes in the PKN and
  values are continuous values that represents log2 fold change or
  t-values from a differential analysis. These values are compared to
  the simulation results (simulated nodes can take value -1, 0 or 1)

- diff_expression_data:

  (optional) numerical vector that represents the results of a
  differential gene expression analysis. Names are gene names using gene
  symbole and values are log fold change or t-values. We use the
  “`diff_exp_threshold`” parameter to decide which genes changed
  significantly. Genes with NA values are considered none expressed and
  they will be removed from the TF-gene expression interactions.

- diff_exp_threshold:

  threshold parameter (default 1) used to binarize the values of
  “`diff_expression_data`”.

- maximum_network_depth:

  integer \> 0 (default: 8). Nodes that are further than
  “`maximum_network_depth`” steps from the signaling nodes on the
  directed graph of the PKN are considered non-reachable and are
  removed.

- expressed_genes:

  character vector. Names of nodes that are expressed. By default we
  consider all the nodes that appear in `diff_expression_data` with a
  numeric value (i.e. nodes with NA are removed)

- remove_unexpressed_nodes:

  if TRUE (default) removes nodes from the PKN that are not expressed,
  see input “`expressed_genes`”.

- filter_tf_gene_interaction_by_optimization:

  (default:TRUE), if TRUE then runs a network optimization that
  estimates TF activity not included in the inputs and checks the
  consistency between the estimated activity and change in gene
  expression. Removes interactions where TF and gene expression are
  inconsistent

- CARNIVAL_options:

  list that controls the options of CARNIVAL. See details in
  [`default_CARNIVAL_options`](default_CARNIVAL_options.md).

## Value

cosmos_data object with the following fields:

- `meta_network`:

  Filtered PKN

- `tf_regulon`:

  TF - target regulatory network

- `signaling_data_bin`:

  Binarised signaling data

- `metabolic_data`:

  Metabolomics data

- `diff_expression_data_bin`:

  Binarized gene expression data

- `optimized_network`:

  Initial optimized network if
  `filter_tf_gene_interaction_by_optimization is TRUE`

## See also

[`meta_network`](meta_network.md) for meta PKN,
[`load_tf_regulon_dorothea`](load_tf_regulon_dorothea.md) for tf
regulon,
[`runCARNIVAL`](https://rdrr.io/pkg/CARNIVAL/man/runCARNIVAL.html).

## Examples

``` r
data(toy_network)
data(toy_signaling_input)
data(toy_metabolic_input)
data(toy_RNA)
test_back <- preprocess_COSMOS_metabolism_to_signaling(
meta_network = toy_network,
signaling_data = toy_signaling_input,
metabolic_data = toy_metabolic_input,
diff_expression_data = toy_RNA,
maximum_network_depth = 15,
remove_unexpressed_nodes = TRUE,
CARNIVAL_options = default_CARNIVAL_options("lpSolve")
)
#> [1] "COSMOS: all 3 signaling nodes from data were found in the meta PKN"
#> [1] "COSMOS: all 2 metabolic nodes from data were found in the meta PKN"
#> [1] "COSMOS: 2975 of the 9300 genes in expression data were found as transcription factor target"
#> [1] "COSMOS: 2975 of the 5321 transcription factor targets were found in expression data"
#> [1] "COSMOS: removing unexpressed nodes from PKN..."
#> [1] "COSMOS: 0 interactions removed"
#> [1] "COSMOS: removing nodes that are not reachable from inputs within 15 steps"
#> [1] "COSMOS: 0 from  101 interactions are removed from the PKN"
#> [1] "COSMOS: removing nodes that are not observable by measurements within 15 steps"
#> [1] "COSMOS: 54 from  101 interactions are removed from the PKN"
#> [1] "COSMOS: 1 input/measured nodes are not in PKN any more: Metab__HMDB0000190_c and 0 more."
#> [1] "COSMOS:  0 interactions are removed from the PKN based on consistency check between TF activity and gene expression"
#> [1] "lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run COSMSO, we recommend using cplex, or cbc solvers."
#> [1] "lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run COSMSO, we recommend using cplex, or cbc solvers."
#> --- Start of the CARNIVAL pipeline ---
#> 10:39:49 10.02.2026 Carnival flavour: vanilla
#> 10:39:49 10.02.2026 Generating variables for lp problem
#> 10:39:49 10.02.2026 Done: generating variables for lp problem
#> Saving preprocessed data.
#> Done: saving parsed data: /__w/cosmosR/cosmosR/docs/reference//parsedData_t10_39_49d10_02_2026n32.RData
#> 10:39:49 10.02.2026 Generating formulation for LP problem
#> 10:39:49 10.02.2026 Done: generating formulation for LP problem.
#> Saving LP file
#> Done: Saving LP file: /__w/cosmosR/cosmosR/docs/reference//lpFile_t10_39_49d10_02_2026n32.lp
#> 10:39:49 10.02.2026 Solving LP problem
#> Parsing .lp file for lpSolve
#> Rows: 842 Columns: 1
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (1): X1
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Done: parsing .lp file for lpSolve
#> 10:39:49 10.02.2026 Done: solving LP problem.
#> 10:39:49 10.02.2026 Getting the solution matrix
#> 10:39:49 10.02.2026 Done: getting the solution matrix.
#> 10:39:49 10.02.2026 Exporting solution matrix
#> 10:39:49 10.02.2026 Done: exporting solution matrix.
#> Cleaning intermediate files
#> Done: cleaning
#> 10:39:49 10.02.2026 All tasks finished.
#> 
#> --- End of the CARNIVAL pipeline --- 
#> [1] "COSMOS:  0 interactions are removed from the PKN based on consistency check between TF activity and gene expression"
#> [1] "COSMOS: all 3 signaling nodes from data were found in the meta PKN"
#> [1] "COSMOS: all 1 metabolic nodes from data were found in the meta PKN"
#> [1] "COSMOS: 2975 of the 9300 genes in expression data were found as transcription factor target"
#> [1] "COSMOS: 2975 of the 5321 transcription factor targets were found in expression data"
```
