# run COSMOS metabolism to signaling

Runs COSMOS from metabolism to signaling. This function uses CARNIVAL to
find a subset of the prior knowledge network based on optimization that
(1) includes the most measured and input nodes and (2) which is in
agreement with the data. Use
[`preprocess_COSMOS_metabolism_to_signaling`](preprocess_COSMOS_metabolism_to_signaling.md)
to prepare the the inputs, measurements and the prior knowledge network.

## Usage

``` r
run_COSMOS_metabolism_to_signaling(
  data,
  CARNIVAL_options = default_CARNIVAL_options("lpSolve")
)
```

## Arguments

- data:

  [`cosmos_data`](cosmos_data.md) object. Use the
  [`preprocess_COSMOS_metabolism_to_signaling`](preprocess_COSMOS_metabolism_to_signaling.md)
  function to create an instance.

- CARNIVAL_options:

  List that controls the options of CARNIVAL. See details in
  [`default_CARNIVAL_options`](default_CARNIVAL_options.md).

## Value

List with the following elements:

- `weightedSIF`:

  The averaged networks found by optimization in a format of a Simple
  Interaction network, i.e. each row codes an edge

- `N_networks`:

  Number of solutions found by the optimization

- `nodesAttributes`:

  Estimated node properties

- `individual_networks`:

  List of optimial networks found

- `individual_networks_node_attributes`:

  Node activity in each network

## See also

[`preprocess_COSMOS_metabolism_to_signaling`](preprocess_COSMOS_metabolism_to_signaling.md),
[`runCARNIVAL`](https://rdrr.io/pkg/CARNIVAL/man/runCARNIVAL.html),
[`cosmos_data`](cosmos_data.md)

## Examples

``` r
data(toy_network)
data(toy_signaling_input)
data(toy_metabolic_input)
data(toy_RNA)
test_back <- preprocess_COSMOS_metabolism_to_signaling(meta_network = toy_network,
                        signaling_data = toy_signaling_input,
                        metabolic_data = toy_metabolic_input,
                        diff_expression_data = toy_RNA,
                        maximum_network_depth = 15,
                        remove_unexpressed_nodes = TRUE,
                        CARNIVAL_options = default_CARNIVAL_options("lpSolve"))
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
#> 13:33:27 10.12.2025 Carnival flavour: vanilla
#> 13:33:27 10.12.2025 Generating variables for lp problem
#> 13:33:27 10.12.2025 Done: generating variables for lp problem
#> Saving preprocessed data.
#> Done: saving parsed data: /__w/cosmosR/cosmosR/docs/reference//parsedData_t13_33_27d10_12_2025n17.RData
#> 13:33:27 10.12.2025 Generating formulation for LP problem
#> 13:33:27 10.12.2025 Done: generating formulation for LP problem.
#> Saving LP file
#> Done: Saving LP file: /__w/cosmosR/cosmosR/docs/reference//lpFile_t13_33_27d10_12_2025n17.lp
#> 13:33:27 10.12.2025 Solving LP problem
#> Parsing .lp file for lpSolve
#> Rows: 842 Columns: 1
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (1): X1
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Done: parsing .lp file for lpSolve
#> 13:33:27 10.12.2025 Done: solving LP problem.
#> 13:33:27 10.12.2025 Getting the solution matrix
#> 13:33:27 10.12.2025 Done: getting the solution matrix.
#> 13:33:27 10.12.2025 Exporting solution matrix
#> 13:33:27 10.12.2025 Done: exporting solution matrix.
#> Cleaning intermediate files
#> Done: cleaning
#> 13:33:27 10.12.2025 All tasks finished.
#> 
#> --- End of the CARNIVAL pipeline --- 
#> [1] "COSMOS:  0 interactions are removed from the PKN based on consistency check between TF activity and gene expression"
#> [1] "COSMOS: all 3 signaling nodes from data were found in the meta PKN"
#> [1] "COSMOS: all 1 metabolic nodes from data were found in the meta PKN"
#> [1] "COSMOS: 2975 of the 9300 genes in expression data were found as transcription factor target"
#> [1] "COSMOS: 2975 of the 5321 transcription factor targets were found in expression data"

test_result_back <- run_COSMOS_metabolism_to_signaling(data = test_back,
                       CARNIVAL_options = default_CARNIVAL_options("lpSolve"))
#> [1] "lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run COSMSO, we recommend using cplex, or cbc solvers."
#> Input nodes should have values from {-1, 0, 1}. We discretize your input with sign().
#> [1] "lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run COSMSO, we recommend using cplex, or cbc solvers."
#> --- Start of the CARNIVAL pipeline ---
#> 13:33:28 10.12.2025 Carnival flavour: vanilla
#> 13:33:28 10.12.2025 Generating variables for lp problem
#> 13:33:28 10.12.2025 Done: generating variables for lp problem
#> Saving preprocessed data.
#> Done: saving parsed data: /__w/cosmosR/cosmosR/docs/reference//parsedData_t13_33_28d10_12_2025n52.RData
#> 13:33:28 10.12.2025 Generating formulation for LP problem
#> 13:33:28 10.12.2025 Done: generating formulation for LP problem.
#> Saving LP file
#> Done: Saving LP file: /__w/cosmosR/cosmosR/docs/reference//lpFile_t13_33_28d10_12_2025n52.lp
#> 13:33:28 10.12.2025 Solving LP problem
#> Parsing .lp file for lpSolve
#> Rows: 842 Columns: 1
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (1): X1
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Done: parsing .lp file for lpSolve
#> 13:33:28 10.12.2025 Done: solving LP problem.
#> 13:33:28 10.12.2025 Getting the solution matrix
#> 13:33:28 10.12.2025 Done: getting the solution matrix.
#> 13:33:28 10.12.2025 Exporting solution matrix
#> 13:33:28 10.12.2025 Done: exporting solution matrix.
#> Cleaning intermediate files
#> Done: cleaning
#> 13:33:28 10.12.2025 All tasks finished.
#> 
#> --- End of the CARNIVAL pipeline --- 
```
