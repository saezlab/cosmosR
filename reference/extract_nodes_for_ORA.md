# Extract COSMOS nodes for ORA analysis

Function to extract the nodes that appear in the COSMOS output network
and the background genes (all genes present in the prior knowledge
network)

## Usage

``` r
extract_nodes_for_ORA(sif, att)
```

## Arguments

- sif:

  df; COSMOS network solution in sif format like the first list element
  returned by the format_cosmos_res function

- att:

  df; attributes of the nodes of the COMSOS network solution like the
  second list element returned by the format_cosmos_res function

## Value

List with 2 objects: the success and the background genes

## Examples

``` r
CARNIVAL_options <- cosmosR::default_CARNIVAL_options("lpSolve")
data(toy_network)
data(toy_signaling_input)
data(toy_metabolic_input)
data(toy_RNA)
test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = toy_network,
signaling_data = toy_signaling_input,
metabolic_data = toy_metabolic_input,
diff_expression_data = toy_RNA,
maximum_network_depth = 15,
remove_unexpressed_nodes = TRUE,
CARNIVAL_options = CARNIVAL_options
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
#> [1] "COSMOS: 52 from  101 interactions are removed from the PKN"
#> [1] "COSMOS: 2 input/measured nodes are not in PKN any more: USF1, SRF and 0 more."
#> [1] "COSMOS:  0 interactions are removed from the PKN based on consistency check between TF activity and gene expression"
#> [1] "lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run COSMSO, we recommend using cplex, or cbc solvers."
#> [1] "lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run COSMSO, we recommend using cplex, or cbc solvers."
#> --- Start of the CARNIVAL pipeline ---
#> 10:39:44 10.02.2026 Carnival flavour: vanilla
#> 10:39:44 10.02.2026 Generating variables for lp problem
#> 10:39:44 10.02.2026 Done: generating variables for lp problem
#> Saving preprocessed data.
#> Done: saving parsed data: /__w/cosmosR/cosmosR/docs/reference//parsedData_t10_39_44d10_02_2026n55.RData
#> 10:39:44 10.02.2026 Generating formulation for LP problem
#> 10:39:44 10.02.2026 Done: generating formulation for LP problem.
#> Saving LP file
#> Done: Saving LP file: /__w/cosmosR/cosmosR/docs/reference//lpFile_t10_39_44d10_02_2026n55.lp
#> 10:39:44 10.02.2026 Solving LP problem
#> Parsing .lp file for lpSolve
#> Rows: 882 Columns: 1
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (1): X1
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Done: parsing .lp file for lpSolve
#> 10:39:45 10.02.2026 Done: solving LP problem.
#> 10:39:45 10.02.2026 Getting the solution matrix
#> 10:39:45 10.02.2026 Done: getting the solution matrix.
#> 10:39:45 10.02.2026 Exporting solution matrix
#> 10:39:45 10.02.2026 Done: exporting solution matrix.
#> Cleaning intermediate files
#> Done: cleaning
#> 10:39:45 10.02.2026 All tasks finished.
#> 
#> --- End of the CARNIVAL pipeline --- 
#> [1] "COSMOS:  0 interactions are removed from the PKN based on consistency check between TF activity and gene expression"
#> [1] "COSMOS: all 1 signaling nodes from data were found in the meta PKN"
#> [1] "COSMOS: all 2 metabolic nodes from data were found in the meta PKN"
#> [1] "COSMOS: 2975 of the 9300 genes in expression data were found as transcription factor target"
#> [1] "COSMOS: 2975 of the 5321 transcription factor targets were found in expression data"
test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
CARNIVAL_options = CARNIVAL_options)
#> [1] "lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run COSMSO, we recommend using cplex, or cbc solvers."
#> Input nodes should have values from {-1, 0, 1}. We discretize your input with sign().
#> [1] "lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run COSMSO, we recommend using cplex, or cbc solvers."
#> --- Start of the CARNIVAL pipeline ---
#> 10:39:45 10.02.2026 Carnival flavour: vanilla
#> 10:39:45 10.02.2026 Generating variables for lp problem
#> 10:39:45 10.02.2026 Done: generating variables for lp problem
#> Saving preprocessed data.
#> Done: saving parsed data: /__w/cosmosR/cosmosR/docs/reference//parsedData_t10_39_45d10_02_2026n43.RData
#> 10:39:45 10.02.2026 Generating formulation for LP problem
#> 10:39:45 10.02.2026 Done: generating formulation for LP problem.
#> Saving LP file
#> Done: Saving LP file: /__w/cosmosR/cosmosR/docs/reference//lpFile_t10_39_45d10_02_2026n43.lp
#> 10:39:45 10.02.2026 Solving LP problem
#> Parsing .lp file for lpSolve
#> Rows: 882 Columns: 1
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (1): X1
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Done: parsing .lp file for lpSolve
#> 10:39:45 10.02.2026 Done: solving LP problem.
#> 10:39:45 10.02.2026 Getting the solution matrix
#> 10:39:45 10.02.2026 Done: getting the solution matrix.
#> 10:39:45 10.02.2026 Exporting solution matrix
#> 10:39:46 10.02.2026 Done: exporting solution matrix.
#> Cleaning intermediate files
#> Done: cleaning
#> 10:39:46 10.02.2026 All tasks finished.
#> 
#> --- End of the CARNIVAL pipeline --- 
test_result_for <- format_COSMOS_res(test_result_for)
extreacted_nodes <- extract_nodes_for_ORA(
sif = test_result_for[[1]],
att = test_result_for[[2]]
)
```
