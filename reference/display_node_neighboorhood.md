# display_node_neighboorhood

display input and measurements within n steps of a given set of nodes

## Usage

``` r
display_node_neighboorhood(central_node, sif, att, n = 100)
```

## Arguments

- central_node:

  character or character vector; node ID(s) around which a network will
  be branched out untill meansurments and input are reached

- sif:

  df; COSMOS network solution in sif format like the first list element
  returned by the format_cosmos_res function

- att:

  df; attributes of the nodes of the COMSOS network solution like the
  second list element returned by the format_cosmos_res function

- n:

  numeric; maximum number of steps in the network to look for inputs and
  measurments

## Value

a visnetwork object

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
#> 13:33:19 10.12.2025 Carnival flavour: vanilla
#> 13:33:19 10.12.2025 Generating variables for lp problem
#> 13:33:19 10.12.2025 Done: generating variables for lp problem
#> Saving preprocessed data.
#> Done: saving parsed data: /__w/cosmosR/cosmosR/docs/reference//parsedData_t13_33_19d10_12_2025n79.RData
#> 13:33:19 10.12.2025 Generating formulation for LP problem
#> 13:33:19 10.12.2025 Done: generating formulation for LP problem.
#> Saving LP file
#> Done: Saving LP file: /__w/cosmosR/cosmosR/docs/reference//lpFile_t13_33_19d10_12_2025n79.lp
#> 13:33:19 10.12.2025 Solving LP problem
#> Parsing .lp file for lpSolve
#> Rows: 882 Columns: 1
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (1): X1
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Done: parsing .lp file for lpSolve
#> 13:33:20 10.12.2025 Done: solving LP problem.
#> 13:33:20 10.12.2025 Getting the solution matrix
#> 13:33:20 10.12.2025 Done: getting the solution matrix.
#> 13:33:20 10.12.2025 Exporting solution matrix
#> 13:33:20 10.12.2025 Done: exporting solution matrix.
#> Cleaning intermediate files
#> Done: cleaning
#> 13:33:20 10.12.2025 All tasks finished.
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
#> 13:33:20 10.12.2025 Carnival flavour: vanilla
#> 13:33:20 10.12.2025 Generating variables for lp problem
#> 13:33:20 10.12.2025 Done: generating variables for lp problem
#> Saving preprocessed data.
#> Done: saving parsed data: /__w/cosmosR/cosmosR/docs/reference//parsedData_t13_33_20d10_12_2025n77.RData
#> 13:33:20 10.12.2025 Generating formulation for LP problem
#> 13:33:20 10.12.2025 Done: generating formulation for LP problem.
#> Saving LP file
#> Done: Saving LP file: /__w/cosmosR/cosmosR/docs/reference//lpFile_t13_33_20d10_12_2025n77.lp
#> 13:33:20 10.12.2025 Solving LP problem
#> Parsing .lp file for lpSolve
#> Rows: 882 Columns: 1
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (1): X1
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Done: parsing .lp file for lpSolve
#> 13:33:20 10.12.2025 Done: solving LP problem.
#> 13:33:20 10.12.2025 Getting the solution matrix
#> 13:33:20 10.12.2025 Done: getting the solution matrix.
#> 13:33:20 10.12.2025 Exporting solution matrix
#> 13:33:20 10.12.2025 Done: exporting solution matrix.
#> Cleaning intermediate files
#> Done: cleaning
#> 13:33:20 10.12.2025 All tasks finished.
#> 
#> --- End of the CARNIVAL pipeline --- 
test_result_for <- format_COSMOS_res(test_result_for)
network_plot <- display_node_neighboorhood(central_node = 'MYC',
sif = test_result_for[[1]],
att = test_result_for[[2]],
n = 7)
#> Warning: At vendor/cigraph/src/paths/unweighted.c:444 : Couldn't reach some vertices.
network_plot

{"x":{"nodes":{"id":["Enzyme1338__AKR1A1","Enzyme7338__LDHB_LDHA_reverse","Enzyme9049__SLC2A1_reverse","FGFR1","LDHA","MAPK14","MYC","Metab__D-Glucitol_c","Metab__D-Glucose_c","Metab__L-Lactic acid_c","SLC2A1"],"NodeType":["","","","","","","P","M","","M",""],"ZeroAct":[0,0,0,0,0,0,0,0,0,0,0],"UpAct":[1,1,1,1,1,1,1,1,1,1,1],"DownAct":[0,0,0,0,0,0,0,0,0,0,0],"AvgAct":[1,1,1,1,1,1,1,1,1,1,1],"measured":[0,0,0,0,0,0,1,1,0,1,0],"Activity":[1,1,1,1,1,1,1,1,1,1,1],"label":["Enzyme1338__AKR1A1","Enzyme7338__LDHB_LDHA_reverse","Enzyme9049__SLC2A1_reverse","FGFR1","LDHA","MAPK14","MYC","Metab__D-Glucitol_c","Metab__D-Glucose_c","Metab__L-Lactic acid_c","SLC2A1"],"color":["green","green","green","green","green","green","green","green","green","green","green"],"shape":["dot","dot","dot","dot","dot","dot","dot","dot","dot","dot","dot"],"shadow":[false,false,false,false,false,false,true,true,false,true,false]},"edges":{"from":["MAPK14","Metab__D-Glucitol_c","FGFR1","Enzyme1338__AKR1A1","Enzyme7338__LDHB_LDHA_reverse","Enzyme9049__SLC2A1_reverse","LDHA","MAPK14","MYC","Metab__D-Glucose_c","SLC2A1"],"to":["MYC","MAPK14","LDHA","Metab__D-Glucitol_c","Metab__L-Lactic acid_c","Metab__D-Glucose_c","Enzyme7338__LDHB_LDHA_reverse","FGFR1","SLC2A1","Enzyme1338__AKR1A1","Enzyme9049__SLC2A1_reverse"],"sign":[1,1,1,1,1,1,1,-1,1,1,1],"weigth":[1,1,1,1,1,1,1,0,1,1,1],"color":["grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey"],"arrows.to.type":["arrow","arrow","arrow","arrow","arrow","arrow","arrow","circle","arrow","arrow","arrow"],"enabled":[true,true,true,true,true,true,true,true,true,true,true],"scaleFactor":[1,1,1,1,1,1,1,1,1,1,1]},"nodesToDataframe":true,"edgesToDataframe":true,"options":{"width":"100%","height":"100%","nodes":{"shape":"dot"},"manipulation":{"enabled":false}},"groups":null,"width":1600,"height":1600,"idselection":{"enabled":true,"style":"width: 200px; height: 26px;\n                                              background: #f8f8f8;\n                                              color: darkblue;\n                                              border:none;\n                                              outline:none;","useLabels":true,"main":"Select by id"},"byselection":{"enabled":false,"style":"width: 150px; height: 26px","multiple":false,"hideColor":"rgba(200,200,200,0.5)","highlight":false},"main":null,"submain":null,"footer":null,"background":"rgba(0, 0, 0, 0)","highlight":{"enabled":true,"hoverNearest":false,"degree":1,"algorithm":"all","hideColor":"rgba(200,200,200,0.5)","labelOnly":true},"collapse":{"enabled":false,"fit":false,"resetHighlight":true,"clusterOptions":null,"keepCoord":true,"labelSuffix":"(cluster)"}},"evals":[],"jsHooks":[]}
```
