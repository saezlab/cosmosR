---
title: "COSMOS-tutorial"
author: "A. Dugourd, A. Gabor and K. Zirngibl"
date: "11/10/2020"
output:
  html_document: 
    toc: true
    toc_depth: 2
    keep_md: true
vignette: >
  %\VignetteIndexEntry{cosmosR tutorial}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
---
  


## Installation

```r
# install from bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("cosmosR")

# install the development version from GitHub
# install.packages("remotes")
remotes::install_github("saezlab/cosmosR")

# Install an older stable version of CARNIVAL (1.3.0)
remotes::install_github("saezlab/CARNIVAL@b3a84c6ba9706547caca02644566d75ee621f568")
```

## Introduction

COSMOS (Causal Oriented Search of Multi-Omic Space) is a method that integrates phosphoproteomics, transcriptomics, and metabolomics data sets. COSMOS leverages extensive prior knowledge of signaling pathways, metabolic networks, and gene regulation with computational methods to estimate activities of transcription factors and kinases as well as network-level causal reasoning. This pipeline can provide mechanistic explanations for experimental observations across multiple omic data sets.

![data_intro_figure](../inst/figures/intro_data.png)

First, we load the package


```r
library(cosmosR)
```

## Tutorial section: signaling to metabolism

In this part, we can set up the options for the CARNIVAL run, such as timelimit and min gap tolerance.

The user should provide a path to its CPLEX/cbc executable.

You can check the CARNIVAL_options variable to see all possible options that can be adjusted

In this example, we will use the built-in solver lpSolve. User should be aware that lpSolve should ONLY be used for TESTS. To obtain meaningful results, best solver is cplex, or cbc if not possible.


```r
CARNIVAL_options <- cosmosR::default_CARNIVAL_options()
# CARNIVAL_options$solverPath <- "C:/Program Files/CPLEX_solver/cplex/bin/x64_win64/cplex.exe"
# CARNIVAL_options$solver <- "cplex" #or cbc
# CARNIVAL_options$solverPath <- "~/Documents/cplex"
# CARNIVAL_options$solver <- "cplex" #or cbc
CARNIVAL_options$solver <- "lpSolve" #or cbc
CARNIVAL_options$timelimit <- 3600
CARNIVAL_options$mipGAP <- 0.05
CARNIVAL_options$threads <- 2
```

In the next section, we prepare the input to run cosmosR.
The signaling inputs are the result of footprint based TF and kinase activity estimation.
For more info on TF activity estimation from transcriptomic data, see:https://github.com/saezlab/transcriptutorial (Especially chapter 4)

Here we use of toy PKN, to see the full meta PKN, you can load it with data(meta_network).

The metabolites in the prior knowledge network are identified as XMetab__PUBCHEMid___compartment____ or XMetab__BIGGid___compartment____
(for example “XMetab__6804___m____”). The compartment code is the BIGG model standard (r, c, e, x, m, l, n, g). 
Thus we will first need to map whatever identifier for metabolite the data has to the one of the network.
Genes are identified as XENTREZid (in the signaling part of network) or XGene####__ENTREZid (in the reaction network part of network).

The maximum network depth will define the maximum number of step downstream of kinase/TF COSMOS will look for deregulated metabolites. Good first guess for max depth could be around 6 (here it is 15 for the toy dataset)

The differential experession data is used to filter out wrong TF-target interactions in this context after a pre-optimisation.

The list of genes in the differential expression data will also be used as a reference to define which genes are expressed or not (all genes in the diff_expression_data are considered expressed, and genes that are no in diff_expression_data are removed from the network)


```r
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
```

```
## [1] "COSMOS: all 3 signaling nodes from data were found in the meta PKN"
## [1] "COSMOS: all 3 metabolic nodes from data were found in the meta PKN"
## [1] "COSMOS: 4660 of the 15919 genes in expression data were found as transcription factor target"
## [1] "COSMOS: 4660 of the 5312 transcription factor targets were found in expression data"
## [1] "COSMOS: removing unexpressed nodes from PKN..."
## [1] "COSMOS: 0 interactions removed"
## [1] "COSMOS: removing nodes that are not reachable from inputs within 15 steps"
## [1] "COSMOS: 86 from  115 interactions are removed from the PKN"
## [1] "COSMOS: 1 input/measured nodes are not in PKN any more: XMetab__439155___c____ and 0 more."
## [1] "COSMOS: removing nodes that are not observable by measurements within 15 steps"
## [1] "COSMOS: 10 from  29 interactions are removed from the PKN"
## [1] "COSMOS: 1 input/measured nodes are not in PKN any more: X1445 and 0 more."
## [1] "COSMOS:  0 interactions are removed from the PKN based on consistency check between TF activity and gene expression"
## [1] "lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run COSMSO, we recommend using cplex, or cbc solvers."
## [1] "lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run COSMSO, we recommend using cplex, or cbc solvers."
## [1] "COSMOS:  0 interactions are removed from the PKN based on consistency check between TF activity and gene expression"
## [1] "COSMOS: all 2 signaling nodes from data were found in the meta PKN"
## [1] "COSMOS: all 2 metabolic nodes from data were found in the meta PKN"
## [1] "COSMOS: 4660 of the 15919 genes in expression data were found as transcription factor target"
## [1] "COSMOS: 4660 of the 5312 transcription factor targets were found in expression data"
```

In this part, we can set up the options for the actual run, such as timelimit and min gap tolerance.

The running time should be much higher here than in pre-optimisation. You cna increase the number of threads to use if you have many available CPUs.


```r
CARNIVAL_options$timelimit <- 14400
CARNIVAL_options$mipGAP <- 0.05
CARNIVAL_options$threads <- 2
```

This is where cosmosR run.


```r
test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
                                                      CARNIVAL_options = CARNIVAL_options)
```

```
## [1] "lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run COSMSO, we recommend using cplex, or cbc solvers."
## [1] "lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run COSMSO, we recommend using cplex, or cbc solvers."
```

Finally, we process the results of the first cosmosR run, to translate gene names and metabolites name.


```r
data(metabolite_to_pubchem)
data(omnipath_ptm)
formated_result_for <- format_COSMOS_res(test_result_for,
                                     metab_mapping = metabolite_to_pubchem,
                     measured_nodes = unique(c(names(toy_metabolic_input),
                                               names(toy_signaling_input))),
                                     omnipath_ptm = omnipath_ptm)
```

## Tutorial section: metabolism to signaling 

Before we run the metabolism to signaling part, we need to prepare again the inputs. 

```r
CARNIVAL_options$timelimit <- 3600
CARNIVAL_options$mipGAP <- 0.05
CARNIVAL_options$threads <- 2
```

Now that the correct time is set up for the pre-optimisation run, we can prepare the inputs.


```r
test_back <- preprocess_COSMOS_metabolism_to_signaling(meta_network = toy_network,
                                        signaling_data = toy_signaling_input,
                                        metabolic_data = toy_metabolic_input,
                                                       diff_expression_data = toy_RNA,
                                                       maximum_network_depth = 15,
                                                       remove_unexpressed_nodes = FALSE,
                                                       CARNIVAL_options = CARNIVAL_options)
```

```
## [1] "COSMOS: all 3 signaling nodes from data were found in the meta PKN"
## [1] "COSMOS: all 3 metabolic nodes from data were found in the meta PKN"
## [1] "COSMOS: 4660 of the 15919 genes in expression data were found as transcription factor target"
## [1] "COSMOS: 4660 of the 5312 transcription factor targets were found in expression data"
## [1] "COSMOS: removing nodes that are not reachable from inputs within 15 steps"
## [1] "COSMOS: 105 from  115 interactions are removed from the PKN"
## [1] "COSMOS: 2 input/measured nodes are not in PKN any more: X4790, X5062 and 0 more."
## [1] "COSMOS: 2 input/measured nodes are not in PKN any more: XMetab__65359___c____, XMetab__107738___m____ and 0 more."
## [1] "COSMOS: removing nodes that are not observable by measurements within 15 steps"
## [1] "COSMOS: 5 from  10 interactions are removed from the PKN"
## [1] "COSMOS:  0 interactions are removed from the PKN based on consistency check between TF activity and gene expression"
## [1] "lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run COSMSO, we recommend using cplex, or cbc solvers."
## [1] "lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run COSMSO, we recommend using cplex, or cbc solvers."
## [1] "COSMOS:  0 interactions are removed from the PKN based on consistency check between TF activity and gene expression"
## [1] "COSMOS: all 1 signaling nodes from data were found in the meta PKN"
## [1] "COSMOS: all 1 metabolic nodes from data were found in the meta PKN"
## [1] "COSMOS: 4660 of the 15919 genes in expression data were found as transcription factor target"
## [1] "COSMOS: 4660 of the 5312 transcription factor targets were found in expression data"
```

Then we can run cosmosR to connect metabolism to signaling. The running time here usually needs to be longer, as this problem seems to be harder to solve for CPLEX.


```
## [1] "lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run COSMSO, we recommend using cplex, or cbc solvers."
## [1] "lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run COSMSO, we recommend using cplex, or cbc solvers."
```

Finally we can format the result of the backward run as well (same as for forward run)


```r
formated_result_back <- format_COSMOS_res(test_result_back,
                                      metab_mapping = metabolite_to_pubchem,
                      measured_nodes = unique(c(names(toy_metabolic_input),
                                                names(toy_signaling_input))),
                                      omnipath_ptm = omnipath_ptm)
```

## Tutorial section: Merge forward and backward networks and visualise network

Here we simply take the union of forward and backward runs to create a full network solution lopping between signaling, gene-regulation and metabolism. Since there is an overlap between the result network of forward and backward run, you may optionally want  to check if there are any node sign that are incoherent in the overlap between the two solutions.


```r
full_sif <- as.data.frame(rbind(formated_result_for[[1]], formated_result_back[[1]]))
full_attributes <- as.data.frame(rbind(formated_result_for[[2]], formated_result_back[[2]]))

full_sif <- unique(full_sif)
full_attributes <- unique(full_attributes)
```

This function will generate a dynamic network plot centered on a given node of the network solution, and connecting it to measured nodes in the given range (here 7 steps).


```r
network_plot <- display_node_neighboorhood(central_node = 'NFKB1', 
                                           sif = full_sif, 
                                           att = full_attributes, 
                                           n = 7)

network_plot
```

```{=html}
<div id="htmlwidget-5af8a9c093b9a27277f4" style="width:1600px;height:1600px;" class="visNetwork html-widget"></div>
<script type="application/json" data-for="htmlwidget-5af8a9c093b9a27277f4">{"x":{"nodes":{"id":["Enzyme0__CRAT_c3","Enzyme0__HADHB","HADHB","Metab__92753___m____","Metab__O-Propanoylcarnitine___m____","NFKB1","TBK1","TNFAIP3"],"ZeroAct":[0,0,0,0,0,0,0,0],"UpAct":[1,1,1,1,1,1,1,1],"DownAct":[0,0,0,0,0,0,0,0],"AvgAct":[1,1,1,1,1,1,1,1],"NodeType":["","","","","T","S","",""],"measured":[0,0,0,0,1,1,0,0],"type":["metab_enzyme","metab_enzyme","Kinase","metabolite","metabolite","TF","Kinase","metab_enzyme"],"Activity":[1,1,1,1,1,1,1,1],"label":["Enzyme0__CRAT_c3","Enzyme0__HADHB","HADHB","Metab__92753___m____","Metab__O-Propanoylcarnitine___m____","NFKB1","TBK1","TNFAIP3"],"color":["green","green","green","green","green","green","green","green"],"shape":["square","square","triangle","dot","dot","diamond","triangle","square"],"shadow":[false,false,false,false,true,true,false,false]},"edges":{"from":["TBK1","Enzyme0__HADHB","Metab__92753___m____","Enzyme0__CRAT_c3","HADHB","NFKB1","TNFAIP3"],"to":["HADHB","Metab__92753___m____","Enzyme0__CRAT_c3","Metab__O-Propanoylcarnitine___m____","Enzyme0__HADHB","TNFAIP3","TBK1"],"sign":[1,1,1,1,1,1,1],"weigth":[1,1,1,1,1,1,1],"color":["grey","grey","grey","grey","grey","grey","grey"],"arrows.to.type":["arrow","arrow","arrow","arrow","arrow","arrow","arrow"],"enabled":[true,true,true,true,true,true,true],"scaleFactor":[1,1,1,1,1,1,1]},"nodesToDataframe":true,"edgesToDataframe":true,"options":{"width":"100%","height":"100%","nodes":{"shape":"dot"},"manipulation":{"enabled":false}},"groups":null,"width":1600,"height":1600,"idselection":{"enabled":true,"style":"width: 200px; height: 26px;\n                                              background: #f8f8f8;\n                                              color: darkblue;\n                                              border:none;\n                                              outline:none;","useLabels":true,"main":"Select by id"},"byselection":{"enabled":false,"style":"width: 150px; height: 26px","multiple":false,"hideColor":"rgba(200,200,200,0.5)","highlight":false},"main":null,"submain":null,"footer":null,"background":"rgba(0, 0, 0, 0)","highlight":{"enabled":true,"hoverNearest":false,"degree":1,"algorithm":"all","hideColor":"rgba(200,200,200,0.5)","labelOnly":true},"collapse":{"enabled":false,"fit":false,"resetHighlight":true,"clusterOptions":null,"keepCoord":true,"labelSuffix":"(cluster)"}},"evals":[],"jsHooks":[]}</script>
```
Here is how this network can be intepreted (this is purely illustrative, as some of those interaction may be incorectly signed because lpsolve can only use positive interactions): 

This network represents the flow of activities that can connect NFKB1 up-regulation with o-propanoylcarnitine accumulation
accumulation. Here, NFKB1 can upregulate the expression of TNFAIP3, which in turn activate TBK1 (potentially through ubiquitin editing mechanism). The activation of TBK1 can lead to the phosphorilation and activation of HADHB kinase, leading
to the increased production of propionyl-Coa, which in turn can be converted to o-propanoylcarnitine by the CRAT enzyme.

It is important to understand that each of this links is hypothetical. The come from a larger pool of 
potential molecular interactions present in multiple online databases and compiled in omnipath, STITCH 
and recon metabolic network. They exist in the literature and are interactions that are known to
potentially exists in other experimental contexts. Thus, COSMOS compile all those potential interactions 
together and proposes a coherent set that can explain the data at hand.

Those links should however be considered only as potential mechanistic connections, and will need to be further 
confirmed experimentally.

## Tutorial section: Over Representation Analysis

Often it is useful to perform an Over Representation Analysis (ORA) on the
resulting nodes of a COSMOS network as a first analysis step to get a more 
functional interpretation on the modeled signaling cascade. A common way to this
is to test whether the selected genes (nodes) in the COSMOS solution network
show statistically significant differences in comparison to the prior-knowledge
network (PKN). 

The differentially expressed genes give information about the cellular
processes that are deregulated and if the proportions in various pathways are
SIGNIFICANTLY different from what is expected.In this way the significant 
differences between two biological conditions (e.g. cancer vs. normal tissue, or
treatment vs. untreated cells) can be shown.

Algorithms that perform an ORA are implemented in other R packages like piano or
decoupleR. In addition to a gene set collection these algorithms require two 
different lists as inputs:
- nodes in the COSMOS solution network which relate back to the input data
  (e.g. transcriptomics, proteomics, metabolomics, fluxomics, or perturbations) 
- all nodes (kinases, transcription factors, metabolites) in the prior-knowledge
  network (which are used as the background in our analysis)
  
In this section we will show how to obtain these two lists from a formated 
COSMOS result object.


```r
sif_forward = formated_result_for[[1]]
att_forward = formated_result_for[[2]]
nodes_ORA = extract_nodes_for_ORA(
    sif = sif_forward, 
    att = att_forward)
```

Now this forground and background set can be used with any ORA anaylsis. 


```r
sessionInfo()
```

```
## R version 4.1.0 (2021-05-18)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] cosmosR_1.1.2
## 
## loaded via a namespace (and not attached):
##   [1] segmented_1.3-4        Category_2.58.0        bitops_1.0-7          
##   [4] bit64_4.0.5            filelock_1.0.2         doParallel_1.0.16     
##   [7] httr_1.4.2             GenomeInfoDb_1.28.0    tools_4.1.0           
##  [10] bslib_0.2.5.1          utf8_1.2.1             R6_2.5.0              
##  [13] KernSmooth_2.23-20     DBI_1.1.1              BiocGenerics_0.38.0   
##  [16] colorspace_2.0-1       tidyselect_1.1.1       curl_4.3.1            
##  [19] bit_4.0.4              compiler_4.1.0         cli_3.0.0             
##  [22] graph_1.70.0           Biobase_2.52.0         CARNIVAL_1.3.0        
##  [25] sass_0.4.0             scales_1.1.1           readr_2.1.0           
##  [28] genefilter_1.74.0      proxy_0.4-26           RBGL_1.68.0           
##  [31] rappdirs_0.3.3         stringr_1.4.0          digest_0.6.27         
##  [34] mixtools_1.2.0         rmarkdown_2.10         XVector_0.32.0        
##  [37] pkgconfig_2.0.3        htmltools_0.5.1.1      bcellViper_1.28.0     
##  [40] dbplyr_2.1.1           fastmap_1.1.0          htmlwidgets_1.5.3     
##  [43] rlang_0.4.11           rstudioapi_0.13        RSQLite_2.2.7         
##  [46] visNetwork_2.0.9       jquerylib_0.1.4        generics_0.1.0        
##  [49] jsonlite_1.7.2         vroom_1.5.6            viper_1.26.0          
##  [52] dplyr_1.0.6            RCurl_1.98-1.3         magrittr_2.0.1        
##  [55] GenomeInfoDbData_1.2.6 Matrix_1.3-4           Rcpp_1.0.7            
##  [58] munsell_0.5.0          S4Vectors_0.30.0       fansi_0.5.0           
##  [61] lifecycle_1.0.0        stringi_1.6.2          yaml_2.2.1            
##  [64] UniProt.ws_2.32.0      MASS_7.3-54            zlibbioc_1.38.0       
##  [67] org.Hs.eg.db_3.13.0    BiocFileCache_2.0.0    grid_4.1.0            
##  [70] blob_1.2.1             parallel_4.1.0         crayon_1.4.1          
##  [73] lattice_0.20-44        Biostrings_2.60.1      splines_4.1.0         
##  [76] annotate_1.70.0        hms_1.1.0              KEGGREST_1.32.0       
##  [79] knitr_1.33             dorothea_1.4.1         pillar_1.6.1          
##  [82] igraph_1.2.6           lpSolve_5.6.15         codetools_0.2-18      
##  [85] stats4_4.1.0           XML_3.99-0.6           glue_1.4.2            
##  [88] evaluate_0.14          png_0.1-7              vctrs_0.3.8           
##  [91] tzdb_0.1.2             foreach_1.5.1          gtable_0.3.0          
##  [94] purrr_0.3.4            kernlab_0.9-29         assertthat_0.2.1      
##  [97] cachem_1.0.5           ggplot2_3.3.3          xfun_0.23             
## [100] xtable_1.8-4           e1071_1.7-7            class_7.3-19          
## [103] survival_3.2-11        tibble_3.1.2           iterators_1.0.13      
## [106] AnnotationDbi_1.54.1   memoise_2.0.0          IRanges_2.26.0        
## [109] ellipsis_0.3.2         GSEABase_1.54.0
```
