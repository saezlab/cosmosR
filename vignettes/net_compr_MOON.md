## Tutorial with NCI60 cohort

# libraries, data loading and feature selection

    # We advise to instal from github to get the latest version of the tool.
    # if (!requireNamespace("devtools", quietly = TRUE))
    #   install.packages("devtools")
    # 
    # devtools::install_github("saezlab/cosmosR")

    library(cosmosR)
    library(reshape2)
    library(readr)

    data("meta_network")

    meta_network <- meta_network_cleanup(meta_network)

    ## Warning: `summarise_each()` was deprecated in dplyr 0.7.0.
    ## ℹ Please use `across()` instead.
    ## ℹ The deprecated feature was likely used in the cosmosR package.
    ##   Please report the issue at <https://github.com/saezlab/COSMOSR/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: `funs()` was deprecated in dplyr 0.8.0.
    ## ℹ Please use a list of either functions or lambdas:
    ## 
    ## # Simple named list: list(mean = mean, median = median)
    ## 
    ## # Auto named with `tibble::lst()`: tibble::lst(mean, median)
    ## 
    ## # Using lambdas list(~ mean(., trim = .2), ~ median(., na.rm = TRUE))
    ## ℹ The deprecated feature was likely used in the cosmosR package.
    ##   Please report the issue at <https://github.com/saezlab/COSMOSR/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    load("data/cosmos/cosmos_inputs.RData")

    names(cosmos_inputs)

    ##  [1] "786-0"       "A498"        "A549/ATCC"   "ACHN"        "BT-549"     
    ##  [6] "CAKI-1"      "CCRF-CEM"    "COLO 205"    "DU-145"      "EKVX"       
    ## [11] "HCC-2998"    "HCT-116"     "HCT-15"      "HL-60(TB)"   "HOP-62"     
    ## [16] "HOP-92"      "HS 578T"     "HT29"        "IGROV1"      "K-562"      
    ## [21] "KM12"        "LOX IMVI"    "M14"         "MALME-3M"    "MCF7"       
    ## [26] "MDA-MB-435"  "MOLT-4"      "NCI-H226"    "NCI-H23"     "NCI-H322M"  
    ## [31] "NCI-H460"    "NCI-H522"    "NCI/ADR-RES" "OVCAR-3"     "OVCAR-4"    
    ## [36] "OVCAR-5"     "OVCAR-8"     "PC-3"        "RPMI-8226"   "SF-268"     
    ## [41] "SF-295"      "SF-539"      "SK-MEL-28"   "SK-MEL-5"    "SK-OV-3"    
    ## [46] "SN12C"       "SNB-19"      "SNB-75"      "SR"          "SW-620"     
    ## [51] "T-47D"       "TK-10"       "U251"        "UACC-257"    "UACC-62"    
    ## [56] "UO-31"

    cell_line <- "786-0"

    #see scripts/prepare_cosmos_inputs.R
    sig_input <- cosmos_inputs[[cell_line]]$TF_scores
    metab_input <- cosmos_inputs[[cell_line]]$metabolomic
    RNA_input <- cosmos_inputs[[cell_line]]$RNA

    #Choose which compartment to assign to the metabolic measurments
    metab_input <- prepare_metab_inputs(metab_input, c("c","m"))

    ## [1] "Adding compartment codes."

    ##Filter significant inputs
    sig_input <- sig_input[abs(sig_input) > 2]
    # metab_input <- metab_input[abs(metab_input) > 2]

# Filter inputs and prior knowledge network

    #Remove genes that are not expressed from the meta_network
    meta_network <- cosmosR:::filter_pkn_expressed_genes(names(RNA_input), meta_pkn = meta_network)

    ## [1] "COSMOS: removing unexpressed nodes from PKN..."
    ## [1] "COSMOS: 20357 interactions removed"

    #Filter inputs and prune the meta_network to only keep nodes that can be found downstream of the inputs
    #The number of step is quite flexible, 7 steps already covers most of the network

    n_steps <- 6

    # in this step we prune the network to keep only the relevant part between upstream and downstream nodes
    sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(sig_input, meta_network)

    ## [1] "COSMOS: 12 input/measured nodes are not in PKN any more: CEBPA, ESR1, FOS, FOXA1, GATA3, HNF4A and 6 more."

    meta_network <- cosmosR:::keep_controllable_neighbours(meta_network, n_steps, names(sig_input))

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 6 steps"
    ## [1] "COSMOS: 26540 from  37740 interactions are removed from the PKN"

    metab_input <- cosmosR:::filter_input_nodes_not_in_pkn(metab_input, meta_network)

    ## [1] "COSMOS: 195 input/measured nodes are not in PKN any more: Metab__HMDB0011747_c, Metab__HMDB0000755_c, Metab__HMDB0000905_c, Metab__HMDB0001191_c, Metab__HMDB0000355_c, Metab__HMDB0000479_c and 189 more."

    meta_network <- cosmosR:::keep_observable_neighbours(meta_network, n_steps, names(metab_input))

    ## [1] "COSMOS: removing nodes that are not observable by measurements within 6 steps"
    ## [1] "COSMOS: 3657 from  11200 interactions are removed from the PKN"

    sig_input <- cosmosR:::filter_input_nodes_not_in_pkn(sig_input, meta_network)

    ## [1] "COSMOS: 4 input/measured nodes are not in PKN any more: CTCF, EPAS1, ETS1, USF1 and 0 more."

    #compress the network
    meta_network_compressed_list <- compress_same_children(meta_network, sig_input = sig_input, metab_input = metab_input)

    meta_network_compressed <- meta_network_compressed_list$compressed_network

    node_signatures <- meta_network_compressed_list$node_signatures

    duplicated_parents <- meta_network_compressed_list$duplicated_signatures

    meta_network_compressed <- meta_network_cleanup(meta_network_compressed)

# run MOON ot score the and contextualise the PKN

    load("support/dorothea_reg.RData")

    meta_network_TF_to_metab <- meta_network_compressed

    before <- 1
    after <- 0
    i <- 1
    while (before != after & i < 10) {
      before <- length(meta_network_TF_to_metab[,1])
      moon_res <- moon(upstream_input = sig_input, 
                                                     downstream_input = metab_input, 
                                                     meta_network = meta_network_TF_to_metab, 
                                                     n_layers = n_steps, 
                                                     statistic = "ulm") 
      
      meta_network_TF_to_metab <- filter_incohrent_TF_target(moon_res, dorothea_reg, meta_network_TF_to_metab, RNA_input)
      after <- length(meta_network_TF_to_metab[,1])
      i <- i + 1
    }

    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5
    ## [1] 6
    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5
    ## [1] 6
    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5
    ## [1] 6

    if(i < 10)
    {
      print(paste("Converged after ",paste(i-1," iterations", sep = ""),sep = ""))
    } else
    {
      print(paste("Interupted after ",paste(i," iterations. Convergence uncertain.", sep = ""),sep = ""))
    }

    ## [1] "Converged after 3 iterations"

    #####
    write_csv(moon_res, file = paste("results/moon/",paste(cell_line, "_ATT_decouplerino_full.csv",sep = ""), sep = ""))

    source("scripts/support_decompression.R")
    moon_res <- decompress_moon_result(moon_res, meta_network_compressed_list, meta_network_TF_to_metab)

    plot(density(moon_res$score))
    abline(v = 1)
    abline(v = -1)

![](net_compr_MOON_files/figure-markdown_strict/extract%20subnetwork%20from%20scored%20MOON%20network-1.png)

    solution_network <- reduce_solution_network(decoupleRnival_res = moon_res, 
                                                meta_network = meta_network,
                                                cutoff = 0.5, 
                                                upstream_input = sig_input, 
                                                RNA_input = RNA_input, 
                                                n_steps = n_steps)

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 6 steps"
    ## [1] "COSMOS: 458 from  1747 interactions are removed from the PKN"

    SIF <- solution_network$SIF
    names(SIF)[3] <- "sign"
    ATT <- solution_network$ATT

    data("HMDB_mapper_vec")

    translated_res <- translate_res(SIF,ATT,HMDB_mapper_vec)

    SIF <- translated_res[[1]]
    ATT <- translated_res[[2]]

    write_csv(SIF, file = paste("results/",paste(cell_line, "_dec_compressed_SIF.csv",sep = ""), sep = ""))
    write_csv(ATT, file = paste("results/",paste(cell_line, "_dec_compressed_ATT.csv",sep = ""), sep = ""))

    sessionInfo()

    ## R version 4.2.0 (2022-04-22)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Monterey 12.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] readr_2.1.4    reshape2_1.4.4 cosmosR_1.5.2 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.10       highr_0.10        plyr_1.8.8        pillar_1.9.0     
    ##  [5] compiler_4.2.0    prettyunits_1.1.1 tools_4.2.0       progress_1.2.2   
    ##  [9] bit_4.0.5         digest_0.6.31     evaluate_0.20     lifecycle_1.0.3  
    ## [13] tibble_3.2.1      lattice_0.20-45   pkgconfig_2.0.3   rlang_1.1.0      
    ## [17] igraph_1.4.2      Matrix_1.5-3      cli_3.6.1         rstudioapi_0.14  
    ## [21] yaml_2.3.7        parallel_4.2.0    xfun_0.42         fastmap_1.1.1    
    ## [25] withr_2.5.0       dplyr_1.1.3       stringr_1.5.0     knitr_1.42       
    ## [29] generics_0.1.3    vctrs_0.6.1       hms_1.1.3         bit64_4.0.5      
    ## [33] grid_4.2.0        tidyselect_1.2.0  glue_1.6.2        R6_2.5.1         
    ## [37] fansi_1.0.4       parallelly_1.34.0 vroom_1.6.1       rmarkdown_2.21   
    ## [41] tzdb_0.3.0        tidyr_1.3.0       purrr_1.0.1       decoupleR_2.9.1  
    ## [45] magrittr_2.0.3    htmltools_0.5.5   utf8_1.2.3        stringi_1.7.12   
    ## [49] crayon_1.5.2
