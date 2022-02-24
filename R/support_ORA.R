# Copyright (C) 2021  Caroline Lohoff, Aurelien Dugourd

# based on "support_enrichment.R" from
# https://github.com/saezlab/transcriptutorial/blob/master/scripts/06_analysis_CARNIVAL_results.md
# Copyright (C) 2020  Aurelien Dugourd, Alberto Valdeolivas, Rosa Hernansaiz-Ballesteros
# Contact : aurelien.dugourd@bioquant.uni-heidelberg.de

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' Convert gmt file to data frame
#' 
#' This function is designed to convert a gmt file (gene set file from MSigDB)
#' into a two column data frame where the first column corresponds to omic
#' features (genes) and the second column to associated terms (pathway the gene
#' belongs to). One gene can belong to several pathways.
#' 
#' @param gmtfile a full path name of the gmt file to be converted
#' @return a two column data frame where the first column corresponds to omic
#' features and the second column to associated terms (pathways).

gmt_to_dataframe <- function(gmtfile)
{
    genesets = GSEABase::getGmt(con = gmtfile)
    genesets = unlist(genesets)
    
    gene_to_term =plyr::ldply(genesets,function(geneset){
    temp <- GSEABase::geneIds(geneset)
    temp2 <- GSEABase::setName(geneset)
    temp3 <- as.data.frame(cbind(temp,rep(temp2,length(temp))))
    
    },.progress = plyr::progress_text())
    names(gene_to_term) <- c("gene","term")
    return(gene_to_term[stats::complete.cases(gene_to_term),])
}


#' Extract COSMOS nodes for ORA analysis
#' 
#' Function to extract the nodes that appear in the COSMOS output network and
#' the background genes (all genes present in the prior knowledge network)
#' 
#' @param sif df; COSMOS network solution in sif format like the first list
#' element returned by the format_cosmos_res function
#' @param att df; attributes of the nodes of the COMSOS network solution like 
#' the second list element returned by the format_cosmos_res function
#' @return List with 2 objects: the success and the background genes
#' @export
#' @examples
#' CARNIVAL_options <- cosmosR::default_CARNIVAL_options()
#' CARNIVAL_options$solver <- "lpSolve"
#' data(toy_network)
#' data(toy_signaling_input)
#' data(toy_metabolic_input)
#' data(toy_RNA)
#' test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = toy_network,
#' signaling_data = toy_signaling_input,
#' metabolic_data = toy_metabolic_input,
#' diff_expression_data = toy_RNA,
#' maximum_network_depth = 15,
#' remove_unexpressed_nodes = TRUE,
#' CARNIVAL_options = CARNIVAL_options
#' )
#' test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
#' CARNIVAL_options = CARNIVAL_options)
#' test_result_for <- format_COSMOS_res(test_result_for)
#' extreacted_nodes <- extract_nodes_for_ORA(
#' sif = test_result_for[[1]],
#' att = test_result_for[[2]]
#' )

extract_nodes_for_ORA <- function(sif, att){
  
  # Extract all nodes from COSMOS's solution network
  CosmosNetwork <- 
    as.data.frame(sif, stringsAsFactors = FALSE)
  colnames(CosmosNetwork) <- c("source", "sign", "target", "weight")
  
  # Remove all metabolites
  CosmosNetwork <- CosmosNetwork[!grepl("Metab", CosmosNetwork$source), ]
  CosmosNetwork <- CosmosNetwork[!grepl("Metab", CosmosNetwork$target), ]
  
  ## We define the set of unique nodes (genes) in the solution network
  sucesses <- unique(c(gsub("_.*","",CosmosNetwork$source), 
                       gsub("_.*","",CosmosNetwork$target)))
  
  # Extract all unique genes from PKN as background nodes 
  Cosmos_attributes <- as.data.frame(att, 
                                     stringsAsFactors = FALSE)
  Cosmos_attributes <- Cosmos_attributes[!grepl("Metab", Cosmos_attributes$Nodes), ]
  bg <- unique(gsub("_.*","",Cosmos_attributes$Nodes)) 

  return(list(sucesses = sucesses, bg = bg))
}
