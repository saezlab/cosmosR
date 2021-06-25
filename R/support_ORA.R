#Copyright (C) 2021  Caroline Lohoff, Aurelien Dugourd

#based on "support_enrichment.R" from
#https://github.com/saezlab/transcriptutorial/blob/master/scripts/06_analysis_CARNIVAL_results.md
#Copyright (C) 2020  Aurelien Dugourd, Alberto Valdeolivas, Rosa Hernansaiz-Ballesteros
#Contact : aurelien.dugourd@bioquant.uni-heidelberg.de

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.


#'\code{gmt_to_csv}
#'
#'This function is designed to convert a gmt file (gene set file from MSigDB)
#'into a two column dataframe where the first column corresponds to omic features
#'(genes) and the second column to associated terms (pathway the gene belongs to).
#'One gene can belong to various pathways.
#'
#'@param gmtfile a full path name of the gmt file to be converted
#'@param outfile an optional output file name. If none is provided, the function
#' will simply return a dataframe. If outfile is provided with a full length
#' path name file, the dataframe will be written as a csv file to the path provided.
#'@return a two column dataframe where the first column corresponds to omic
#'features and the second column to associated terms (pathways).
#'@import GSEABase

gmt_to_csv <- function(gmtfile, fast = T)
{
  if(fast)
  {
    genesets = GSEABase::getGmt(con = gmtfile)
    genesets = unlist(genesets)
    
    gene_to_term =plyr::ldply(genesets,function(geneset){
      temp <- geneIds(geneset)
      temp2 <- setName(geneset)
      temp3 <- as.data.frame(cbind(temp,rep(temp2,length(temp))))
      
    },.progress = plyr::progress_text())
    names(gene_to_term) <- c("gene","term")
    return(gene_to_term[complete.cases(gene_to_term),])
  }
  else
  {
    genesets = getGmt(con = gmtfile)
    genesets = unlist(genesets)
    
    gene_to_term <- data.frame(NA,NA)
    names(gene_to_term) <- c("gene","term")
    for (geneset in genesets)
    {
      temp <- geneIds(geneset)
      temp2 <- setName(geneset)
      temp3 <- as.data.frame(cbind(temp,rep(temp2,length(temp))))
      names(temp3) <- c("gene","term")
      gene_to_term <- rbind(gene_to_term,temp3)
    }
    
    return(gene_to_term[complete.cases(gene_to_term),])
  }
}



#'\code{extract_COSMOS_nodes}
#'
#'Function to extract the nodes that appear in the COSMOS output network and
#'the background genes (all genes present in the prior knowledge network)
#'
#'@param CosmosResults  COSMOS output
#'@return List with 2 objects: the success and the background genes
#'@import magrittr

extract_COSMOS_nodes <- function(CosmosResults){
  
  # Extract all nodes from COSMOS's solution network
  CosmosNetwork <- 
    as.data.frame(CosmosResults$weightedSIF, stringsAsFactors = FALSE)
  colnames(CosmosNetwork) <- c("source", "sign", "target", "weight")
  
  # Remove all metabolites
  CosmosNetwork <- CosmosNetwork[!grepl("Metab", CosmosNetwork$source), ]
  CosmosNetwork <- CosmosNetwork[!grepl("Metab", CosmosNetwork$target), ]
  
  ## We define the set of unique nodes (genes) in the solution network
  sucesses <- unique(c(gsub("_.*","",CosmosNetwork$source), 
                       gsub("_.*","",CosmosNetwork$target)))
  
  # Extract all unique genes from PKN as background nodes 
  Cosmos_attributes <- as.data.frame(CosmosResults$nodesAttributes, 
                                     stringsAsFactors = FALSE)
  Cosmos_attributes <- Cosmos_attributes[!grepl("Metab", Cosmos_attributes$Nodes), ]
  bg <- unique(gsub("_.*","",Cosmos_attributes$Nodes)) 

  return(list(sucesses = sucesses, bg = bg))
}



#'\code{plot_pathways}
#'
#'Function to display the most significantly over-expressed pathways sorted by
#'their source (e.g. Reactome, Biocarta, etc.)
#'
#'@param sigPathwaysDf  Data frame containing all significant pathways as a result of the ORA
#'@param cutoff number used as a filter criteria for p-value
#'@return plot of the 5 most significant pathways per source
#'@import ggplot2
#'@import magrittr
#'@import scales
#'@export

plot_pathways <- function(sigPathwaysDf, cutoff){

  # Prepare data for plotting
  PathwaysSelect <- sigPathwaysDf %>%
    dplyr::select(pathway, `p-value`, `Adjusted p-value`) %>%
    dplyr::rename(pvalue = `p-value`, AdjPvalue = `Adjusted p-value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
  
  PathwaysSelect <- data.frame(t(apply(PathwaysSelect, 1, function(r){
    aux = unlist(strsplit(sub("_",";", r["pathway"]), ";" ))
    r["pathway"] = gsub("_", " ", aux[2])
    return(c(r, "source" = aux[1]))
  })))
  
  colnames(PathwaysSelect) <- c("pathway", "pvalue", "AdjPvalue", "source")
  PathwaysSelect$AdjPvalue = as.numeric(PathwaysSelect$AdjPvalue)
  PathwaysSelect$pvalue = as.numeric(PathwaysSelect$pvalue)
  
  PathwaysSelect$pathway <- tolower(PathwaysSelect$pathway)
  PathwaysSelect$pathway <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",
                                 PathwaysSelect$pathway, perl = TRUE)
  
  ggdata = PathwaysSelect %>% 
    dplyr::filter(pvalue <= cutoff) %>% 
    dplyr::group_by(source) %>% 
    dplyr::arrange(pvalue) %>%
    dplyr::slice(1:5)
  
  # Visualize top results
  plot <- ggplot(ggdata, aes(y = reorder(pathway, pvalue),
                             x = -log10(pvalue),
                         color = source)) + 
                 geom_bar(stat = "identity") +
                 facet_grid(source ~ ., scales="free_y") +
                 ggtitle("Most Significant Pathways By Source") +
                 scale_x_continuous(
                   expand = c(0.01, 0.01),
                   limits = c(0, ceiling(max(-log10(PathwaysSelect$pvalue)))),
                   breaks = seq(floor(min(-log10(PathwaysSelect$pvalue))),
                                ceiling(max(-log10(PathwaysSelect$pvalue))), 3),
                   labels = math_format(10^-.x)
                 ) +
                 annotation_logticks(sides = "bt") +
                 theme_bw() +
                 theme(axis.title = element_text(face = "bold", size = 12),
                       axis.text.y = element_text(size = 10)) +
                 ylab("")
  
#  ggsave("significant_pathways.png", plot = pathway_plot,
#         path = "results/", scale = 1, dpi = 300, limitsize = TRUE)
  
  return(plot)
}



#'\code{top_pathways}
#'
#'Function to plot the most significantly over-expressed pathways (max. 20)
#'(regardless of source, e.g. Reactome, Biocarta)
#'
#'@param sigPathwaysDf  Data frame containing all significant pathways as a result of the ORA
#'@return Bar plot of the 20 most significant pathways by p-value
#'@import ggplot2
#'@export

top_pathways <- function(sigPathwaysDf){

  # Filtering out top 20 pathways by p-value
  sigPathwaysDf <- sigPathwaysDf[order(sigPathwaysDf$'p-value'),]
  top_pathway <- sigPathwaysDf[1:20,c(1,2,3)] 
  top_pathway <- top_pathway[order(top_pathway$'p-value', decreasing = TRUE),]

  # Write pathways with capital letter and without underscore
  top_pathway$pathway <- gsub("_"," ",top_pathway$pathway)
  top_pathway$pathway <- tolower(top_pathway$pathway)
  top_pathway$pathway <- gsub("(^|\\p{P})(.)", "\\1\\U\\2",    
                              top_pathway$pathway, perl = TRUE)

  # Calculate -log10(p-value)
  top_pathway$`p-value` <- -log10(top_pathway$`p-value`)
  top_pathway$pathway <- factor(top_pathway$pathway, levels = top_pathway$pathway)
  names(top_pathway)[2] <- "-log10(p-value)"

  plot <- ggplot(top_pathway, aes(x = pathway, y = `-log10(p-value)`,
                                   fill = `-log10(p-value)`)) + 
    geom_bar(stat = "identity") + 
    coord_flip() + 
    theme_minimal() + 
    ggtitle("ORA TOP 20 PATHWAYS") +
    scale_fill_gradient(low="grey", high="darkred")
  
#  ggsave("top_pathways.png", plot = plot,
#         path = "results/", scale = 1, dpi = 300, limitsize = TRUE)
  
  return(plot)
}
