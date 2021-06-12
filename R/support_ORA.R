##################### Supporting functions for ORA COSMOS ######################

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

library(piano)
library(biomaRt)  # interface to databases, e.g. Ensembl
library(parallel)
library(GSEABase)
library(snowfall)
library(readr)
library(ggplot2)
library(scales) 
library(magrittr)
library("org.Hs.eg.db")  # Genome wide annotation for human, based on mapping
                         # using Entrez Gene identifiers


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
#'features and the second column to associated terms.

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


#'\code{geneMapping}
#'
#'Function to map all genes from a list that was extracted from a network to
#'Ensembl or ENTREZID identifiers.
#'
#'@param GenList  List of genes extracted from a network (genes as identifiers)
#'@return List of genes (genes as names)

geneMapping <- function(GenList){
  
  prots <- GenList[!grepl("XMetab",GenList)] # filter out metabolites
  prots <- gsub("^X","",prots)
  prots <- gsub("Gene[0-9]+__","",prots)
  prots <- gsub("_reverse","",prots)
  prots <- gsub("EXCHANGE.*","",prots)
  prots <- unique(prots)
  prots <- unlist(sapply(prots, function(x){unlist(strsplit(x,split = "_"))}))
  
  gene_mapping <- "else"
  
  if(gene_mapping == "ensembl")
  {
    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    
    G_list <- getBM(filters = "entrezgene_id",
                    attributes = c('hgnc_symbol','entrezgene_id', "description"),
                    values = metabs, mart = ensembl)
    gene_mapping <- G_list[,1]
    names(gene_mapping) <- G_list[,2]
    
  } else
  {
    entrezid <- prots
    gene_mapping <- mapIds(org.Hs.eg.db, entrezid, 'SYMBOL', 'ENTREZID')
    gene_mapping <- unlist(gene_mapping)
    gene_mapping <- gene_mapping[!is.na(gene_mapping)]
  }
  
  output_list <- sapply(prots, function(x, translation_dictionary){
    return(translation_dictionary[x])
  },translation_dictionary = gene_mapping)
  names(output_list) <- 1:length(output_list)
  
  output_list <- unique(output_list)
  return(output_list)
}



#'\code{extractCOSMOSnodes}
#'
#'Function to extract the nodes that appear in the COSMOS output network and
#'the background genes (all genes present in the prior knowledge network)
#'
#'@param CosmosResults  COSMOS output.
#'@return List with 2 objects: the success and the background genes.

extractCOSMOSnodes <- function(CosmosResults){
  
  ## 1. Extract all nodes from COSMOS's solution network
  # extract network dataframe from Cosmos's output list
  CosmosNetwork <- 
    as.data.frame(CosmosResults$weightedSIF, stringsAsFactors = FALSE)
  # renaming of column names in df
  colnames(CosmosNetwork) <- c("source", "sign", "target", "Weight")
  
  ## We define the set of nodes interesting for our condition
  #  (if node is inhibited or activated is in this case not important, we only
  #   want to get all nodes from Cosmos's solution network as a list)
  # unique: we get every node exactly one time
  sucesses <- unique(c(gsub("_.*","",CosmosNetwork$source), 
                       gsub("_.*","",CosmosNetwork$target)))
  
  # Map all genes (kinase, enzyme, TF) to identifiers (code to name)
  #sucesses_sol_net <- unique(c(CosmosNetwork$source,CosmosNetwork$target))
  #sucesses_genes = geneMapping(sucesses_sol_net)
  
  # Map all metabolites
  # sucesses_metabs = metaboliteMapping(sucesses_sol_net)
  # sucesses = append(sucesses_genes, sucesses_metabs)
  
  
  # Now we extract for each node the detailed activation or downregulation
  # score and if the node is a TF or a perturbed node
  CosmosAttributes <- as.data.frame(CosmosResults$nodesAttributes, 
                                      stringsAsFactors = FALSE)
  ## 2. Extract all nodes from PKN
  # We define the background as all the genes in our prior knowledge network.
  bg <- unique(gsub("_.*","",CosmosAttributes$Nodes)) 
  
  # Map all genes (kinase, enzyme, TF) to identifiers (code to name) 
  #bg_sol_net <- unique(CosmosAttributes$Node) 
  #bg_genes = geneMapping(bg_sol_net)
  
  # Map all metabolites
  # bg_metabs = metaboliteMapping(bg_sol_net)
  # bg = append(bg_genes, bg_metabs)
  
  return(list(sucesses = sucesses, bg = bg))
  # If mapping of identifiers was required
  #return(list(sucesses = sucesses_genes, bg = bg_genes))
}



#'\code{plotPathways}
#'
#'Function to plot the pathways which are significantly overexpressed
#'
#'@param sigPathwaysDf  Dataframe containing all significant pathways as a result of the ORA
#'@return Bar plot of the most significant pathways

plotPathways <- function(sigPathwaysDf){

  #Prepare data for plotting
  PathwaysSelect <- sigPathwaysDf %>%
    dplyr::select(pathway, `p-value`, `Adjusted p-value`) %>%
    dplyr::filter(`Adjusted p-value` <= 0.001) %>%  
    dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))

  PathwaysSelect <- data.frame(t(apply(PathwaysSelect, 1, function(r){
    aux = unlist(strsplit( sub("_",";", r["pathway"]), ";" ))
    r["pathway"] = gsub("_", " ", aux[2])
    return(c(r, "source" = aux[1]))
  })))

  colnames(PathwaysSelect) = c("pathway", "pvalue", "AdjPvalu", "source")
  PathwaysSelect$AdjPvalu = as.numeric(PathwaysSelect$AdjPvalu)

  ggdata = PathwaysSelect %>% 
    dplyr::filter(AdjPvalu <= 0.05) %>% 
    dplyr::group_by(source) %>% 
    dplyr::arrange(AdjPvalu) %>%
    dplyr::slice(1:5)

  #Visualize top results
  pathway_plot <- ggplot(ggdata, aes(y = reorder(pathway, AdjPvalu),
                                     x = -log10(AdjPvalu)),
                         color = source) + 
                        geom_bar(stat = "identity") +
                        facet_grid(source ~ ., scales="free_y") +
                        ggtitle("Most Significant Pathways") +
                        scale_x_continuous(
                          expand = c(0.01, 0.01),
                          limits = c(0, ceiling(max(-log10(PathwaysSelect$AdjPvalu)))),
                          breaks = seq(floor(min(-log10(PathwaysSelect$AdjPvalu))),
                                   ceiling(max(-log10(PathwaysSelect$AdjPvalu))), 3),
                          labels = math_format(10^-.x)
                        ) +
                        annotation_logticks(sides = "bt") +
                        theme_bw() +
                        theme(axis.title = element_text(face = "bold", size = 12),
                              axis.text.y = element_text(size = 6)) +
                        ylab("")
  ggsave("significant_pathways.png", plot = pathway_plot,
         path = "results/", scale = 1, dpi = 300, limitsize = TRUE)
  
  return(pathway_plot)
}


#'\code{plotPathways2}
#'
#'Function to plot the pathways which are significantly overexpressed
#'
#'@param sigPathwaysDf  Dataframe containing all significant pathways as a result of the ORA
#'@return Bar plot of the 20 most significant pathways

plotPathways2 <- function(sigPathwaysDf){

  # filtering out top 20 pathways by Adjusted p-value
  sigPathwaysDf <- sigPathwaysDf[order(sigPathwaysDf$'Adjusted p-value'),]
  top_hallmark <- sigPathwaysDf[1:20,c(7,1,2)] 
  top_hallmark <- top_hallmark[order(top_hallmark$'p-value', decreasing = TRUE),]
  
  # write pathways with capital letter and without underscore & HALLMARK
  top_hallmark$pathway <- gsub("HALLMARK","",top_hallmark$pathway)
  top_hallmark$pathway <- tolower(top_hallmark$pathway)
  top_hallmark$pathway <- gsub("(^|\\p{P})(.)",
                                   "\\1\\U\\2",    
                               top_hallmark$pathway, perl = TRUE)
  top_hallmark$pathway <- gsub("_"," ",top_hallmark$pathway)
  
  top_hallmark$`p-value` <- -log10(top_hallmark$`p-value`)
  top_hallmark$pathway <- factor(top_hallmark$pathway, levels = top_hallmark$pathway)
  names(top_hallmark)[3] <- "-log10(p-value)"

  write_csv(top_hallmark,'results/top_pathways.csv')

  plot <- ggplot(top_hallmark, aes(x = pathway, y = `-log10(p-value)`,
                                   fill = `-log10(p-value)`)) + 
    geom_bar(stat = "identity") + 
    coord_flip() + 
    theme_minimal() + 
    ggtitle("ORA TOP 20 PATHWAYS") +
    scale_fill_gradient(low="grey", high="darkred")
  
  ggsave("top_pathways.png", plot = plot,
         path = "results/", scale = 1, dpi = 300, limitsize = TRUE)
  
  return(plot)
}
