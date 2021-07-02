#Copyright (C) 2021  Caroline Lohoff, Aurelien Dugourd
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


#'\code{create_metabolite_enzyme_df}
#'
#'Function to create a data frame where metabolites are in one column and their 
#'target or source enzyme in another column.
#'All metabolites from the cosmos solution network are used.
#'This function is a prerequisite for mapping pathways to metabolites.
#'
#'@import sjmisc

create_metabolite_enzyme_df <- function(){

  ##Create data frame by only using metabolite rows from cosmos solution network
  metabolites_nodes <- cosmos_result$nodesAttributes[cosmos_result$nodesAttributes$type == "metabolite", ]

  metabolites <- as.vector(metabolites_nodes$Nodes)
  metabolites_enzymes <- cosmos_result$weightedSIF[cosmos_result$weightedSIF$Node1 %in% metabolites, ]
  metabolites_enzymes <- rbind(metabolites_enzymes,
                             cosmos_result$weightedSIF[cosmos_result$weightedSIF$Node2 %in% metabolites, ])
  
  ##Rearrange data frame
  metabolites_enzymes[, c(2, 4)] <- NULL

  new_cols = c("enzymes", "metabolites")
  metabolites_enzymes[new_cols] <- NA

  metabolites_enzymes <- metabolites_enzymes[c("enzymes", "metabolites", "Node1", "Node2")]

  ##Copy entries from Node1 and Node2 column to metabolite or enzyme column

  for (i in 3:ncol(metabolites_enzymes)) #iterate through Node1 and Node2 column
  {
    for (j in 1:nrow(metabolites_enzymes)) #iterate through all rows
    {
      if (str_contains(metabolites_enzymes[j,i],"Metab")){
        metabolites_enzymes$metabolites[j] <- metabolites_enzymes[j,i]
      } else
      {
        metabolites_enzymes$enzymes[j] <- metabolites_enzymes[j,i]
      }
    }
  }

  metabolites_enzymes <- metabolites_enzymes[, -c(3:4)]  #delete columns Node1 & Node2
  metabolites_enzymes <- na.omit(metabolites_enzymes)    #remove rows containing NA

 return(metabolites_enzymes)
}



#'\code{map_pathways_to_metabolites}
#'
#'Function to map pathways to enzymes and subsequently to metabolites.
#'
#'@param  metab_df data frame with metabolites and their target/source enzymes 
#'@return data frame containing metabolites and their pathways
#'@import sjmisc
#'@import dplyr

map_pathways_to_metabolites <- function(metab_df){

  #Get pure gene names by removing "_reverse" and "Enzyme[0-9]+__"
  metab_df$enzymes <- gsub("Enzyme[0-9]+__","",metab_df$enzymes)
  metab_df$enzymes <- gsub("_reverse","",metab_df$enzymes)
  
  ##Split complexes and map metabolites to every enzyme of complex
  for (i in 1:nrow(metab_df))     #iterate through all rows in "enzymes" column
  {
    if(str_contains(metab_df$enzymes[i],"_"))
    {
      split_str <- as.vector(unlist(strsplit(as.character(metab_df$enzymes[i]),
                                             split = "_")))
      
      #create new row for every element in split_str (from the 2nd element on) 
      for (j in 2:length(split_str))
      {
        new_row <- c(split_str[j], metab_df$metabolites[i])
        metab_df <- rbind(metab_df, new_row)
      }
      metab_df$enzymes[i] <- split_str[1]
    }
  }
  
  ##Keep only unique rows (if a complex consisted of two or more similar enzymes)
  metab_df <- distinct(metab_df)
  
  ##Reorder rows by two conditions: 1. column "enzymes", 2. column "metabolites"
  metab_df <- metab_df[order(metab_df$enzymes, metab_df$metabolites), ]
  
  ##Add new column "pathways" by merging data frame with pathways df (by enzymes)
  metab_pathway_df = merge(metab_df[, c("enzymes", "metabolites")], 
                           pathways[, c("gene", "term")], 
                           by.x = "enzymes",         
                           by.y = "gene",
                           all.x = TRUE)      
  
  colnames(metab_pathway_df)[colnames(metab_pathway_df) == "term"] <- "pathway"
  
  metab_pathway_df$enzymes <- NULL                    #remove column "enzymes"
  
  ##Reorder rows by column "metabolites"
  metab_pathway_df <- metab_pathway_df[order(metab_pathway_df$metabolites), ]
  metab_pathway_df <- na.omit(metab_pathway_df)       #remove rows containing NA
  metab_pathway_df <- distinct(metab_pathway_df)      #keep only unique rows
  
  return(metab_pathway_df)
}



#'\code{highlight_pathways}
#'
#'Function to map the most significant pathway to every node in the Cosmos solution network.
#'This function can be used after an ORA has been performed.
#'
#'@param  pathwaysMSigDB Canonical pathways with annotated genes as Gene Symbols from MSigDB
#'                       (example: c2.cp.v7.4.symbols.gmt)
#'@param  sigPathwaysDf  Dataframe containing all significant pathways as a result of the ORA
#'@param  cosmosResults  COSMOS output (list containing two dataframes: network and nodes attributes)
#'@param  pathwaysInt    Vector containing interesting pathways which should be
#'                       highlighted in a visualized network (using Cytoscape)
#'@param  filterMethod   "p-value" or "Adjusted p-value" can be used for filtering the hits                     
#'@return Dataframe containing the nodes attributes of the Cosmos output and
#'        additionally pathways for every node if applicable and colors
#'@export

highlight_pathways <- function(pathwaysMSigDB, sigPathwaysDf, cosmosResults,
                              pathwaysInt, filterMethod){
  
  ##Function for adding new column for p-value by merging pathways data frame 
  # with result of ORA (sigPathwaysDf) by pathway column
  map_most_sig_pathway_to_node <- function(pathways, sigPathways){
    
    pathways_pvalue = merge(pathways[, c(1, 2)], 
                              sigPathways[, c("pathway", filterMethod)], 
                              by.x = "pathway", 
                              by.y = "pathway",
                              all.x = TRUE)    
    
    pathways_pvalue <- pathways_pvalue[, c(2, 1, 3)]
    
    ##Sort by two conditions: 1. metabolite 2. Increasing p-value
    pathways_pvalue <- pathways_pvalue[order(pathways_pvalue[,1],
                                             pathways_pvalue[,3]), ]
    ##Keep only most significant pathway per node (the one with lowest p-value)
    pathways_unique <- pathways_pvalue[!duplicated(pathways_pvalue[, 1]),]
    
    return(pathways_unique)
  }
  
  ### Enzymes
  
  ##Select most significant pathway per enzyme
  colnames(pathwaysMSigDB)[colnames(pathwaysMSigDB) == "term"] <- "pathway"
  pathways_nodes_unique <- map_most_sig_pathway_to_node(pathwaysMSigDB, sigPathwaysDf)

  ##Create Nodes data frame and remove all metabolite rows
  nodes_attributes_enzymes <- cosmosResults$nodesAttributes
  nodes_attributes_enzymes <- nodes_attributes_enzymes[nodes_attributes_enzymes$type != "metabolite", ]
  
  ##Make sure that also enzymes with modified names can be mapped to pathways
  Nodes <- nodes_attributes_enzymes$Nodes
  nodes_attributes_enzymes$Nodes <- gsub("Enzyme[0-9]+__","",nodes_attributes_enzymes$Nodes)
  nodes_attributes_enzymes$Nodes <- gsub("_reverse","",nodes_attributes_enzymes$Nodes)
  nodes_attributes_enzymes$ID <- 1:nrow(nodes_attributes_enzymes)
  
  ##Map most significant pathway per gene to COSMOS output
  full_att_pathway_e = merge(nodes_attributes_enzymes[, c("Nodes", "measured",
                                                          "type", "Activity", "ID")], 
                           pathways_nodes_unique[, c("gene", "pathway")], 
                           by.x = "Nodes",  #merging by nodes column
                           by.y = "gene",
                           all.x = TRUE)    
  
  ##Change back to previous names in "Nodes" column, e.g. enzyme_reverse
  full_att_pathway_e <- full_att_pathway_e[order(full_att_pathway_e$ID), ]
  full_att_pathway_e[, c("Nodes", "ID")] <- NULL
  full_att_pathway_e <- cbind(Nodes, full_att_pathway_e)
  
  
  ### Metabolites
  
  ##Create data frame with metabolites and their source/target enzymes
  metabolites_enzymes_df <- create_metabolite_enzyme_df()
  
  #Map pathways to metabolites 
  metabolites_pathways <- map_pathways_to_metabolites(metabolites_enzymes_df)
  
  ##Select most significant pathway per metabolite
  metabs_pathways_unique <- map_most_sig_pathway_to_node(metabolites_pathways,
                                                         sigPathwaysDf)
  
  ##Create Nodes data frame and keep only metabolite rows
  nodes_attributes_metabs <- cosmosResults$nodesAttributes
  nodes_attributes_metabs <- nodes_attributes_metabs[nodes_attributes_metabs$type == "metabolite", ]

  ##Map most significant pathway per metabolite to COSMOS output
  full_att_pathway_m = merge(nodes_attributes_metabs[, c("Nodes", "measured",
                                                         "type", "Activity")], 
                             metabs_pathways_unique[, c("metabolites", "pathway")], 
                             by.x = "Nodes",  #merging by nodes column
                             by.y = "metabolites",
                             all.x = TRUE)  
  
  
  ##Combine enzymes and metabolites pathways df
  full_att_pathway <- rbind(full_att_pathway_e, full_att_pathway_m)

  ##Map colors to interesting pathways
  colors <- as.vector(topo.colors(20))  #20 different colors as hex values
  length(colors) <- length(pathwaysInt) #adopt length of color vector to pathways
  full_att_pathway["color"] <- NA       #add new column "color" to data frame
  
  ##Map one unique color to every interesting pathway
  for(i in 1:length(pathwaysInt)) {
    full_att_pathway$color[full_att_pathway$pathway == pathwaysInt[i]] <- 
      (full_att_pathway$color[full_att_pathway$pathway == pathwaysInt[i]] <- colors[i])
  }
  
  return(full_att_pathway)
}
