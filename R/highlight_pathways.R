################ Supporting function for highlighting pathways #################

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

# library(dplyr)
# library(magrittr)


#'\code{highlightPathways}
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
#'@return Dataframe containing the nodes attributes of the Cosmos output and
#'        additionally pathways for every node if applicable and colors

highlightPathways <- function(pathwaysMSigDB, sigPathwaysDf, cosmosResults,
                              pathwaysInt){
  
  # Todo Caroline: replace "p-value" two times with "Adjusted p-value"
  
  colnames(pathwaysMSigDB) <- c("Node", "pathway")
  # add new column "Adjusted p-value" by merging pathwaysMSigDB df with 
  # result of ORA (sigPathwaysDf) by thhe column pathway
  pathways_nodes_df = merge(pathwaysMSigDB[, c("Node", "pathway")], 
                            sigPathwaysDf[, c("pathway", "p-value")], 
                            by.x = "pathway", #merging to be done by pathway col
                            by.y = "pathway",
                            all.x = TRUE)    #keep all rows from pathwaysMSigDB df
  #even if there is no matching data in 2nd df
  # change order of columns
  pathways_nodes_df <- pathways_nodes_df[c("Node", "pathway", "p-value")]
  
  # Filter by two conditions: 1. Nodes, 2. Increasing p-value
  pathways_nodes_df <- pathways_nodes_df[order(pathways_nodes_df[,1],
                                               pathways_nodes_df[,3]), ]
  # keep only most significant pathway per node (the one with lowest p-value)
  pathways_nodes_unique <- pathways_nodes_df[!duplicated(pathways_nodes_df$Nodes),]
  
  
  ## map most significant pathway per gene to COSMOS output
  
  # add new column "pathway" to COSMOS attribute output
  full_att_pathway = merge(cosmosResults$nodesAttributes[, c("Nodes",
                                                             "measured", "type",
                                                             "Activity")], 
                           pathways_nodes_unique[, c("Node", "pathway")], 
                           by.x = "Nodes",   #merging to be done by Nodes column
                           by.y = "Node",
                           all.x = TRUE)    
  
  ### Metabolites (to be adjusted)
  ## map pathway to metabolites in COSMOS attribute output
  #metab_pathway <- readRDS("C:/Users/Admin/Documents/MASTER/4Semester/Internship_1_Saez-Rodriguez/Project_2_Highlight_interesting_pathways/metabolite_pathways.RData")
  #colnames(metab_pathway)[colnames(metab_pathway) == "pathway"] <- "pathway_m"
  #full_att_pathway = merge(cosmosResults$full_att[, c("Nodes", "measured", "type",
  #                                                    "Activity", "pathway")], 
  #                         metab_pathway[, c("metabolites", "pathway_m")], 
  #                         by.x = "Nodes",   #merging to be done by Nodes column
  #                         by.y = "metabolites",
  #                         all.x = TRUE)
  #full_att_pathway <- full_att_pathway(coalesce("pathway", "pathway_m"))
  
  
  ## map colors to interesting pathways
  
  colors <- as.vector(topo.colors(20))  #20 different colors as hex values
  length(colors) <- length(pathwaysInt) #adopt length of col vector to length of pathway vector
  full_att_pathway["border_color"] <- NA       #add new column 'color' to dataframe
  
  # mapping colors to selected pathways
  for(i in 1:length(pathwaysInt)) {
    full_att_pathway$border_color[full_att_pathway$pathway == pathwaysInt[i]] <- 
      (full_att_pathway$border_color[full_att_pathway$pathway == pathwaysInt[i]] <- colors[i])
  }
  
  return(full_att_pathway)
}