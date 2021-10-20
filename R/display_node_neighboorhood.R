#' display_node_neighboorhood
#' 
#' display input and measurements within n steps of a given set of nodes
#' @param central_node character or character vector; node ID(s) around which a 
#' network will be branched out untill meansurments and input are reached
#' @param sif df; COSMOS network solution in sif format like the first list
#' element returned by the format_cosmos_res function
#' @param att df; attributes of the nodes of the COMSOS network solution like 
#' the second list element returned by the format_cosmos_res function
#' @param n numeric; maximum number of steps in the network to look for inputs 
#' and measurments
#' @return a visnetwork object
#' @importFrom magrittr %>%
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
#' data(metabolite_to_pubchem)
#' data(omnipath_ptm)
#' test_result_for <- format_COSMOS_res(test_result_for,
#' metab_mapping = metabolite_to_pubchem,
#' measured_nodes = unique(c(names(toy_metabolic_input),
#'                           names(toy_signaling_input))),
#' omnipath_ptm = omnipath_ptm)
#' network_plot <- display_node_neighboorhood(central_node = 'NFKB1',
#' sif = test_result_for[[1]],
#' att = test_result_for[[2]],
#' n = 7)
#' network_plot
#' @export
display_node_neighboorhood <- function(central_node,sif, att, n = 100)
{
  # require(igraph)
  full_sif <- sif
  if(sum(names(full_sif) %in% c("Node1", "Node2", "Sign", "Weight")) != 4)
  {
    print('sif should be a dataframe with 4 column named "Node1", "Node2", "Sign" and "Weight"')
    return('Error:bad sif')
  }
  full_sif <- full_sif[,c("Node1", "Node2", "Sign", "Weight")]
  full_att <- att
  
  ig_net <- igraph::graph_from_data_frame(full_sif) 
  
  ig_net <- igraph::make_ego_graph(ig_net, nodes = central_node, order = n, mode = "all")
  
  to_keep <- unlist(sapply(ig_net,function(x){igraph::V(x)$name}))
  
  full_sif <- full_sif[full_sif$Node1 %in% to_keep & full_sif$Node2 %in% to_keep,]
  full_att <- full_att[full_att$Nodes %in% to_keep,]
  
  ig_net <- igraph::graph_from_data_frame(full_sif) 
  
  center_node <- sapply(central_node, function(x, ig_net, full_att) {
    names(unlist(igraph::get.shortest.paths(ig_net, from = x, to = full_att[full_att$measured == 1,1])$vpath))
  }, ig_net= ig_net, full_att = full_att, simplify = FALSE, USE.NAMES = TRUE)
  
  center_node <- unlist(center_node)
  
  center_node_out <- full_sif[full_sif$Node1 %in% center_node & full_sif$Node2 %in% center_node,]
  
  center_node <- sapply(central_node, function(x, ig_net, full_att) {
    names(unlist(igraph::get.shortest.paths(ig_net, from = x, to = full_att[full_att$measured == 1,1], mode = "in")$vpath))
  }, ig_net= ig_net, full_att = full_att, simplify = FALSE, USE.NAMES = TRUE)
  
  center_node <- unlist(center_node)
  
  center_node_in <- full_sif[full_sif$Node1 %in% center_node & full_sif$Node2 %in% center_node,]

  center_node_net <- as.data.frame(rbind(center_node_in,center_node_out))
  
  nodes <- full_att[full_att$Nodes %in% center_node_net$Node1 | full_att$Nodes %in% center_node_net$Node2,]
  edges <- center_node_net
  
  names(edges) <- c("from","to","sign","weigth")
  edges$color <- ifelse(edges$sign == 1, "grey","grey")
  
  # edges$arrows <- "to"
  
  edges$arrows.to.type <- ifelse(edges$sign == 1, "arrow","circle")
  edges$enabled <- TRUE
  edges$scaleFactor <- 1
  
  edges <- unique(edges)
  
  names(nodes)[1] <- "id"
  nodes$label <- nodes$id
  nodes$color <- ifelse(nodes$Activity > 0, "green","red")
  nodes <- nodes[!duplicated(nodes$id),]
  nodes$shape <- "dot" 
  nodes[nodes$type == "metab_enzyme","shape"] <- "square"
  nodes[nodes$type == "protein","shape"] <- "square"
  nodes[nodes$type == "Kinase","shape"] <- "triangle"
  nodes[nodes$type == "TF","shape"] <- "diamond"
  nodes <- nodes[order(nodes$id),]
  nodes$shadow <- (nodes$measured == 1)
  
  return(visNetwork::visNetwork(nodes, edges, width = 1600, height = 1600) %>% 
             visNetwork::visOptions(highlightNearest = TRUE, 
                      nodesIdSelection = list(enabled = TRUE,
                                              style = 'width: 200px; height: 26px;
                                              background: #f8f8f8;
                                              color: darkblue;
                                              border:none;
                                              outline:none;')))
  
}