# Central_node = a single node or a vector of nodes in character format
# 
# N = number of max steps to reach measured nodes

display_node_neighbourhood <- function(central_node,sif, att, n = 100)
{
    full_sif <- sif
    full_att <- att
    
    ig_net <- graph_from_data_frame(full_sif) 
    
    ig_net <- make_ego_graph(ig_net, nodes = central_node, order = n, mode = "all")
    
    to_keep <- unlist(sapply(ig_net,function(x){V(x)$name}))
    
    full_sif <- full_sif[full_sif$Node1 %in% to_keep & full_sif$Node2 %in% to_keep,]
    full_att <- full_att[full_att$Nodes %in% to_keep,]
    
    ig_net <- graph_from_data_frame(full_sif) 
    
    center_node <- sapply(central_node, function(x, ig_net, full_att) {
        names(unlist(get.shortest.paths(ig_net, from = x, to = full_att[full_att$measured == 1,1])$vpath))
    }, ig_net= ig_net, full_att = full_att, simplify = F, USE.NAMES = T)
    
    center_node <- unlist(center_node)
    
    center_node_out <- full_sif[full_sif$Node1 %in% center_node & full_sif$Node2 %in% center_node,]
    
    center_node <- sapply(central_node, function(x, ig_net, full_att) {
        names(unlist(get.shortest.paths(ig_net, from = x, to = full_att[full_att$measured == 1,1], mode = "in")$vpath))
    }, ig_net= ig_net, full_att = full_att, simplify = F, USE.NAMES = T)
    
    center_node <- unlist(center_node)
    
    center_node_in <- full_sif[full_sif$Node1 %in% center_node & full_sif$Node2 %in% center_node,]
    
    # center_node <- get.shortest.paths(ig_net, from = central_node, to = full_att[full_att$measured == 1,1])
    # center_node_out <- full_sif[full_sif$Node1 %in% names(unlist(center_node$vpath)) & full_sif$Node2 %in% names(unlist(center_node$vpath)),]
    # 
    # # write_csv(center_node_net,"center_node_sif_newDoro.csv")
    # 
    # center_node <- get.shortest.paths(ig_net, from = central_node[2], to = full_att[full_att$measured == 1,1], mode = "in",  output = "vpath")
    # center_node_in <- full_sif[full_sif$Node1 %in% names(unlist(center_node$vpath)) & full_sif$Node2 %in% names(unlist(center_node$vpath)),]
    # 
    center_node_net <- as.data.frame(rbind(center_node_in,center_node_out))
    
    nodes <- full_att[full_att$Nodes %in% center_node_net$Node1 | full_att$Nodes %in% center_node_net$Node2,]
    edges <- center_node_net
    
    names(edges) <- c("from","to","sign","weigth")
    edges$color <- ifelse(edges$sign == 1, "green","red")
    edges$arrows <- "to"
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
    nodes$shadow <- ifelse(nodes$measured == 1, T, F)
    
    return(visNetwork(nodes, edges, width = 1600, height = 1600) %>% 
               visOptions(highlightNearest = TRUE, 
                          nodesIdSelection = list(enabled = TRUE,
                                                  style = 'width: 200px; height: 26px;
                                             background: #f8f8f8;
                                             color: darkblue;
                                             border:none;
                                             outline:none;')))
    
}