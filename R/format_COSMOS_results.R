cleanup_protein_node <- function(node)
{
  prefixe <- stringr::str_extract(node, "Gene[0-9]+__")
  prefixe <- gsub("Gene","Enzyme",prefixe)
  
  suffixe_direction <- stringr::str_extract(node, "_reverse")
  suffixe_exchange <- stringr::str_extract(node, "EXCHANGE[0-9]+")
  
  node <- gsub("Gene[0-9]+__","",node)
  node <- gsub("_reverse","",node)
  node <- gsub("EXCHANGE[0-9]+","",node)
  
  node <- stringr::str_split(node,'_')[[1]]
  
  return(list("prefixe" = prefixe, "node" = node,"suffixe_direction" = suffixe_direction, "suffixe_exchange" = suffixe_exchange))
}

cleanup_metab_node <- function(node)
{
  prefixe <- 'Metab__'
  
  suffixe <- stringr::str_extract(node,'___[a-z]____')
  
  node <- gsub('Metab__','',node)
  node <- gsub('___[a-z]____','',node)
  
  return(list("prefixe" = prefixe, "node" = node,"suffixe_compartment" = suffixe))
}

translate_column <- function(sif_column,
                             metab_mapping,
                             gene_mapping)
{
  # sif_column_length <- length(sif_column)
  sif_column <- vapply(sif_column, function(x,
                                            metab_mapping,
                                            gene_mapping)
    {

    node <- x
    
    type_node <- ifelse(grepl("Metab",node),'metabolite','protein')
    
    if(type_node == 'protein')
    {
      node <- cleanup_protein_node(node)
      
      node$node <- unlist(vapply(node$node, function(x,
                                              gene_mapping)
        {
        return(ifelse(x %in% names(gene_mapping),
                      gene_mapping[x],
                      x))
      }, character(1),
      gene_mapping = gene_mapping, USE.NAMES = FALSE))
      
      node <- node[!is.na(node)]
      
      node$node <- paste(node$node, collapse = '_')
      
      node <- paste(node, collapse = '')
  
    } else
    {
      node <- cleanup_metab_node(node)
      
      node$node <- unlist(vapply(node$node, function(x,
                                                     metab_mapping)
      {
        return(ifelse(x %in% names(metab_mapping),
                      metab_mapping[x],
                      x))
      }, character(1),
      metab_mapping = metab_mapping, USE.NAMES = FALSE))
      
      node <- node[!is.na(node)]
      node <- paste(node, collapse = '')
    }
    
    
      
  }, character(1), 
  metab_mapping = metab_mapping,
  gene_mapping = gene_mapping, USE.NAMES = FALSE)
  
  return(sif_column)
}

translate_sif <- function(sif,
                          metab_mapping,
                          gene_mapping)
{
  sif$Node1 <- translate_column(sif$Node1,
                                metab_mapping,
                                gene_mapping)
  
  sif$Node2 <- translate_column(sif$Node2,
                                metab_mapping,
                                gene_mapping)
  
  return(sif)
}


#' format_COSMOS_res
#'
#' formats the network with readable names
#'
#' @param cosmos_res  results of CARNIVAL run
#' @param metab_mapping a named vector with pubchem cid as names and desired metabolite names as values.
#' @param gene_mapping by default, use the 'org.Hs.eg.db' to map gene names. Can also be a named vector with entrez gene id as names and desired gene names as values.
#' @param measured_nodes vector of nodes that are measured or inputs
#' @param omnipath_ptm ptms database from OmnipathR
#' @return list with network and attribute tables.
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
#' @export
format_COSMOS_res <- function(cosmos_res,
                              metab_mapping,
                              gene_mapping = 'org.Hs.eg.db',
                              measured_nodes,
                              omnipath_ptm)
{
  # require(dorothea)
  
  
  sif <- as.data.frame(cosmos_res$weightedSIF)
  sif$Node1 <- gsub("^X","",sif$Node1)
  sif$Node2 <- gsub("^X","",sif$Node2)
  att <- as.data.frame(cosmos_res$nodesAttributes)#[,c(1,2)]
  names(att)[1] <- "Nodes"
  att$measured <- ifelse(att$Nodes %in% measured_nodes, 1, 0)
  att$Nodes <- gsub("^X","",att$Nodes)
  att$type <- ifelse(grepl("Metab",att$Nodes), "metabolite","protein")
  att <- att[abs(as.numeric(att$AvgAct)) > 0,]
  
  ########################
  prots <- unique(att$Nodes)
  prots <- prots[!(grepl("Metab",prots))]
  prots <- gsub("Gene[0-9]+__","",prots)
  prots <- gsub("_reverse","",prots)
  prots <- gsub("EXCHANGE.*","",prots)
  prots <- unique(prots)
  prots <- unlist(sapply(prots, function(x){unlist(strsplit(x,split = "_"))}))
  
  if(is.null(names(gene_mapping)))
  {
    entrezid <- prots
    gene_mapping <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, entrezid, 'SYMBOL', 'ENTREZID')
    gene_mapping <- unlist(gene_mapping)
    gene_mapping <- gene_mapping[!is.na(gene_mapping)]
  } else
  {
    if(sum(names(gene_mapping) %in% prots) == 0)
    {
      print('Error : none of gene mapping identifiers are not found in the network solution. Please make sure you inputed a properlly named vector. Else use the default value for this argument : "org.Hs.eg.db"')
      return('Bad identifer mapping')
    }
  }
  
  sif <- translate_sif(sif,
                       metab_mapping,
                       gene_mapping)
  
  att$Nodes <- translate_column(att$Nodes,
                                metab_mapping,
                                gene_mapping)
  
  ########################
  omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation","phosphorylation"),]
  KSN <- omnipath_ptm[,c(4,3)]
  KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_type, sep ="_")
  KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_offset, sep = "")
  KSN$sign <- ifelse(omnipath_ptm$modification == "phosphorylation", 1, -1)
  att$type <- ifelse(att$type == 'protein', 'metab_enzyme', att$type)
  
  att$type <- ifelse(att$Nodes %in% KSN$enzyme_genesymbol, "Kinase",att$type)
  dorothea <- as.data.frame(dorothea::dorothea_hs[dorothea::dorothea_hs$confidence %in% c("A","B","C","D","E"),c(3,1,4)])
  names(dorothea) <- c("target_genesymbol","source_genesymbol","sign")
  att$type <- ifelse(att$Nodes %in% dorothea$source_genesymbol, "TF",att$type)

  att$Activity <- sign(as.numeric(as.character(att$AvgAct)))

  sif <- sif[sif$Node1 != sif$Node2,]

  sif <- sif[sif$Node1 %in% att$Nodes & sif$Node2 %in% att$Nodes,]
  
  measured <- att[att$measured == 1,"Nodes"]
  att$measured <- ifelse(att$Nodes %in% measured, 1, 0)
  
  return(list(sif,att))
}


#' graphs from COSMOS results
#'
#' formats the network with readable names
#'
#' @param cosmos_res  results of CARNIVAL run
#' @param metab_mapping mapping table between metabolite naming conventions,
#' a two column dataframe with names: c("name","pubchem")
#' @param gene_mapping -- not sure !!!
#' @param measured_nodes vector of nodes that are measured or inputs
#' @param omnipath_ptm ptms database from OmnipathR
#' @noRd
format_COSMOS_results_deprecated <- function(cosmos_res,
                              metab_mapping,
                              gene_mapping = "ensembl",
                              measured_nodes,
                              omnipath_ptm)
{
    
    
    sif <- as.data.frame(cosmos_res$weightedSIF)
    
    sif$Node1 <- gsub("^X","",sif$Node1)
    sif$Node2 <- gsub("^X","",sif$Node2)
    
    
    att <- as.data.frame(cosmos_res$nodesAttributes)#[,c(1,2)]
    names(att)[1] <- "Nodes"
    att$measured <- ifelse(att$Nodes %in% measured_nodes, 1, 0)
    att$Nodes <- gsub("^X","",att$Nodes)
    
    
    att$type <- ifelse(grepl("Metab",att$Nodes), "metabolite","protein")
    
    att <- att[abs(as.numeric(att$AvgAct)) > 0,]
    
    ########################
    prots <- unique(att$Nodes)
    prots <- prots[!(grepl("Metab",prots))]
    prots <- gsub("Gene[0-9]+__","",prots)
    prots <- gsub("_reverse","",prots)
    prots <- gsub("EXCHANGE.*","",prots)
    prots <- unique(prots)
    prots <- unlist(sapply(prots, function(x){unlist(strsplit(x,split = "_"))}))
    
    
    if(gene_mapping == "ensembl")
    {
        ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
        
        G_list <- biomaRt::getBM(filters = "entrezgene_id",
                        attributes = c('hgnc_symbol','entrezgene_id', "description"),
                        values = prots, mart = ensembl)
        
    } else
    {
        G_list <- as.data.frame(readr::read_csv(gene_mapping))
    }
    
    gene_mapping <- G_list[,1]
    names(gene_mapping) <- G_list[,2]
    
    
    sif$Node1 <- gsub("_reverse","",sif$Node1)
    sif$Node2 <- gsub("_reverse","",sif$Node2)
    
    att$Nodes <- gsub("_reverse","",att$Nodes)
    
    
    sif$Node1 <- gsub("EXCHANGE.*","",sif$Node1)
    sif$Node2 <- gsub("EXCHANGE.*","",sif$Node2)
    
    att$Nodes <- gsub("EXCHANGE.*","",att$Nodes)
    
    
    metabs <- unique(c(att$Nodes))
    metabs <- metabs[(grepl("Metab",metabs))]
    
    metab_to_pubchem <- as.data.frame(metab_mapping)
    metab_to_pubchem_vec <- metab_to_pubchem$name
    names(metab_to_pubchem_vec) <- metab_to_pubchem$pubchem
    
    for(i in seq_along(sif$Node1))
    {
        if(gsub("Gene[0-9]+__","",sif[i,1]) %in% names(gene_mapping))
        {
            if(grepl("Gene",sif[i,1]))
            {
                prefix <- gsub("__.*","",sif[i,1])
                sif[i,1] <- paste(prefix,gene_mapping[gsub("Gene[0-9]+__","",sif[i,1])], sep = "__")
            }
            else
            {
                sif[i,1] <- gene_mapping[gsub("Gene[0-9]+__","",sif[i,1])]
            }
        }
        if(gsub("Gene[0-9]+__","",sif[i,3]) %in% names(gene_mapping))
        {
            if(grepl("Gene",sif[i,3]))
            {
                prefix <- gsub("__.*","",sif[i,3])
                sif[i,3] <- paste(prefix,gene_mapping[gsub("Gene[0-9]+__","",sif[i,3])], sep = "__")
            }
            else
            {
                sif[i,3] <- gene_mapping[gsub("Gene[0-9]+__","",sif[i,3])]
            }
        }
        if(gsub("Metab__","",gsub("___[a-z]____","",sif[i,1])) %in% names(metab_to_pubchem_vec))
        {
            suffix <- stringr::str_match(sif[i,1],"___.____")
            metab <- metab_to_pubchem_vec[gsub("Metab__","",gsub("___[a-z]____","",sif[i,1]))]
            sif[i,1] <- paste(metab,suffix,sep = "")
        }
        if(gsub("Metab__","",gsub("___[a-z]____","",sif[i,3])) %in% names(metab_to_pubchem_vec))
        {
            suffix <- stringr::str_match(sif[i,3],"___.____")
            metab <- metab_to_pubchem_vec[gsub("Metab__","",gsub("___[a-z]____","",sif[i,3]))]
            sif[i,3] <- paste(metab,suffix,sep = "")
        }
        if(length(intersect(unlist(strsplit(gsub("Gene[0-9]+__","",sif[i,1]),split = "_")),names(gene_mapping))) > 0)
        {
            genes <- unlist(strsplit(gsub("Gene[0-9]+__","",sif[i,1]),split = "_"))
            genes <- unlist(sapply(genes, function(x){gene_mapping[x]}))
            genes <- paste0(genes, collapse = "_")
            sif[i,1] <- genes
        }
        if(length(intersect(unlist(strsplit(gsub("Gene[0-9]+__","",sif[i,3]),split = "_")),names(gene_mapping))) > 0)
        {
            genes <- unlist(strsplit(gsub("Gene[0-9]+__","",sif[i,3]),split = "_"))
            genes <- unlist(sapply(genes, function(x){gene_mapping[x]}))
            genes <- paste0(genes, collapse = "_")
            sif[i,3] <- genes
        }
    }
    
    
    for(i in seq_along(att$Nodes))
    {
        if(gsub("Gene[0-9]+__","",att[i,1]) %in% names(gene_mapping))
        {
            if(grepl("Gene",att[i,1]))
            {
                prefix <- gsub("__.*","",att[i,1])
                att[i,1] <- paste(prefix,gene_mapping[gsub("Gene[0-9]+__","",att[i,1])], sep = "__")
            }
            else
            {
                att[i,1] <- gene_mapping[gsub("Gene[0-9]+__","",att[i,1])]
            }
        }
        if(gsub("Metab__","",gsub("___[a-z]____","",att[i,1])) %in% names(metab_to_pubchem_vec))
        {
            suffix <- stringr::str_match(att[i,1],"___.____")
            metab <- metab_to_pubchem_vec[gsub("Metab__","",gsub("___[a-z]____","",att[i,1]))]
            att[i,1] <- paste(metab,suffix,sep = "")
        }
        if(length(intersect(unlist(strsplit(gsub("Gene[0-9]+__","",att[i,1]),split = "_")),names(gene_mapping))) > 0)
        {
            genes <- unlist(strsplit(gsub("Gene[0-9]+__","",att[i,1]),split = "_"))
            genes <- unlist(sapply(genes, function(x){gene_mapping[x]}))
            genes <- paste0(genes, collapse = "_")
            att[i,1] <- genes
        }
    }
    
    
    
    ########################
    
    
    omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation","phosphorylation"),]
    KSN <- omnipath_ptm[,c(4,3)]
    KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_type, sep ="_")
    KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_offset, sep = "")
    KSN$sign <- ifelse(omnipath_ptm$modification == "phosphorylation", 1, -1)
    
    att$type <- ifelse(att$Nodes %in% KSN$enzyme_genesymbol, "Kinase",att$type)
    
    dorothea <- as.data.frame(dorothea::dorothea_hs[dorothea::dorothea_hs$confidence %in% c("A","B","C","D"),c(3,1,4)])
    names(dorothea) <- c("target_genesymbol","source_genesymbol","sign")
    
    att$type <- ifelse(att$Nodes %in% dorothea$source_genesymbol, "TF",att$type)
    
    i <- 1
    for(node in att$Nodes) #check if a node name is basename appears more than once in nodes. If so, it is a metabolic enzyme
    {
        if(sum(gsub("Gene.*_","",node) == gsub("Gene.*_","",att$Nodes)) > 1)
        {
            att[i,"type"] <- "metab_enzyme"
        }
        i <- i+1
    }
    
    att$type <- ifelse(grepl("Gene.*_",att$Nodes), "metab_enzyme",att$type)
    
    
    att$Activity <- sign(as.numeric(as.character(att$AvgAct)))
    
    sif$Node1 <- gsub("Gene[0-9]+__","",sif$Node1)
    sif$Node2 <- gsub("Gene[0-9]+__","",sif$Node2)
    att$Nodes <- gsub("Gene[0-9]+__","",att$Nodes)
    
    
    sif <- sif[sif$Node1 != sif$Node2,]
    
    sif$Node1 <- gsub("[_]{3,5}","_",sif$Node1)
    sif$Node2 <- gsub("[_]{3,5}","_",sif$Node2)
    att$Nodes <- gsub("[_]{3,5}","_",att$Nodes)
    
    sif <- sif[sif$Node1 %in% att$Nodes & sif$Node2 %in% att$Nodes,]
    
    return(list(sif,att))
}