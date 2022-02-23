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
#' @param metab_mapping a named vector with HMDB Ids as names and desired metabolite names as values.
#' @import stringr
#' @return list with network and attribute tables.
#' @export
format_COSMOS_res <- function(cosmos_res,
                              metab_mapping = NULL)
{
  # require(dorothea)
  if(is.null(metab_mapping))
  {
    data(HMDB_mapper_vec)
  }
  SIF <- as.data.frame(cosmos_res$weightedSIF)
  ATT <- as.data.frame(cosmos_res$nodesAttributes)
  names(ATT)[1] <- "Nodes"
  ATT$measured <- ifelse(ATT$NodeType %in% c("T","S"),1,0)
  ATT$Activity <- ATT$AvgAct
  
  for(i in c(1,3))
  {
    SIF[,i] <- sapply(SIF[,i], function(x, HMDB_mapper_vec){
      x <- gsub("Metab__","",x)
      x <- gsub("^Gene","Enzyme",x)
      suffixe <- str_extract(x,"_[a-z]$")
      x <- gsub("_[a-z]$","",x)
      if(x %in% names(HMDB_mapper_vec))
      {
        x <- HMDB_mapper_vec[x]
        x <- paste("Metab__",paste(x,suffixe,sep = ""),sep = "")
      }
      return(x)
    },HMDB_mapper_vec = HMDB_mapper_vec)
  }
  
  ATT[,1] <- sapply(ATT[,1], function(x, HMDB_mapper_vec){
    x <- gsub("Metab__","",x)
    x <- gsub("^Gene","Enzyme",x)
    suffixe <- str_extract(x,"_[a-z]$")
    x <- gsub("_[a-z]$","",x)
    if(x %in% names(HMDB_mapper_vec))
    {
      x <- HMDB_mapper_vec[x]
      x <- paste("Metab__",x,sep = "")
    }
    if(!is.na(suffixe))
    {
      x <- paste(x,suffixe,sep = "")
    }
    return(x)
  },HMDB_mapper_vec = HMDB_mapper_vec)
  
  return(list(SIF,ATT))
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