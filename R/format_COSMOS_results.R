
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
format_COSMOS_results <- function(cosmos_res,
                              metab_mapping,
                              gene_mapping = "ensembl",
                              measured_nodes,
                              omnipath_ptm)
{
    require(dorothea)
    
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
        ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
        
        G_list <- getBM(filters = "entrezgene_id",
                        attributes = c('hgnc_symbol','entrezgene_id', "description"),
                        values = prots, mart = ensembl)
        
    } else
    {
        G_list <- as.data.frame(read_csv(gene_mapping))
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
    
    for(i in 1:length(sif$Node1))
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
            suffix <- str_match(sif[i,1],"___.____")
            metab <- metab_to_pubchem_vec[gsub("Metab__","",gsub("___[a-z]____","",sif[i,1]))]
            sif[i,1] <- paste(metab,suffix,sep = "")
        }
        if(gsub("Metab__","",gsub("___[a-z]____","",sif[i,3])) %in% names(metab_to_pubchem_vec))
        {
            suffix <- str_match(sif[i,3],"___.____")
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
    
    
    for(i in 1:length(att$Nodes))
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
            suffix <- str_match(att[i,1],"___.____")
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