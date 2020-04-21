library(CARNIVAL) # load CARNIVAL library
library(biomaRt)
library(gsubfn)
library(readr)
library(rpubchem)
library(stringr)

setwd("~/Dropbox/COSMOS/")

CARNIVAL_Result <- runCARNIVAL(CplexPath="~/Documents/cplex",
                               Result_dir="results/multi_omic/carnival/",
                               CARNIVAL_example=NULL,
                               UP2GS=F,
                               netFile = "support/meta_network_with_X_nobadchar.tsv", 
                               measFile = "results/metabolomic/metab_input_carnival.tsv", 
                               inputFile = "results/phospho/signaling_input_carnival.tsv",
                               weightFile = NULL, 
                               timelimit = 7200,
                               mipGAP = 0.12
                               
) #11.24

CARNIVAL_Result2 <- runCARNIVAL(CplexPath="~/Documents/cplex",
                                Result_dir="results/multi_omic/carnival/",
                                CARNIVAL_example=NULL,
                                UP2GS=F,
                                netFile = "support/meta_network_with_X_nobadchar.tsv", 
                                measFile = "results/phospho/signaling_input_carnival.tsv", 
                                inputFile = "results/metabolomic/metab_input_carnival.tsv",
                                weightFile = NULL, 
                                timelimit = 7200
)# 2.44

signaling_input_carnival <- as.data.frame(
  read_delim("results/phospho/signaling_input_carnival.tsv",
             "\t", escape_double = FALSE, trim_ws = TRUE))

metab_input_carnival <- as.data.frame(
  read_delim("results/metabolomic/metab_input_carnival.tsv",
             "\t", escape_double = FALSE, trim_ws = TRUE))

sif <- as.data.frame(CARNIVAL_Result[[1]]$weightedSIF)
sif2 <- as.data.frame(CARNIVAL_Result2[[1]]$weightedSIF)

sif$Node1 <- gsub("^X","",sif$Node1)
sif$Node2 <- gsub("^X","",sif$Node2)

sif2$Node1 <- gsub("^X","",sif2$Node1)
sif2$Node2 <- gsub("^X","",sif2$Node2)

att <- as.data.frame(CARNIVAL_Result[[1]]$attributesAll)#[,c(1,2)]
att$measured <- ifelse(att$Nodes %in% names(signaling_input_carnival) | att$Nodes %in% names(metab_input_carnival), 1, 0)
att$Nodes <- gsub("^X","",att$Nodes)

att2 <- as.data.frame(CARNIVAL_Result2[[1]]$attributesAll)#[,c(1,2)]
att2$measured <- ifelse(att2$Nodes %in% names(signaling_input_carnival) | att2$Nodes %in% names(metab_input_carnival), 1, 0)
att2$Nodes <- gsub("^X","",att2$Nodes)

att$type <- ifelse(grepl("Metab",att$Nodes), "metabolite","protein")
att2$type <- ifelse(grepl("Metab",att2$Nodes), "metabolite","protein")


########################
prots <- unique(c(att$Nodes, att2$Nodes))
prots <- prots[!(grepl("Metab",prots))]
prots <- gsub("Gene[0-9]+__","",prots)
prots <- gsub("_reverse","",prots)
prots <- gsub("EXCHANGE.*","",prots)
prots <- unique(prots)

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

G_list <- getBM(filters = "entrezgene_id",
                attributes = c('hgnc_symbol','entrezgene_id', "description"),
                values = prots, mart = ensembl)

gene_mapping <- G_list[,1]
names(gene_mapping) <- G_list[,2]

# sif$Node1 <- gsub("Gene[0-9]+__","",sif$Node1)
# sif$Node2 <- gsub("Gene[0-9]+__","",sif$Node2)
# sif2$Node1 <- gsub("Gene[0-9]+__","",sif2$Node1)
# sif2$Node2 <- gsub("Gene[0-9]+__","",sif2$Node2)
# att$Nodes <- gsub("Gene[0-9]+__","",att$Nodes)
# att2$Nodes <- gsub("Gene[0-9]+__","",att2$Nodes)

sif$Node1 <- gsub("_reverse","",sif$Node1)
sif$Node2 <- gsub("_reverse","",sif$Node2)
sif2$Node1 <- gsub("_reverse","",sif2$Node1)
sif2$Node2 <- gsub("_reverse","",sif2$Node2)
att$Nodes <- gsub("_reverse","",att$Nodes)
att2$Nodes <- gsub("_reverse","",att2$Nodes)

sif$Node1 <- gsub("EXCHANGE.*","",sif$Node1)
sif$Node2 <- gsub("EXCHANGE.*","",sif$Node2)
sif2$Node1 <- gsub("EXCHANGE.*","",sif2$Node1)
sif2$Node2 <- gsub("EXCHANGE.*","",sif2$Node2)
att$Nodes <- gsub("EXCHANGE.*","",att$Nodes)
att2$Nodes <- gsub("EXCHANGE.*","",att2$Nodes)

metabs <- unique(c(att$Nodes, att2$Nodes))
metabs <- metabs[(grepl("Metab",metabs))]

metab_to_pubchem <- as.data.frame(read_csv("support/metab_to_pubchem.csv"))
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
}

for(i in 1:length(sif2$Node1))
{
  if(gsub("Gene[0-9]+__","",sif2[i,1]) %in% names(gene_mapping))
  {
    if(grepl("Gene",sif2[i,1]))
    {
      prefix <- gsub("__.*","",sif2[i,1])
      sif2[i,1] <- paste(prefix,gene_mapping[gsub("Gene[0-9]+__","",sif2[i,1])], sep = "__")
    }
    else
    {
      sif2[i,1] <- gene_mapping[gsub("Gene[0-9]+__","",sif2[i,1])]
    }
  }
  if(gsub("Gene[0-9]+__","",sif2[i,3]) %in% names(gene_mapping))
  {
    if(grepl("Gene",sif2[i,3]))
    {
      prefix <- gsub("__.*","",sif2[i,3])
      sif2[i,3] <- paste(prefix,gene_mapping[gsub("Gene[0-9]+__","",sif2[i,3])], sep = "__")
    }
    else
    {
      sif2[i,3] <- gene_mapping[gsub("Gene[0-9]+__","",sif2[i,3])]
    }
  }
  if(gsub("Metab__","",gsub("___[a-z]____","",sif2[i,1])) %in% names(metab_to_pubchem_vec))
  {
    suffix <- str_match(sif2[i,1],"___.____")
    metab <- metab_to_pubchem_vec[gsub("Metab__","",gsub("___[a-z]____","",sif2[i,1]))]
    sif2[i,1] <- paste(metab,suffix,sep = "")
  }
  if(gsub("Metab__","",gsub("___[a-z]____","",sif2[i,3])) %in% names(metab_to_pubchem_vec))
  {
    suffix <- str_match(sif2[i,3],"___.____")
    metab <- metab_to_pubchem_vec[gsub("Metab__","",gsub("___[a-z]____","",sif2[i,3]))]
    sif2[i,3] <- paste(metab,suffix,sep = "")
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
}

for(i in 1:length(att2$Nodes))
{
  if(gsub("Gene[0-9]+__","",att2[i,1]) %in% names(gene_mapping))
  {
    if(grepl("Gene",att2[i,1]))
    {
      prefix <- gsub("__.*","",att2[i,1])
      att2[i,1] <- paste(prefix,gene_mapping[gsub("Gene[0-9]+__","",att2[i,1])], sep = "__")
    }
    else
    {
      att2[i,1] <- gene_mapping[gsub("Gene[0-9]+__","",att2[i,1])]
    }
  }
  if(gsub("Metab__","",gsub("___[a-z]____","",att2[i,1])) %in% names(metab_to_pubchem_vec))
  {
    suffix <- str_match(att2[i,1],"___.____")
    metab <- metab_to_pubchem_vec[gsub("Metab__","",gsub("___[a-z]____","",att2[i,1]))]
    att2[i,1] <- paste(metab,suffix,sep = "")
  }
}

########################
url <- paste0(
  'http://omnipathdb.org/ptms?',
  'fields=sources,references&genesymbols=1'
)

download_omnipath <- function(){
  
  read.table(url, sep = '\t', header = TRUE)
  
}

omnipath_ptm <- download_omnipath()
omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation","phosphorylation"),]
KSN <- omnipath_ptm[,c(4,3)]
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_type, sep ="_")
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_offset, sep = "")
KSN$sign <- ifelse(omnipath_ptm$modification == "phosphorylation", 1, -1)

att$type <- ifelse(att$Nodes %in% KSN$enzyme_genesymbol, "Kinase",att$type)

url <- paste0(
  'http://omnipathdb.org/interactions?',
  'datasets=tfregulons&tfregulons_levels=A,B,C&genesymbols=1&fields=sources,tfregulons_level'
)

dorothea <- download_omnipath()
dorothea <- dorothea[,c(4,3,6,7)]
dorothea$sign <- dorothea$is_stimulation - dorothea$is_inhibition
dorothea$sign <- ifelse(dorothea$sign == 0, 1, dorothea$sign)
dorothea <- dorothea[,c(1,2,5)]

att$type <- ifelse(att$Nodes %in% dorothea$source_genesymbol, "TF",att$type)

i <- 1
for(node in att$Nodes) #check if a node name is basename appears more than once in nodes. If so, it is a metabolic enzyme
{
  if(sum(gsub("Gene.*_","",node) == gsub("Gene.*_","",att$Nodes)) > 1)
  {
    att[i,4] <- "metab_enzyme"
  }
  i <- i+1
}

att$type <- ifelse(grepl("Gene.*_",att$Nodes), "metab_enzyme",att$type)

att2$type <- ifelse(att2$Nodes %in% KSN$enzyme_genesymbol, "Kinase",att2$type)

att2$type <- ifelse(att2$Nodes %in% dorothea$source_genesymbol, "TF",att2$type)

i <- 1
for(node in att2$Nodes) #check if a node name is basename appears more than once in nodes. If so, it is a metabolic enzyme
{
  if(sum(gsub("Gene.*_","",node) == gsub("Gene.*_","",att2$Nodes)) > 1)
  {
    att2[i,4] <- "metab_enzyme"
  }
  i <- i+1
}

att2$type <- ifelse(grepl("Gene.*_",att2$Nodes), "metab_enzyme",att2$type)

att$Activity <- sign(as.numeric(as.character(att$Activity)))
att2$Activity <- sign(as.numeric(as.character(att2$Activity)))

sif$Node1 <- gsub("Gene[0-9]+__","",sif$Node1)
sif$Node2 <- gsub("Gene[0-9]+__","",sif$Node2)
att$Nodes <- gsub("Gene[0-9]+__","",att$Nodes)

sif2$Node1 <- gsub("Gene[0-9]+__","",sif2$Node1)
sif2$Node2 <- gsub("Gene[0-9]+__","",sif2$Node2)
att2$Nodes <- gsub("Gene[0-9]+__","",att2$Nodes)

sif <- sif[sif$Node1 != sif$Node2,]
sif2 <- sif2[sif2$Node1 != sif2$Node2,]

sif$Node1 <- gsub("[_]{3,5}","_",sif$Node1)
sif$Node2 <- gsub("[_]{3,5}","_",sif$Node2)
att$Nodes <- gsub("[_]{3,5}","_",att$Nodes)

sif2$Node1 <- gsub("[_]{3,5}","_",sif2$Node1)
sif2$Node2 <- gsub("[_]{3,5}","_",sif2$Node2)
att2$Nodes <- gsub("[_]{3,5}","_",att2$Nodes)

########################

write_csv(sif, "results/multi_omic/carnival/subnet_sif_full_wTF.csv")

write_csv(att, "results/multi_omic/carnival/subnet_sif_full_att_wTF.csv")

write_csv(sif2, "results/multi_omic/carnival/subnet_sif2_full_wTF.csv")

write_csv(att2, "results/multi_omic/carnival/subnet_sif_full_att2_wTF.csv")

sif3 <- as.data.frame(rbind(sif[as.numeric(as.character(sif$Weight)) >= 80,], sif2[as.numeric(as.character(sif2$Weight)) >= 80,]))
write_csv(sif3, "results/multi_omic/carnival/subnet_sif3_wTF.csv")

save.image("results/multi_omic/carnival/results_CARNIVAL.Rdata")
