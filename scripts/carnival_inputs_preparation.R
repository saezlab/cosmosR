library(readr)
library(stringr)
library('org.Hs.eg.db')

setwd("~/Dropbox/COSMOS/")

#### Preparation des input metabolomic

ttop_tumour_vs_healthy <- as.data.frame(
  read_csv("results/metabolomic/ttop_tumour_vs_healthy.csv"))
metab_to_kegg <- as.data.frame(
  read_csv("support/metab_to_kegg.txt"))
meta_network_with_X <- as.data.frame(
  read_delim("support/meta_network_carnival_ready_exch_solved.tsv", 
             "\t", escape_double = FALSE, trim_ws = TRUE))
kegg_to_pubchem <- as.data.frame(
  read_csv("support/kegg_to_pubchem.txt", 
           col_names = FALSE))

kegg_to_pubchem$X2 <- paste("XMetab__",kegg_to_pubchem$X2, sep = "")

compartment_codes <- unique(c(meta_network_with_X$source,meta_network_with_X$target))
compartment_codes <- compartment_codes[grepl("Metab",compartment_codes)]
compartment_codes <- unique(str_match(compartment_codes,"___.____"))

compartment_mapping <- list()
for(i in 1:length(compartment_codes))
{
  df <- kegg_to_pubchem
  df$X2 <- paste(df$X2, compartment_codes[i], sep = "")
  compartment_mapping[[i]] <- df
}

compartment_mapping <- as.data.frame(do.call(rbind, compartment_mapping))

compartment_mapping <- compartment_mapping[
  compartment_mapping$X2 %in% meta_network_with_X$source |
    compartment_mapping$X2 %in% meta_network_with_X$target,
  ]

kegg_to_pubchem_with_comp <- compartment_mapping
names(kegg_to_pubchem_with_comp) <- c("KEGG","pubchem")

full_mapping <- merge(metab_to_kegg, kegg_to_pubchem_with_comp, by = "KEGG")

names(ttop_tumour_vs_healthy)[1] <- "metab_name"

ttop_tumour_vs_healthy <- merge(ttop_tumour_vs_healthy, full_mapping, by = "metab_name")
ttop_tumour_vs_healthy <- ttop_tumour_vs_healthy[,c(9,2:7)]

#### Preparation des input kinase/TF

##KINASE
kinase_activities <- as.data.frame(
  read_csv("results/phospho/kinase_activities.csv"))

symbols <- kinase_activities$X1

# use mapIds method to obtain Entrez IDs
mapping_symbole_to_entrez <- mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
for(i in 1:length(kinase_activities[,1]))
{
  kinase_activities[i,1] <- mapping_symbole_to_entrez[kinase_activities[i,1]]
}

kinase_activities <- kinase_activities[complete.cases(kinase_activities),]
kinase_activities$X1 <- paste("X", kinase_activities$X1, sep = "")

##TF
TF_scores <- as.data.frame(
  read_csv("results/transcriptomic/TF_scores.csv"))

symbols <- TF_scores$ID

# use mapIds method to obtain Entrez IDs
mapping_symbole_to_entrez <- mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
for(i in 1:length(TF_scores[,1]))
{
  TF_scores[i,1] <- mapping_symbole_to_entrez[TF_scores[i,1]]
}

TF_scores <- TF_scores[complete.cases(TF_scores),]
TF_scores$ID <- paste("X", TF_scores$ID, sep = "")

##Combine TF and kinase
names(TF_scores) <- c("X1","NES")
kinase_activities <- as.data.frame(rbind(kinase_activities, TF_scores))

#### Input formatting

metab_input_carnival <- ttop_tumour_vs_healthy[ttop_tumour_vs_healthy$P.Value < 0.05,c(1,4)]
metab_input_carnival <- metab_input_carnival[metab_input_carnival$pubchem %in% meta_network_with_X$source | metab_input_carnival$pubchem %in% meta_network_with_X$target,]
metabs <- metab_input_carnival$pubchem
metab_input_carnival <- as.data.frame(sign(t(metab_input_carnival[,2])))
names(metab_input_carnival) <- metabs

signaling_input_carnival <- kinase_activities[abs(kinase_activities$NES) > 1.7,]
signaling_input_carnival <- signaling_input_carnival[signaling_input_carnival$X1 %in% meta_network_with_X$source | signaling_input_carnival$X1 %in% meta_network_with_X$target,]
signaling_enzymes <- signaling_input_carnival$X1
signaling_input_carnival <- as.data.frame(sign(t(signaling_input_carnival[,2])))
names(signaling_input_carnival) <- signaling_enzymes

write_csv(metab_input_carnival, "results/metabolomic/metab_input_carnival.csv")
write_tsv(metab_input_carnival, "results/metabolomic/metab_input_carnival.tsv")

write_csv(signaling_input_carnival, "results/phospho/signaling_input_carnival.csv")
write_tsv(signaling_input_carnival, "results/phospho/signaling_input_carnival.tsv")

meta_network_with_X$source <- gsub("[-+{},;() ]","______",meta_network_with_X$source)
meta_network_with_X$target <- gsub("[-+{},;() ]","______",meta_network_with_X$target)

isdup <- rep(FALSE,length(meta_network_with_X[,1]))
for(i in 1:length(meta_network_with_X[,1]))
{
  isdup[i] <- meta_network_with_X[i,1] == meta_network_with_X[i,3]
}
meta_network_with_X <- meta_network_with_X[!isdup,]
meta_network_with_X <- unique(meta_network_with_X)

write_tsv(meta_network_with_X, "support/meta_network_with_X_nobadchar.tsv")

sif_for_tiedie <- meta_network_with_X
sif_for_tiedie$interaction <- as.character(sif_for_tiedie$interaction)
sif_for_tiedie$interaction <- ifelse(sif_for_tiedie$interaction == "1","activates>","inhibits>")

write_tsv(sif_for_tiedie,"~/Documents/TieDIE/examples/kidney_cancer/pathway.sif",col_names = F)

tiedie_up <- as.data.frame(t(signaling_input_carnival))
tiedie_up$node <- row.names(tiedie_up)
tiedie_up$weight <- abs(tiedie_up$V1)*100
tiedie_up$effect <- ifelse(tiedie_up$V1 > 0,"+","-")
tiedie_up <- tiedie_up[,-1]

tiedie_down <- as.data.frame(t(metab_input_carnival))
tiedie_down$node <- row.names(tiedie_down)
tiedie_down$weight <- abs(tiedie_down$V1)*100
tiedie_down$effect <- ifelse(tiedie_down$V1 > 0, "+", "-")
tiedie_down <- tiedie_down[,-1]

write_tsv(tiedie_up,"~/Documents/TieDIE/examples/kidney_cancer/upstream.input",col_names = F)
write_tsv(tiedie_down,"~/Documents/TieDIE/examples/kidney_cancer/downstream.input",col_names = F)
