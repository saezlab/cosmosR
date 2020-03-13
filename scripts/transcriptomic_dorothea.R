library(readr)
library(viper)

setwd("~/Dropbox/COSMOS/")
source("scripts/limma_functions.R")
source("scripts/viper_functions.R")

dorothea <- as.data.frame(read_csv("support/DOROTHEA_20200205.csv"))
dorothea <- dorothea[,c(4,3,6,7)]
dorothea$sign <- dorothea$is_stimulation - dorothea$is_inhibition
dorothea$sign <- ifelse(dorothea$sign == 0, 1, dorothea$sign)
dorothea <- dorothea[,c(1,2,5)]
dorothea <- unique(dorothea)

dorothea_viper <- df_to_viper_regulon(dorothea)

ttop_tumorvshealthy <- as.data.frame(
  read_csv("results/transcriptomic/ttop_tumorvshealthy.csv"))

RNAseq_entrez_to_symbol <- as.data.frame(read_delim("support/RNAseq_entrez_to_symbol", 
                                                    "\t", escape_double = FALSE, col_types = cols(`yourlist:M20191127216DA2B77BFBD2E6699CA9B6D1C41EB259129CL` = col_character()), 
                                                    trim_ws = TRUE)) #from uniprot 20191127

names(RNAseq_entrez_to_symbol)[1] <- "ID"

ttop_tumorvshealthy <- merge(ttop_tumorvshealthy, RNAseq_entrez_to_symbol[,c(1,6)])
ttop_tumorvshealthy <- ttop_tumorvshealthy[,c(8,2:7)]
names(ttop_tumorvshealthy)[1] <- "ID"
ttop_tumorvshealthy$ID <- gsub(" .*","",ttop_tumorvshealthy$ID)
ttop_tumorvshealthy <- unique(ttop_tumorvshealthy)

eset <- ttop_tumorvshealthy$t
names(eset) <- ttop_tumorvshealthy$ID

viperRes <- as.data.frame(viper(eset = eset, regulon = dorothea_viper, minsize = 25, adaptive.size = F, eset.filter = F, pleiotropy = T))
viperRes$TF <- row.names(viperRes)

viperRes <- viperRes[,c(2,1)]
names(viperRes) <- c("ID","NES")

write_csv(viperRes, "results/transcriptomic/TF_scores.csv")


