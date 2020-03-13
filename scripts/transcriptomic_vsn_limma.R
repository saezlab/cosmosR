library(vsn)
library(limma)
library(readr)

setwd("~/Dropbox/COSMOS/")
source("scripts/limma_functions.R")

count_df_integrated <- as.data.frame(read_csv("data/transcriptomic_raw.csv"))
targets <- as.data.frame(read_csv("support/transcriptomic_targets.csv"))

row.names(count_df_integrated) <- count_df_integrated$gene
count_df_integrated <- count_df_integrated[,-1]

count_df_integrated[is.na(count_df_integrated)] <- 0
count_df_integrated <- count_df_integrated[rowMeans(count_df_integrated) > 50,]

count_df_integrated[count_df_integrated == 0] <- 0.5

fit <- vsnMatrix(as.matrix(count_df_integrated))
meanSdPlot(fit)
count_df_integrated_vsn <- as.data.frame(vsn::predict(fit,as.matrix(count_df_integrated)))

limmaRes <- runLimma(count_df_integrated_vsn, targets, comparisons = list(c(2,-1)))

ttop <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = 15919, adjust.method = "fdr"))

write_csv(ttop, "results/transcriptomic/ttop_tumorvshealthy.csv")

counts_to_write <- count_df_integrated_vsn

RNAseq_entrez_to_symbol <- as.data.frame(read_delim("support/RNAseq_entrez_to_symbol", 
                                                    "\t", escape_double = FALSE, col_types = cols(`yourlist:M20191127216DA2B77BFBD2E6699CA9B6D1C41EB259129CL` = col_character()), 
                                                    trim_ws = TRUE)) #from uniprot 20191127

mapping_vec <- gsub(" .*","",RNAseq_entrez_to_symbol$`Gene names`)
names(mapping_vec) <- RNAseq_entrez_to_symbol$`yourlist:M20191127216DA2B77BFBD2E6699CA9B6D1C41EB259129CL`
counts_to_write <- counts_to_write[row.names(counts_to_write) %in% names(mapping_vec),]
for(i in 1:length(counts_to_write[,1]))
{
  row.names(counts_to_write)[i] <- mapping_vec[row.names(counts_to_write)[i]] 
}

counts_to_write$ID <- row.names(counts_to_write)
counts_to_write <- counts_to_write[,c(23,1:22)]

write_csv(counts_to_write, "data/counts_vsn.csv")
