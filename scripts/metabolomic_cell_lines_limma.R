library(readr)
library(limma)
library(vsn)
library(pheatmap)

setwd("~/Dropbox/COSMOS/")
source("scripts/limma_functions.R")

metabolomic_cell_inhibition <- as.data.frame(read_delim("data/metabolomic_cell_inhibition.csv", 
 ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
 grouping_mark = "."), trim_ws = TRUE))

samples <- metabolomic_cell_inhibition$sample
samples <- paste(samples,rep(c(1:3),8),sep = "_")

metabolomic_cell_inhibition <- metabolomic_cell_inhibition[,-1]
metabolomic_cell_inhibition <- as.data.frame(t(metabolomic_cell_inhibition))
names(metabolomic_cell_inhibition) <- samples

targets <- as.data.frame(matrix(NA,length(metabolomic_cell_inhibition[1,]),2))
names(targets) <- c("sample","condition")
targets$sample <- names(metabolomic_cell_inhibition)
targets$condition <- gsub("_[0-9]","",targets$sample)

metabolomic_cell_inhibition[metabolomic_cell_inhibition == 0] <- NA 

metabolomic_cell_inhibition <- log2(metabolomic_cell_inhibition)

# magicPlotMaker(metabolomic_cell_inhibition,"~/Desktop/test/",targets)

comparisons <- list(c(5,-1),
                    c(2,-1),
                    c(3,-1),
                    c(4,-1),
                    c(6,-5),
                    c(7,-5),
                    c(8,-5))

limmaRes <- runLimma(metabolomic_cell_inhibition, targets, comparisons = comparisons)
ttop_786_vs_HK2 <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = length(metabolomic_cell_inhibition[,1]), adjust.method = "fdr"))
ttop_ROS_HK2 <- ttopFormatter(topTable(limmaRes[[1]], coef = 2, number = length(metabolomic_cell_inhibition[,1]), adjust.method = "fdr"))
ttop_FIL_HK2 <- ttopFormatter(topTable(limmaRes[[1]], coef = 3, number = length(metabolomic_cell_inhibition[,1]), adjust.method = "fdr"))
ttop_GSK_HK2 <- ttopFormatter(topTable(limmaRes[[1]], coef = 4, number = length(metabolomic_cell_inhibition[,1]), adjust.method = "fdr"))
ttop_ROS_786 <- ttopFormatter(topTable(limmaRes[[1]], coef = 5, number = length(metabolomic_cell_inhibition[,1]), adjust.method = "fdr"))
ttop_FIL_786 <- ttopFormatter(topTable(limmaRes[[1]], coef = 6, number = length(metabolomic_cell_inhibition[,1]), adjust.method = "fdr"))
ttop_GSK_786 <- ttopFormatter(topTable(limmaRes[[1]], coef = 7, number = length(metabolomic_cell_inhibition[,1]), adjust.method = "fdr"))

t_table <- merge(ttop_786_vs_HK2[,c(1,4)],ttop_ROS_HK2[,c(1,4)], by = "ID") 
t_table <- merge(t_table,ttop_ROS_786[,c(1,4)], by = "ID")
t_table <- merge(t_table,ttop_FIL_HK2[,c(1,4)], by = "ID")
t_table <- merge(t_table,ttop_FIL_786[,c(1,4)], by = "ID")
t_table <- merge(t_table,ttop_GSK_HK2[,c(1,4)], by = "ID")
t_table <- merge(t_table,ttop_GSK_786[,c(1,4)], by = "ID")

names(t_table) <- c("ID","786DMSO_vs_HK2DMSO","HK2ROS_vs_HK2DMSO","786ROS_vs_786DMSO","HK2FIL_vs_HK2DMSO","786FIL_vs_786DMSO","HK2GSK_vs_HK2DMSO","786GSK_vs_786DMSO")
t_table$ID <- tolower(t_table$ID)

ttop_plasmax_2_786 <- as.data.frame(
  read_csv("~/Dropbox/marco_metactivity/results/ttop_plasmax_2_786.csv"))
ttop_plasmax_2_786$ID <- tolower(ttop_plasmax_2_786$ID)
t_table <- merge(ttop_plasmax_2_786[,c(1,4)],t_table, by = "ID", all = T)
names(t_table)[2] <- c("786_vs_HK2_PLASMAX")


ttop_tumour_vs_healthy <- as.data.frame(
  read_csv("results/metabolomic/ttop_tumour_vs_healthy.csv"))
ttop_tumour_vs_healthy$ID <- tolower(ttop_tumour_vs_healthy$ID)

t_table <- merge(ttop_tumour_vs_healthy[,c(1,4)],t_table, by = "ID", all = T)
names(t_table)[2] <- c("tumor_tissue")
row.names(t_table) <- t_table$ID
t_table <- t_table[,-1]

t_table[is.na(t_table)] <- 0
t_table[abs(t_table) <= 2] <- 0
t_table[t_table >= 10] <- 10
t_table[t_table <= -10] <- -10
t_table <- t_table[rowMax(as.matrix(abs(t_table))) >= 2,]



pheatmap(t_table, cluster_cols = F)

t_table_reduced <- t_table[t_table$tumor_tissue != 0,]

pheatmap(t_table_reduced, cluster_cols = F)

comparisons_2 <- list(c(5,-1),
                    c(6,-2),
                    c(7,-3),
                    c(8,-4))

limmaRes <- runLimma(metabolomic_cell_inhibition, targets, comparisons = comparisons_2)
ttop_786_vs_HK2 <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = length(metabolomic_cell_inhibition[,1]), adjust.method = "fdr"))
ttop_786ROS_vs_HK2ROS <- ttopFormatter(topTable(limmaRes[[1]], coef = 2, number = length(metabolomic_cell_inhibition[,1]), adjust.method = "fdr"))
ttop_786FIL_vs_HK2FIL <- ttopFormatter(topTable(limmaRes[[1]], coef = 3, number = length(metabolomic_cell_inhibition[,1]), adjust.method = "fdr"))
ttop_786GSK_vs_HK2GSK <- ttopFormatter(topTable(limmaRes[[1]], coef = 4, number = length(metabolomic_cell_inhibition[,1]), adjust.method = "fdr"))

t_table <- merge(ttop_786_vs_HK2[,c(1,4)],ttop_786ROS_vs_HK2ROS[,c(1,4)], by = "ID") 
t_table <- merge(t_table,ttop_786FIL_vs_HK2FIL[,c(1,4)], by = "ID")
t_table <- merge(t_table,ttop_786GSK_vs_HK2GSK[,c(1,4)], by = "ID")


names(t_table) <- c("ID","786_vs_HK2","786ROS_vs_HK2ROS","786FIL_vs_HK2FI","786GSK_vs_HK2GSK")
t_table$ID <- tolower(t_table$ID)

ttop_plasmax_2_786 <- as.data.frame(
  read_csv("~/Dropbox/marco_metactivity/results/ttop_plasmax_2_786.csv"))
ttop_plasmax_2_786$ID <- tolower(ttop_plasmax_2_786$ID)
t_table <- merge(ttop_plasmax_2_786[,c(1,4)],t_table, by = "ID", all.y = T)
names(t_table)[2] <- c("786_vs_HK2_PLASMAX")


ttop_tumour_vs_healthy <- as.data.frame(
  read_csv("results/metabolomic/ttop_tumour_vs_healthy.csv"))
ttop_tumour_vs_healthy$ID <- tolower(ttop_tumour_vs_healthy$ID)

t_table <- merge(ttop_tumour_vs_healthy[,c(1,4)],t_table, by = "ID", all = T)
names(t_table)[2] <- c("tumor_tissue")
row.names(t_table) <- t_table$ID
t_table <- t_table[,-1]

t_table[is.na(t_table)] <- 0
t_table[abs(t_table) <= 2] <- 0
t_table[t_table >= 10] <- 10
t_table[t_table <= -10] <- -10
t_table <- t_table[rowMax(as.matrix(abs(t_table))) >= 2,]



pheatmap(t_table, cluster_cols = F)

t_table_reduced <- t_table[t_table$tumor_tissue != 0,]

pheatmap(t_table_reduced, cluster_cols = F)

comparisons_3 <- list(c(5,-1),
                      c(6,-2),
                      c(7,-3),
                      c(8,-4))
