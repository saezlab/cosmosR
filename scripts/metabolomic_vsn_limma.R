library(readr)
library(limma)
library(vsn)

setwd("~/Dropbox/COSMOS/")
source("scripts/limma_functions.R")

targets_1 <- as.data.frame(read_csv("support/metabolomic_targets.csv"))
raw_metabolomic <- as.data.frame(read_csv("data/raw_metabolomic.csv"))
row.names(raw_metabolomic) <- raw_metabolomic[,1]
raw_metabolomic <- raw_metabolomic[,-1]

fit <- vsnMatrix(as.matrix(raw_metabolomic))
meanSdPlot(fit)
raw_metabolomic_vsn <- as.data.frame(vsn::predict(fit,as.matrix(raw_metabolomic)))

limmaRes <- runLimma(raw_metabolomic_vsn, targets_1, comparisons = list(c(1,-2)))
ttop_tumour_vs_healthy <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = length(raw_metabolomic_vsn[,1]), adjust.method = "fdr"))

write_csv(ttop_tumour_vs_healthy, "results/metabolomic/ttop_tumour_vs_healthy.csv")
