library(readr)
library(vsn)
library(limma)
library(viper)

setwd("~/Dropbox/COSMOS/")
source("scripts/limma_functions.R")
source("scripts/viper_functions.R")

psites_1 <- as.data.frame(read_csv("data/phosphoproteomic_raw.csv"))
targets_1 <- as.data.frame(read_csv("support/phosphoproteomic_targets.csv"))
row.names(psites_1) <- psites_1$ID
psites_1 <- psites_1[,-1]

### NEED TO NORMALISE PER BATCH
sub_df1 <- psites_1[,c("38KI","38TU","15KI","15TU","29KI","29TU","16KI","16TU","32KI","32TU","35KI","35TU")]
sub_df2 <- psites_1[,c("40KI","40TU","24KI","24TU","11KI","11TU")]
sub_df1 <- sub_df1[rowSums(is.na(sub_df1)) < (length(sub_df1[1,])-1),]
sub_df2 <- sub_df2[rowSums(is.na(sub_df2)) < (length(sub_df2[1,])-1),]

sub_df_list <- list(sub_df1,sub_df2)

for(i in 1:length(sub_df_list))
{
  df <- sub_df_list[[i]]
  fit <- vsnMatrix(as.matrix(df))
  meanSdPlot(fit)
  df <- as.data.frame(vsn::predict(fit,as.matrix(df)))
  sub_df_list[[i]] <- df
}

psites_1_vsn <- as.data.frame(do.call(merge,c(sub_df_list, all = T, by = "row.names")))
row.names(psites_1_vsn) <- psites_1_vsn$Row.names
psites_1_vsn <- psites_1_vsn[,-1]
psites_1_vsn <- psites_1_vsn[rowSums(is.na(psites_1_vsn)) < (length(psites_1_vsn[1,])-3),]

psites_1_vsn_bcorrect <- as.data.frame(removeBatchEffect(psites_1_vsn, batch = c(rep("A",2),rep("B",6), rep("C",4), rep("D",6)))) #Based on correlation heatmap and PCA

targets_1 <- targets_1[targets_1$sample %in% names(psites_1),]

limmaRes <- runLimma(psites_1_vsn_bcorrect, targets_1, comparisons = list(c(2,-1)))

ttop <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = 14243, adjust.method = "fdr"))

write_csv(ttop, "results/phospho/ttop_tumorVsHealthy.csv")


############

omnipath_ptm <- as.data.frame(read_csv("support/omnipath_ptm_20200205.csv"))

omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation","phosphorylation"),]
KSN <- omnipath_ptm[,c(4,3)]
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_type, sep ="_")
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_offset, sep = "")
KSN$sign <- ifelse(omnipath_ptm$modification == "phosphorylation", 1, -1)

KSN_viper <- df_to_viper_regulon(KSN)

#Prepare the measurments for viper
viper_expression <- ttop$t
names(viper_expression) <- ttop$ID
###This is just to handle a bug in viper
viper_expression <- as.data.frame(viper_expression)
viper_expression$dummy <- viper_expression$viper_expression
###COmpute TF activity scores
Kin_activity <- as.data.frame(viper(eset = viper_expression, regulon = KSN_viper, minsize = 5, adaptive.size = F, eset.filter = F))
kinases <- row.names(Kin_activity) 
Kin_activity <- as.data.frame(Kin_activity[,-1])
row.names(Kin_activity) <- kinases
names(Kin_activity) <- "NES"

write.csv(Kin_activity, "results/phospho/kinase_activities.csv")
