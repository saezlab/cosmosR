library(readr)
library(ggplot2)

setwd("~/Dropbox/COSMOS/")
##Dorothea/viper
dorothea <- as.data.frame(read_csv("support/DOROTHEA_20200205.csv"))
dorothea <- dorothea[,c(4,3,6,7)]
dorothea$sign <- dorothea$is_stimulation - dorothea$is_inhibition
dorothea$sign <- ifelse(dorothea$sign == 0, 1, dorothea$sign)
dorothea <- dorothea[,c(1,2,5)]

ttop_tumorvshealthy <- as.data.frame(read_csv("results/transcriptomic/ttop_tumorvshealthy.csv"))

RNAseq_entrez_to_symbol <- as.data.frame(read_delim("support/RNAseq_entrez_to_symbol", 
                                                    "\t", escape_double = FALSE, col_types = cols(`yourlist:M20191127216DA2B77BFBD2E6699CA9B6D1C41EB259129CL` = col_character()), 
                                                    trim_ws = TRUE)) #from uniprot 20191127

names(RNAseq_entrez_to_symbol)[1] <- "ID"

ttop_tumorvshealthy <- merge(ttop_tumorvshealthy, RNAseq_entrez_to_symbol[,c(1,6)])
ttop_tumorvshealthy <- ttop_tumorvshealthy[,c(8,2:7)]
names(ttop_tumorvshealthy)[1] <- "ID"
ttop_tumorvshealthy$ID <- gsub(" .*","",ttop_tumorvshealthy$ID)


subnet_sif3_wTF <- as.data.frame(read_csv("results/multi_omic/carnival/subnet_sif3_wTF.csv"))
subnet_sif_full_att_wTF <- as.data.frame(read_csv("results/multi_omic/carnival/subnet_sif_full_att_wTF.csv"))
subnet_sif_full_att2_wTF <- as.data.frame(read_csv("results/multi_omic/carnival/subnet_sif_full_att2_wTF.csv"))

subnet_sif_full_att <- as.data.frame(rbind(subnet_sif_full_att_wTF, subnet_sif_full_att2_wTF))
subnet_sif_full_att <- unique(subnet_sif_full_att)

TPR_all <- c()
TPR_restricted <- c()

n_all <- c()
n_restricted <- c()
baseline <- c()
for( i in seq(0,2,0.01))
{
  TF_carnival_net <- subnet_sif_full_att[subnet_sif_full_att$type == "TF",]
  TF_targets_carnival_net <- subnet_sif_full_att[subnet_sif_full_att$Nodes %in% dorothea$target_genesymbol,]
  
  TF_targets_carnival_net$true_carnival_target <- unlist(lapply(TF_targets_carnival_net$Nodes,function(x, TF_carnival_net) {
    ifelse(sum(subnet_sif3_wTF[subnet_sif3_wTF$Node2 == x,"Node1"] %in% TF_carnival_net$Nodes), T, F) #If the target is directly downstream of a TF in carnival network
  },TF_carnival_net = TF_carnival_net))
  
  TF_targets_carnival_net <- TF_targets_carnival_net[TF_targets_carnival_net$true_carnival_target,]
  
  TF_targets_carnival_net$TF_is_dorothea <- unlist(lapply(TF_targets_carnival_net$Nodes,function(x, TF_carnival_net) {
    ifelse(sum(subnet_sif3_wTF[subnet_sif3_wTF$Node2 == x,"Node1"] %in% TF_carnival_net[TF_carnival_net$measured == 1, 1]), T, F) #If the target is directly downstream of a TF in carnival network
  },TF_carnival_net = TF_carnival_net))
  
  TF_targets_carnival_net$RNA_t_val <- unlist(lapply(TF_targets_carnival_net$Nodes, function(x, ttop_tumorvshealthy) {
    if(x %in% ttop_tumorvshealthy$ID)
    {
      return(ttop_tumorvshealthy[ttop_tumorvshealthy$ID == x,"t"])
    } else 
    {
      return(NA)
    }
    
  }, ttop_tumorvshealthy = ttop_tumorvshealthy))
  
  TF_targets_carnival_net <- TF_targets_carnival_net[complete.cases(TF_targets_carnival_net),]
  TF_targets_carnival_net$RNA_t_val <- ifelse(abs(TF_targets_carnival_net$RNA_t_val) > i, TF_targets_carnival_net$RNA_t_val, 0)
  row.names(TF_targets_carnival_net) <- 1:length(TF_targets_carnival_net[,1])
  
  TPR_all <- c(TPR_all,sum(sign(TF_targets_carnival_net$Activity) == sign(TF_targets_carnival_net$RNA_t_val)) / length(TF_targets_carnival_net[,1])) #TPR
  sum(sign(TF_targets_carnival_net$Activity) != sign(TF_targets_carnival_net$RNA_t_val)) / length(TF_targets_carnival_net[,1]) #FPR
  
  TF_targets_carnival_net_not_dorothea <- TF_targets_carnival_net[!(TF_targets_carnival_net$TF_is_dorothea),]
  TPR_restricted <- c(TPR_restricted,sum(sign(TF_targets_carnival_net_not_dorothea$Activity) == sign(TF_targets_carnival_net_not_dorothea$RNA_t_val)) / length(TF_targets_carnival_net_not_dorothea[,1])) #TPR
  
  n_all <- c(n_all,length(TF_targets_carnival_net[,1]))
  n_restricted <- c(n_restricted,length(TF_targets_carnival_net_not_dorothea[,1]))
  
  tvals <- TF_targets_carnival_net_not_dorothea$RNA_t_val
  tvals <- ifelse(abs(tvals) > i, tvals,0)
  baseline <- c(baseline,(sum(sign(tvals) == 1)/length(tvals) + sum(sign(tvals) == -1)/length(tvals))/2)
}

TPR_all <- as.data.frame(cbind(seq(0,2,0.01),TPR_all))
TPR_all$group <- "all"
TPR_restricted <- as.data.frame(cbind(seq(0,2,0.01),TPR_restricted))
TPR_restricted$group <- "CARNIVAL"
names(TPR_all)[2] <- "TPR" 
names(TPR_restricted)[2] <- "TPR"

TPR <- TPR_restricted
names(TPR)[1] <- "threshold"
TPR$n <- n_restricted

baseline <- as.data.frame(baseline)
baseline$threshold <- seq(0,2,0.01)
baseline$n <- n_restricted
baseline$group <- "baseline"
names(baseline)[1]  <- "TPR"
baseline <- baseline[,c(2,1,4,3)]

TPR <- as.data.frame(rbind(TPR,baseline))

ggplot(TPR, aes(x = threshold, y = TPR, group = group, color = group)) + geom_line() + theme_minimal()

###########################
cor_net_tumor <- as.data.frame(read_csv("support/cor_net_tumor.csv"))

cor_net <- cor_net_tumor

corregs_carnival <- list()
k <- 1
for(node in unique(subnet_sif3_wTF$Node1))
{
  df <- subnet_sif3_wTF[subnet_sif3_wTF$Node1 == node,]
  df <- df[df$Node2 %in% subnet_sif_full_att[subnet_sif_full_att$measured == 1,1],]
  row_list <- list()
  l <- 1
  if(length(df[,1]) > 1)
  {
    print(k)
    for(i in 1:length(df[,1]))
    {
      for(j in i:length(df[,1]))
      {
        row_list[[l]] <- c(df[i,3],df[j,3],df[i,2]*df[j,2])
        l <- l+1
      }
    }
    corregs_carnival[[k]] <- as.data.frame(do.call(rbind,row_list))
    names(corregs_carnival)[k] <- node
    k <- k+1
  }
}
corregs_carnival <- as.data.frame(do.call(rbind,corregs_carnival))
corregs_carnival <- corregs_carnival[corregs_carnival$V1 != corregs_carnival$V2,]
row.names(corregs_carnival) <- 1:length(corregs_carnival[,1])
corregs_carnival_inverse <- corregs_carnival
names(corregs_carnival_inverse) <- names(corregs_carnival_inverse)[c(2,1,3)]
corregs_carnival <- as.data.frame(rbind(corregs_carnival,corregs_carnival_inverse))
corregs_carnival$ID <- paste(corregs_carnival[,1],corregs_carnival[,2], sep = "_")

cor_net$ID <- paste(cor_net[,1],cor_net[,2], sep = "_")

comparison_regnets <- merge(cor_net[,c(4,3)], corregs_carnival[c(4,3)], by = "ID")
names(comparison_regnets) <- c("ID","correlation","carnival_correlation")
comparison_regnets$correlation <- as.numeric(as.character(comparison_regnets$correlation))
comparison_regnets$carnival_correlation <- as.numeric(as.character(comparison_regnets$carnival_correlation))
comparison_regnets <- unique(comparison_regnets)
comparison_regnets$correlation <- ifelse(abs(comparison_regnets$correlation) > 0.4, comparison_regnets$correlation, 0) 

sum(sign(comparison_regnets$correlation) == sign(comparison_regnets$carnival_correlation)) / length(comparison_regnets[,1])
mean(abs(comparison_regnets[(sign(comparison_regnets$correlation) != sign(comparison_regnets$carnival_correlation)),2]))

TPR_coreg <- list()
baseline_coreg <- list()
j <- 1
for( i in seq(0,1,0.01))
{
  comparison_regnets <- merge(cor_net[,c(4,3)], corregs_carnival[c(4,3)], by = "ID")
  names(comparison_regnets) <- c("ID","correlation","carnival_correlation")
  comparison_regnets$correlation <- as.numeric(as.character(comparison_regnets$correlation))
  comparison_regnets$carnival_correlation <- as.numeric(as.character(comparison_regnets$carnival_correlation))
  comparison_regnets <- unique(comparison_regnets)
  comparison_regnets$correlation <- ifelse(abs(comparison_regnets$correlation) > i, comparison_regnets$correlation, 0) 
  
  TPR_coreg[[j]] <- c(i,length(comparison_regnets$ID),sum(sign(comparison_regnets$correlation) == sign(comparison_regnets$carnival_correlation)) / length(comparison_regnets[,1]))
  baseline_coreg[[j]] <- c(i,length(comparison_regnets$ID),(sum(sign(comparison_regnets$correlation) == 1)/length(comparison_regnets$correlation) + sum(sign(comparison_regnets$correlation) == -1)/length(comparison_regnets$correlation))/2)
  j <- j+1
}

TPR_coreg <- as.data.frame(do.call(rbind, TPR_coreg))
names(TPR_coreg) <- c("threshold","n","TPR")
TPR_coreg$group <- "CARNIVAL"

baseline_coreg <- as.data.frame(do.call(rbind, baseline_coreg))
names(baseline_coreg) <- c("threshold","n","TPR")
baseline_coreg$group <- "baseline"

TPR_coreg <- as.data.frame(rbind(TPR_coreg,baseline_coreg))

ggplot(TPR_coreg, aes(x = threshold, y = TPR, group = group, color = group)) + geom_line() + theme_minimal()
