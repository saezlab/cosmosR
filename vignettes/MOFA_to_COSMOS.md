## MOFA-COSMOS pipeline

### Installation and dependency

In this tutorial, we use the new Meta-Footprint Analysis MOON to score a
mechanistic network.

``` r
#install cosmosR from github directly to obtain the newest version
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("cosmosR", quietly = TRUE)) devtools::install_github("saezlab/cosmosR")
```

DecoupleR is used to score various regulatory interactions.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("decoupleR", quietly = TRUE)) BiocManager::install("saezlab/decoupleR")
```

We are using the LIANA package to process the LR interactions.

``` r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("liana", quietly = TRUE)) remotes::install_github('saezlab/liana')
```

We are using MOFA2 (Argelaguet et al., 2018) to find decompose variance
across different omics and use COSMOS to find mechanistic
interpretations of its factors. However, since we are using the python
version of MOFA2, please make sure to have
[mofapy2](https://github.com/bioFAM/mofapy2) as well as panda and numpy
installed in your (conda) environment (e.g. by using reticulate:
reticulate::use_condaenv(“base”, required = T) %\>%
reticulate::conda_install(c(“mofapy2”,“panda”,“numpy”)). For downstream
analysis, the R version of MOFA2 is used.

``` r
# Install MOFA2 from bioconductor (stable version)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("MOFA2", quietly = TRUE)) BiocManager::install("MOFA2")
```

General information (tutorials) on how to use COSMOS and MOFA2, see
[COSMOS](https://saezlab.github.io/cosmosR/articles/tutorial.html) and
[MOFA2](https://biofam.github.io/MOFA2/tutorials.html).

For data manipulation, visualization, etc. further packages are needed
which are loaded here along with MOFA2 and COSMOS.

``` r
#imports
library(readr)

#core methods
library(cosmosR)
library(MOFA2)
library(liana)
library(decoupleR)

#data manipulations
library(dplyr)
library(reshape2)
library(GSEABase)
library(tidyr)

#plotting
library(ggplot2)
library(ggfortify)
library(pheatmap)
library(gridExtra)
library(RColorBrewer)
# library(RCy3)
```

### From omics data to MOFA ready input

Since extensive pre-processing is beyond the scope of this tutorial,
here only the transformation of pre-processed single omics data to a
long data.frame (columns: “sample”, “group”, “feature”, “view”, “value”)
is shown. For more details and information regarding appropriate input,
please refer to the MOFA tutorial [MOFA2: training a model in
R](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/getting_started_R.html).
The raw data is accessed through [NCI60
cellminer](https://discover.nci.nih.gov/cellminer/home.do) and the
scripts for pre-processing can be found in the associated folder.

Here, we have the following views (with samples as columns and genes as
rows): transcriptomics (RNA-seq
[PMID31113817](https://pubmed.ncbi.nlm.nih.gov/31113817/)),
proteomics (SWATH (Mass spectrometry)
[PMID31733513](https://pubmed.ncbi.nlm.nih.gov/31733513/%3E)} and
metabolomics (LC/MS & GC/MS (Mass spectrometry) [DTP NCI60
data](https://wiki.nci.nih.gov/display/NCIDTPdata/Molecular+Target+Data)).
The group information is stored inside the RNA metadata and based on
transcription factor clustering (see scripts/RNA/Analyze_RNA.R).
Further, only samples are kept for which each omic is available (however
this is not necessary since MOFA2 is able to impute data).

``` r
# Transcriptomics
## Load data
RNA <- as.data.frame(read_csv("data/RNA/RNA_log2_FPKM_clean.csv"))
rownames(RNA) <- RNA[,1]
RNA <- RNA[,-1]

## Remove genes with excessive amount of NAs (only keep genes with max. amount 
#of NAs = 33.3 % across cell lines)
RNA <- RNA[rowSums(is.na(RNA))<(dim(RNA)[2]/3),]

## Only keep highly variable genes (as suggested by MOFA): here top ≈70% genes 
#to keep >75% of TFs (103 out of 133). For stronger selection, we can further 
#reduce the number of genes (e.g. top 50%)
RNA_sd <- sort(apply(RNA, 1, function(x) sd(x,na.rm = T)), decreasing = T)
hist(RNA_sd, breaks = 100)
```

![](MOFA_to_COSMOS_files/figure-gfm/Preparing%20MOFA%20input-1.png)<!-- -->

``` r
hist(RNA_sd[1:6000], breaks = 100)
```

![](MOFA_to_COSMOS_files/figure-gfm/Preparing%20MOFA%20input-2.png)<!-- -->

``` r
RNA_sd <- RNA_sd[1:6000]
RNA <- RNA[names(RNA_sd),]

# Proteomics
## Load data
proteo <- as.data.frame(read_csv("data/proteomic/Prot_log10_SWATH_clean.csv"))
rownames(proteo) <- proteo[,1]
proteo <- proteo[,-1]

## Only keep highly variable proteins (as suggested by MOFA): here top ≈60%. 
#For stronger selection, we can further reduce the number of proteins (e.g. only keep top 50%)
proteo_sd <- sort(apply(proteo, 1, function(x) sd(x,na.rm = T)), decreasing = T)
hist(proteo_sd, breaks = 100)
```

![](MOFA_to_COSMOS_files/figure-gfm/Preparing%20MOFA%20input-3.png)<!-- -->

``` r
hist(proteo_sd[1:round(dim(proteo)[1]*3/5)], breaks = 100)
```

![](MOFA_to_COSMOS_files/figure-gfm/Preparing%20MOFA%20input-4.png)<!-- -->

``` r
proteo_sd <- proteo_sd[1:round(dim(proteo)[1]*3/5)]
proteo <- proteo[names(proteo_sd),]

# Metabolomics
## Load data
metab <- as.data.frame(read_csv("data/metabolomic/metabolomic_clean_vsn.csv"))
rownames(metab) <- metab[,1]
metab <- metab[,-1]

## Since we have only a limited number of metabolite measurements, 
#all measurements are kept here
metab_sd <- apply(metab, 1, function(x) sd(x,na.rm=T))
hist(metab_sd, breaks = 100)
```

![](MOFA_to_COSMOS_files/figure-gfm/Preparing%20MOFA%20input-5.png)<!-- -->

``` r
# Create long data frame
## Only keep samples with each view present 
overlap_patients <- intersect(intersect(names(RNA),names(proteo)),names(metab))
RNA <- RNA[,overlap_patients]
proteo <- proteo[,overlap_patients]
metab <- metab[,overlap_patients]

## Create columns required for MOFA
### RNA
RNA <- melt(as.data.frame(cbind(RNA,row.names(RNA))))
RNA$view <- "RNA"
RNA <- RNA[,c(2,1,4,3)]                 
names(RNA) <- c("sample","feature","view","value")

### Proteomics
proteo <- melt(as.data.frame(cbind(proteo,row.names(proteo))))
proteo$view <- "proteo"
proteo <- proteo[,c(2,1,4,3)]                 
names(proteo) <- c("sample","feature","view","value")

### Metabolomics
metab <- melt(as.data.frame(cbind(metab,row.names(metab))))
metab$view <- "metab"
metab <- metab[,c(2,1,4,3)]                 
names(metab) <- c("sample","feature","view","value")

## Merge long data frame to one
mofa_ready_data <- as.data.frame(do.call(rbind,list(RNA,proteo,metab)))

# Save MOFA metadata information and MOFA input
write_csv(mofa_ready_data, file = "data/mofa/mofa_data.csv")
```

The final data frame has the following structure:

``` r
head(mofa_ready_data)
```

    ##   sample  feature view  value
    ## 1  786-0     KRT8  RNA  7.261
    ## 2  786-0      FN1  RNA  7.823
    ## 3  786-0   S100A4  RNA     NA
    ## 4  786-0 RNA5-8S5  RNA 11.045
    ## 5  786-0     MT1E  RNA  7.010
    ## 6  786-0      VIM  RNA  9.366

We have successfully combined single omic data sets to a long
multi-omics data frame that can be used as a MOFA input.

In this part we assess the correlation between the RNA and proteomic
data. We both look for correlation across samples AND across features.
When we check across sample, we interrogate the correlation at the level
of single genes, that is for given gene, how is a change is transcript
abundance is reflected in the the level of its protein. This correlation
can be different between genes, as there are all subject to different
regulatory mechanisms. when we check across features, we interrogate the
correlation at the level of single samples. that is, we check if overall
the protein abundance are good proxy of transcript abundance, and
vice-versa. The imperfect correlation between trnaxcripts and proteins
reflects the complexity of the translation process, as well as protein
and RNA stability.

``` r
condition_prot_RNA_correlation <- sapply(unique(mofa_ready_data$sample), function(condition){
  merged_data <- merge(mofa_ready_data[which(mofa_ready_data$sample == condition & mofa_ready_data$view == "RNA"),-c(1,3)],
                       mofa_ready_data[which(mofa_ready_data$sample == condition & mofa_ready_data$view == "proteo"),-c(1,3)],
                       by = "feature")
  return(cor(merged_data[,2],merged_data[,3], use = "pairwise.complete.obs"))
})

hist(condition_prot_RNA_correlation, breaks = 20)
```

![](MOFA_to_COSMOS_files/figure-gfm/Supp%20figure%20condition_prot_RNA_correlation-1.png)<!-- -->

``` r
overlap_genes <- unique(intersect(proteo$feature,RNA$feature))
single_prot_RNA_correlation <- sapply(overlap_genes, function(gene){
  merged_data <- merge(mofa_ready_data[which(mofa_ready_data$feature == gene & mofa_ready_data$view == "RNA"),-c(2,3)],
                       mofa_ready_data[which(mofa_ready_data$feature == gene & mofa_ready_data$view == "proteo"),-c(2,3)],
                       by = "sample")
  return(cor(merged_data[,2],merged_data[,3], use = "pairwise.complete.obs"))
})

hist(single_prot_RNA_correlation, breaks = 100)
```

![](MOFA_to_COSMOS_files/figure-gfm/Supp%20figure%20single_prot_RNA_correlation-1.png)<!-- -->

### MOFA run

After pre-processing the data into MOFA appropriate input, let’s perform
the MOFA analysis. To change specific MOFA options the function
write_MOFA_options_file is used. For more information regarding option
setting see [MOFA option
tutorial](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/getting_started_R.html).
Since we don’t know how many factors are needed to explain our data
well, various maximum numbers of factors are tested. The final results
can be found in the results/mofa/ folder.

Only run this chunk if you want to re-perform the mofa optimisation. The
results are already available in the results/mofa folder. We are testing
various amount of maximum number of factors that MOFA is allowed to
explore.

You need to instal mofapy2 in python3 to run this chunk. On windows, you
can simply run “python3” the terminal to make sure its installed or be
automatically redirected to the windows store page if it’s not. Then,
you can run “python3 -m pip install –user mofapy2” in a terminal to
install the module to python3.

``` r
for(i in c(5:15,20,length(unique(mofa_ready_data$sample)))){
wd <- getwd()
write_MOFA_options_file(input_file = '/data/mofa/mofa_data.csv', output_file = paste0('/results/mofa/mofa_res_',i,'factor.hdf5'), factors = i, convergence_mode = "slow", likelihoods = c('gaussian','gaussian','gaussian'))
system(paste0(paste('python3', wd, sep = " "), paste('/scripts/mofa/mofa2.py', "options_MOFA.csv")))
}
```

### Investigate MOFA output

We first load the different MOFA models containing the MOFA results
together with the metadata information. Since we don’t know a-priori
which is the ideal number of factors, we can settle on the highest
number of factor where the model stabilises.

``` r
# Load MOFA output
for(i in c(5:15,20,length(unique(mofa_ready_data$sample)))){
  filepath <- paste0('results/mofa/mofa_res_',i,'factor.hdf5')
  model <- load_model(filepath)
  meta_data <- read_csv("data/metadata/RNA_metadata_cluster.csv")[,c(1,2)]
  colnames(meta_data) <- c("sample","cluster")
  mofa_metadata <- merge(samples_metadata(model), meta_data, by = "sample")
  names(mofa_metadata) <- c("sample","group","cluster")
  samples_metadata(model) <- mofa_metadata
  assign(paste0("model",i),model)
}
```

Then, we can compare the factors from different MOFA models and check on
consistency (here: model 7-13).

``` r
compare_factors(list(model7,model8,model9,model10,model11,model12,model13), cluster_rows = F, cluster_cols = F,)
```

![](MOFA_to_COSMOS_files/figure-gfm/Compare%20MOFA%20models-1.png)<!-- -->

Here, we can see that the inferred MOFA factor weights do not change
drastically around 9. For the sake of this tutorial, we choose the
9-factor MOFA model.

``` r
model <- model10
```

The next part deals with the analysis and interpretation of the MOFA
factors with specific focus on the factor weights. The classic tutorial
on how to perform downstream analysis after MOFA model training can be
found here [MOFA+: downstream analysis (in
R)](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/downstream_analysis.html).
In this context, it is important to emphasize that the factors can be
interpreted similarly to principal components (PCs) in a principal
component analysis (PCA).

An initial metric we can explore is the total variance ($R^2$) explained
per view by our MOFA model. This helps us understand how well the model
represents the data.

``` r
# Investigate factors: Explained total variance per view
plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
```

![](MOFA_to_COSMOS_files/figure-gfm/Total%20variance%20per%20factor-1.png)<!-- -->

``` r
calculate_variance_explained(model)
```

    ## $r2_total
    ## $single_group
    ##      RNA    metab   proteo 
    ## 59.59992 18.70369 23.29474 
    ## 
    ## attr(,"class")
    ## [1] "relistable" "list"      
    ## 
    ## $r2_per_factor
    ## $single_group
    ##                RNA        metab       proteo
    ## Factor1 24.0201682  0.002979224  0.214741826
    ## Factor2  4.4849965 14.093737072  2.236180823
    ## Factor3  0.2777266  2.193923605 12.875906515
    ## Factor4  7.3683342  2.461682919  4.752507250
    ## Factor5 13.2330359  0.008854678  0.004884149
    ## Factor6  2.3102555  0.001537887  2.692284962
    ## Factor7  4.5839551  0.007795120  0.390045680
    ## Factor8  2.8171650  0.277818450  0.210418704
    ## Factor9  1.9469859  0.004198608  0.003992617
    ## 
    ## attr(,"class")
    ## [1] "relistable" "list"

The bar plot demonstrates that the nine factors for the view RNA explain
more than 50% of the variance. In contrast, for both metabolomics as
well as proteomics around 20% of the variance can be explained by all
factors.

It isn’t unexpected that more variance of the RNA dataset can be
reconstructed. Usually, more features available in a given view will
lead to a higher % of variance explained. It is interesting to notice
though that similar amount of variance are explained for metabolic and
proteomic views, despite having very different number of features. This
may be due to the fact that metabolite abundances are in general more
cross-correlated than protein abundances (due to the underlying
regulatory structures, e.i. reaction networks vs protein synthesis and
degratation regulatory mechanisms).

``` r
#Get the cross correlation for RNA and proteomic data
RNA <- dcast(mofa_ready_data[which(mofa_ready_data$view == "RNA"),-3], feature~sample, value.var = "value")
row.names(RNA) <- RNA[,1]
RNA <- RNA[,-1]

RNA_crosscor <- cor(t(RNA), use = "pairwise.complete.obs")
# hist(RNA_crosscor[upper.tri(RNA_crosscor)], breaks = 1000)

prot <- dcast(mofa_ready_data[which(mofa_ready_data$view == "proteo"),-3], feature~sample, value.var = "value")
row.names(prot) <- prot[,1]
prot <- prot[,-1]

prot_crosscor <- cor(t(prot), use = "pairwise.complete.obs")
# hist(prot_crosscor[upper.tri(prot_crosscor)], breaks = 1000)

metabo <- dcast(mofa_ready_data[which(mofa_ready_data$view == "metab"),-3], feature~sample, value.var = "value")
row.names(metabo) <- metabo[,1]
metabo <- metabo[,-1]

metabo_crosscor <- cor(t(metabo), use = "pairwise.complete.obs")
# hist(metabo_crosscor[upper.tri(metabo_crosscor)], breaks = 1000)

# plot(calculate_variance_explained(model)$r2_total$single_group, c(mean(RNA_crosscor),mean(metabo_crosscor),mean(prot_crosscor)))
plot(calculate_variance_explained(model)$r2_total$single_group, c(dim(RNA_crosscor)[1],dim(metabo_crosscor)[1],dim(prot_crosscor)[1]))
```

![](MOFA_to_COSMOS_files/figure-gfm/Chekc%20cross%20correlation%20of%20omics%20compared%20to%20variance%20epxlained-1.png)<!-- -->

``` r
df_variance_crosscor <- as.data.frame(cbind(calculate_variance_explained(model)$r2_total$single_group,
                                           c(mean(RNA_crosscor),mean(metabo_crosscor),mean(prot_crosscor)),
                                           c(dim(RNA_crosscor)[1],dim(metabo_crosscor)[1],dim(prot_crosscor)[1])))
names(df_variance_crosscor) <- c("var_expl","mean_crosscor","n_features")
df_variance_crosscor$prod <- sqrt(df_variance_crosscor$mean_crosscor) * df_variance_crosscor$n_features

#Just plot them as scatter plot with regression lines
mofa_solved <- ggplot(df_variance_crosscor, aes(x = prod, y = var_expl)) + 
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  theme_minimal() +
  xlab("sqrt(mean cross-corelation of view) * number of features") +
  ylab("variance explained for each view") +
  ggtitle("MOFA size and cross cor influence")

mofa_solved
```

![](MOFA_to_COSMOS_files/figure-gfm/Chekc%20cross%20correlation%20of%20omics%20compared%20to%20variance%20epxlained-2.png)<!-- -->

``` r
mofa_nfeatures <- ggplot(df_variance_crosscor, aes(x = n_features, y = var_expl)) + 
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +
  theme_minimal() +
  xlab("sqrt(mean cross-corelation of view) * number of features") +
  ylab("variance explained for each view") +
  ggtitle("MOFA size influence")

mofa_nfeatures
```

![](MOFA_to_COSMOS_files/figure-gfm/Chekc%20cross%20correlation%20of%20omics%20compared%20to%20variance%20epxlained-3.png)<!-- -->

``` r
#Get the actual coeficient signficance
summary(lm(data= df_variance_crosscor, var_expl~prod))
```

    ## 
    ## Call:
    ## lm(formula = var_expl ~ prod, data = df_variance_crosscor)
    ## 
    ## Residuals:
    ##      RNA    metab   proteo 
    ## -0.06207 -0.77553  0.83760 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept) 1.885e+01  8.542e-01   22.06   0.0288 *
    ## prod        1.434e-02  5.176e-04   27.70   0.0230 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.143 on 1 degrees of freedom
    ## Multiple R-squared:  0.9987, Adjusted R-squared:  0.9974 
    ## F-statistic: 767.1 on 1 and 1 DF,  p-value: 0.02297

``` r
summary(lm(data= df_variance_crosscor, var_expl~n_features))
```

    ## 
    ## Call:
    ## lm(formula = var_expl ~ n_features, data = df_variance_crosscor)
    ## 
    ## Residuals:
    ##    RNA  metab proteo 
    ##  1.480  3.426 -4.905 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)
    ## (Intercept) 14.366419   5.255059   2.734    0.223
    ## n_features   0.007292   0.001446   5.043    0.125
    ## 
    ## Residual standard error: 6.163 on 1 degrees of freedom
    ## Multiple R-squared:  0.9622, Adjusted R-squared:  0.9243 
    ## F-statistic: 25.43 on 1 and 1 DF,  p-value: 0.1246

Since we would like to use the found correlation for downstream COSMOS
analysis, we first investigate the variance ($R^2$) each factor can
explain per view as well as investigate the factor explaining the most
variance across all views.

``` r
# Investigate factors: Explained variance per view for each factor 
pheatmap(model@cache$variance_explained$r2_per_factor[[1]], display_numbers = T, angle_col = "0", legend_labels = c("0","10", "20", "30", "40", "Variance\n\n"), legend = T, main = "", legend_breaks = c(0,10, 20, 30, 40, max(model@cache$variance_explained$r2_per_factor[[1]])), cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white","red"))(100), fontsize_number = 10)
```

![](MOFA_to_COSMOS_files/figure-gfm/Variance%20per%20view%20per%20factor-1.png)<!-- -->

``` r
pheatmap(model@cache$variance_explained$r2_per_factor[[1]], display_numbers = T, angle_col = "0", legend_labels = c("0","10", "20", "30", "40", "Variance\n\n"), legend = T, main = "", legend_breaks = c(0,10, 20, 30, 40, max(model@cache$variance_explained$r2_per_factor[[1]])), cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white","red"))(100), fontsize_number = 10,filename = "results/mofa/variance_heatmap.pdf",width = 4, height = 2.5)
```

By analyzing this heatmap, we can observe that Factor 1 and 5 explain a
large proportion of variance of the RNA, Factor 2 and 4 mainly explains
the variance of the metabolomics and Factor 3 highlights variability
from the proteomics view. Consequently, no factor explains a large
portion of the variance across all views, but factor 2 and 4 capture
variance across all views. We can also analyze the correlation between
factors. The next plot highlights the spearman correlation between each
factor.

``` r
# Investigate factors: Correlation between factors 
pdf(file="results/mofa/plot_factor_cor.pdf", height = 5, width = 5)
plot_factor_cor(model, type = 'upper', method = "spearman", addCoef.col = "black")
```

Next, we can characterize the factors with the available metadata of
samples. Indeed, it can be interesting to see if some factors are
associated with metadata classes, as they can help to guide our
intepretation of the factors. This part rely on an interpretation of the
Z matrix of the mofa model, while ignoring the factor weight themselves.
The factor weights will be analysed subsequently.

To begin simply, we can check if some tissue of origin are associated
with any of the factors. There are different ways thsi can be done, but
we find that using a simple linear association between a given tissues
category in a factor compared to any other tissue categories
(e.g. breast tissue vs every other tissue types) provide an easy to
interpret insight that serves our purpose well. To do so, we use the
Univariate Linear Model (ULM) function of decoupleR.

``` r
NCI_60_metadata <- as.data.frame(read_csv(file = "data/metadata/RNA_metadata_cluster.csv"))
```

    ## Rows: 60 Columns: 16
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (12): cell_line, tissue of origin a, sex a, prior treatment a,b, Epithel...
    ## dbl  (4): cluster, age a, mdr f, doubling time g
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#discretize age category so that regression can be applied
NCI_60_metadata$`age over 50` <- NCI_60_metadata$`age a` > 50
tissues <-NCI_60_metadata[,c(3,1)]
names(tissues) <- c("source","target")

Z_matrix <- as.data.frame(model@expectations$Z$single_group)

#Check specifically tissue meta-data association with Z matrix
tissue_enrichment <- decoupleR::run_ulm(Z_matrix, tissues)
tissue_enrichment <- reshape2::dcast(tissue_enrichment, source~condition, value.var = "score")
row.names(tissue_enrichment) <- tissue_enrichment$source
tissue_enrichment <- tissue_enrichment[,-1]

#plot them as heatmaps
palette <- make_heatmap_color_palette(tissue_enrichment)
pheatmap(t(tissue_enrichment), show_rownames = T, cluster_cols = F, cluster_rows = F,color = palette, angle_col = 315,display_numbers = F, filename = "results/mofa/mofa_tissue_enrichment.pdf", width = 3, height = 2.6)
pheatmap(tissue_enrichment, show_rownames = T, cluster_cols = F, cluster_rows = F,color = palette, angle_col = 315,display_numbers = T)
```

Here, we can see that factor 4 seems to be strongly associated with
eukemia, and factor 2 with melanoma.

``` r
#format meta-data names
all_metadata <- reshape2::melt(NCI_60_metadata[,c(1,3,17,5,7,8,10,13,14)], id.vars = "cell_line")
all_metadata$variable <- gsub(" a$","",all_metadata$variable)
all_metadata$variable <- gsub(" a,c$","",all_metadata$variable)
all_metadata$variable <- gsub(" d$","",all_metadata$variable)
all_metadata$variable <- gsub("histology","histo",all_metadata$variable)
all_metadata$variable <- gsub("tissue of origin","tissue",all_metadata$variable)
all_metadata$value <- gsub("[(]81-103[)]","",all_metadata$value)
all_metadata$value <- gsub("[(]35-57[)]","",all_metadata$value)
all_metadata$value <- gsub("Memorial Sloan Kettering Cancer Center","Memorial SKCC",all_metadata$value)
all_metadata$value <- gsub("Malignant melanotic melanoma","Malig. melanotic melan.",all_metadata$value)
all_metadata$value <- gsub("Ductal carcinoma-mammary gland","Duct. carci-mam. gland",all_metadata$value)

all_metadata$target <- paste(all_metadata$variable, all_metadata$value, sep = ": ")
all_metadata <- all_metadata[,c(4,1)]
names(all_metadata) <- c("source","target")
  
Z_matrix <- as.data.frame(model@expectations$Z$single_group)

#run ULM on the Z matrix meta data to find associated metadata features
all_metadata_enrichment <- decoupleR::run_ulm(Z_matrix, all_metadata, minsize = 3)
all_metadata_enrichment <- reshape2::dcast(all_metadata_enrichment, source~condition, value.var = "score")
row.names(all_metadata_enrichment) <- all_metadata_enrichment$source
all_metadata_enrichment <- all_metadata_enrichment[,-1]
all_metadata_enrichment_top <- all_metadata_enrichment[
  apply(all_metadata_enrichment, 1, function(x){max(abs(x)) > 2}),
]

#Plot the heatmap
palette <- make_heatmap_color_palette(all_metadata_enrichment_top)
pheatmap(t(all_metadata_enrichment_top), show_rownames = T, cluster_cols = F, cluster_rows = F,color = palette, angle_col = 315,display_numbers = F, filename = "results/mofa/mofa_metadata_enrichment.pdf", width = 4.5, height = 4)
pheatmap(t(all_metadata_enrichment_top), show_rownames = T, cluster_cols = F, cluster_rows = F,color = palette, angle_col = 315,display_numbers = T)
```

Further, we can inspect the composition of the factors by plotting the
top 10 weights per factor and view.

``` r
## Top 20 factor weights per view
grid.arrange(
  ### RNA
  plot_top_weights(model, 
                   view = "RNA",
                   factor = 4,
                   nfeatures = 10) +
    ggtitle("RNA"),
  ### Proteomics
  plot_top_weights(model, 
                   view = "proteo",
                   factor = 4,
                   nfeatures = 10) + 
    ggtitle("Proteomics"),
  ### Metabolomics
  plot_top_weights(model, 
                   view = "metab",
                   factor = 4,
                   nfeatures = 10) + 
    ggtitle("Metabolomics"), ncol =3
)
```

![](MOFA_to_COSMOS_files/figure-gfm/Factor%20weights%20per%20view-1.png)<!-- -->

Using this visualization, features with strong association with the
factor (large absolute values) can be easily identified and further
inspected by literature investigation.

Further tools to investigate the latent factors are available inside the
MOFA framework ([MOFA+: downstream analysis (in
R)](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/downstream_analysis.html)).

### From MOFA output to TF activity, ligand-receptor activity

In the next step, using the weights of Factor 4 (highest shared $R^2$
across views), different databases (liana and dorothea) and decoupleR,
the transcription factor activities, ligand-receptor scores are
estimated.

First, we have to extract the weights from our MOFA model and keep the
weights from Factor 4 in the RNA view.

``` r
weights <- get_weights(model, views = "all", factors = "all")

save(weights, file = "results/mofa/mofa_weights.RData")
# Extract Factor with highest explained variance in RNA view
RNA <- data.frame(weights$RNA[,4]) 
row.names(RNA) <- gsub("_RNA","",row.names(RNA))
```

Then we load the consensus networks from our databases, as well as the
TF-targets regulons from collectri.

``` r
# Load LIANA (receptor and ligand) consensus network
# ligrec_ressource <- distinct(liana::decomplexify(liana::select_resource("Consensus")[[1]]))
# save(ligrec_ressource, file = "support/ligrec_ressource.RData")

#process the LR interaction prior knowledge to make suitable input for decoupleR
load(file = "support/ligrec_ressource.RData")
ligrec_geneset <- format_LR_ressource(ligrec_ressource)

# Load Dorothea (TF) network
# dorothea_df <- decoupleR::get_collectri()
# save(dorothea_df, file = "support/dorothea_df.RData")
load(file = "support/dorothea_df.RData")
```

We can display what are the top weights of the model and to which factor
they are associated.

``` r
RNA_all <- as.data.frame(weights$RNA) 
row.names(RNA_all) <- gsub("_RNA","",row.names(RNA_all))

#Extract the top weights
RNA_top <- RNA_all[order(apply(RNA_all,1,function(x){max(abs(x))}), decreasing = T)[1:25],]

#plot them as heatmaps
palette <- make_heatmap_color_palette(RNA_top)
pheatmap(RNA_top, show_rownames = T, cluster_cols = F, cluster_rows = F,color = palette, angle_col = 315, filename = "results/mofa/mofa_top_RNA.pdf", width = 4, height = 4.3)
pheatmap(RNA_top, show_rownames = T, cluster_cols = F, cluster_rows = F,color = palette, angle_col = 315)
```

We can do the same for metabolites

``` r
metabs_all <- as.data.frame(weights$metab) 
row.names(metabs_all) <- gsub("_metab","",row.names(metabs_all))

#extract top metab weights
metabs_top <- metabs_all[order(apply(metabs_all,1,function(x){max(abs(x))}), decreasing = T)[1:10],]

#plot them as heatmaps
palette <- make_heatmap_color_palette(metabs_top)
pheatmap(metabs_top, show_rownames = T, cluster_cols = F, cluster_rows = F,color = palette, angle_col = 315, filename = "results/mofa/mofa_top_metabs.pdf", width = 4, height = 2.2)
pheatmap(metabs_top, show_rownames = T, cluster_cols = F, cluster_rows = F,color = palette, angle_col = 315)
```

Then, by using decoupleR and prior knowledge networks, we calculate the
different regulatory activities inferred by the ULM approach. Depending
on the task, the minimum number of targets per source varies required to
maintain the source in the output (e.g. two targets by definition for
ligand-receptor interactions). Here can display them as heatmaps.

``` r
# Calculate regulatory activities from Receptor and Ligand network
ligrec_factors <- run_ulm(mat = as.matrix(RNA_all), network = ligrec_geneset, .source = set, .target = gene, minsize = 2) 
ligrec_factors_df <- wide_ulm_res(ligrec_factors)
ligrec_factors_df <- ligrec_factors_df[order(apply(ligrec_factors_df,1,function(x){max(abs(x))}), decreasing = T)[1:25],]


#plot them as heatmaps
palette <- make_heatmap_color_palette(ligrec_factors_df)
pheatmap(ligrec_factors_df, show_rownames = T, cluster_cols = F, cluster_rows = F,color = palette, angle_col = 315, filename = "results/mofa/mofa_top_ligrec.pdf", width = 4, height = 4.3)
pheatmap(ligrec_factors_df, show_rownames = T, cluster_cols = F, cluster_rows = F,color = palette, angle_col = 315)
```

``` r
# Calculate regulatory activities from TF network
TF_factors <- run_ulm(mat = as.matrix(RNA_all), network = dorothea_df, minsize = 10)
TF_factors <- wide_ulm_res(TF_factors)
TF_factors <- TF_factors[order(apply(TF_factors,1,function(x){max(abs(x))}), decreasing = T)[1:25],]

#Plot them as heatmaps
palette <- make_heatmap_color_palette(TF_factors)
pheatmap(TF_factors, show_rownames = T, cluster_cols = F, cluster_rows = F,color = palette, angle_col = 315, filename = "results/mofa/mofa_top_TF.pdf", width = 4, height = 4.3)
pheatmap(TF_factors, show_rownames = T, cluster_cols = F, cluster_rows = F,color = palette, angle_col = 315)
```

These activity estimations are saved in an input file for further
downstream analysis as well.

``` r
# Calculate regulatory activities from Receptor and Ligand network
ligrec_high_vs_low <- run_ulm(mat = as.matrix(RNA), network = ligrec_geneset, .source = set, .target = gene, minsize = 2) 
ligrec_high_vs_low <- ligrec_high_vs_low[ligrec_high_vs_low$statistic == "ulm",]
ligrec_high_vs_low_vector <- ligrec_high_vs_low$score
names(ligrec_high_vs_low_vector) <- ligrec_high_vs_low$source

ligrec_high_vs_low_top <- ligrec_high_vs_low[order(abs(ligrec_high_vs_low$score),decreasing = T)[1:15],c(2,4)]
ligrec_high_vs_low_top$source <- factor(ligrec_high_vs_low_top$source, levels = ligrec_high_vs_low_top$source)
ggplot(ligrec_high_vs_low_top, aes(x=source, y = score)) + geom_bar(stat= "identity",position = "dodge") + theme_minimal() + theme(axis.text.x = element_text(angle = 315, vjust = 0.5, hjust=0))
```

![](MOFA_to_COSMOS_files/figure-gfm/Activity%20estimations%20for%20factor%204-1.png)<!-- -->

``` r
# Calculate regulatory activities from TF network
TF_high_vs_low <- run_ulm(mat = as.matrix(RNA), network = dorothea_df, minsize = 10)
TF_high_vs_low <- TF_high_vs_low[TF_high_vs_low$statistic == "ulm",]
TF_high_vs_low_vector <- TF_high_vs_low$score
names(TF_high_vs_low_vector) <- TF_high_vs_low$source

# Prepare moon analysis input
TF_high_vs_low <- as.data.frame(TF_high_vs_low[,c(2,4)])
row.names(TF_high_vs_low) <- TF_high_vs_low[,1]

TF_high_vs_low_top_10 <- TF_high_vs_low[order(abs(TF_high_vs_low$score), decreasing = T)[1:10],]
TF_high_vs_low_top_10$source <- factor(TF_high_vs_low_top_10$source, levels = TF_high_vs_low_top_10$source)
ggplot(TF_high_vs_low_top_10, aes(x=source, y = score)) + geom_bar(stat= "identity",position = "dodge") + theme_minimal()
```

![](MOFA_to_COSMOS_files/figure-gfm/Activity%20estimations%20for%20factor%204-2.png)<!-- -->

``` r
# Combine results
ligrec_TF_moon_inputs <- list("ligrec" = ligrec_high_vs_low_vector,
                              "TF" = TF_high_vs_low_vector)

save(ligrec_TF_moon_inputs, file = "data/cosmos/ligrec_TF_moon_inputs.Rdata")
```

We have now successfully used the MOFA weights (Factor 4, RNA view) to
infer not only TF activities but also the ligand-receptor interactions.

### Extracting network with COSMOS to formulate mechanistic hypotheses

In this next parts the COSMOS run is going to be performed. We can use
the output of decoupleR (the activities) next to the factor weights as
an input for COSMOS to specifically focus the analysis on pre-computed
data-driven results.

Here, in the first step, the data is loaded. Besides the activity
estimations and the factor weights, the previous filtered out expressed
gene names are loaded.

``` r
# Load data
load(file = "data/cosmos/ligrec_TF_moon_inputs.Rdata")

signaling_input <- ligrec_TF_moon_inputs$TF
ligrec_input <- ligrec_TF_moon_inputs$ligrec

RNA_input <- weights$RNA[,4] 
prot_input <- weights$proteo[,4]
metab_inputs <- weights$metab[,4]

#remove MOFA tags from feature names
names(RNA_input) <- gsub("_RNA","",names(RNA_input))
names(prot_input) <- gsub("_proteo","",names(prot_input))

RNA_log2_FPKM <- as.data.frame(read_csv("data/RNA/RNA_log2_FPKM_clean.csv"))$Genes #since only complete cases were considered in the first part, here the full gene list is loaded
```

We can be interested in looking how RNA and corresponding proteins
correlate in each factor. We can plot this using scatter plots for each
factors.

``` r
weights_RNA <- as.data.frame(weights$RNA)
weights_prot <- as.data.frame(weights$proteo)

weights_prot$ID <- gsub("_proteo$","",row.names(weights_prot))
weights_RNA$ID <- gsub("_RNA$","",row.names(weights_RNA))

plot_list <- list()
r2_list <- list()
for(i in 1:(dim(weights_prot)[2]-1))
{
  merged_weights <- merge(weights_RNA[,c(i,10)],weights_prot[,c(i,10)], by = "ID")
  names(merged_weights) <- c("ID","RNA","prot")
  r2_list[[i]] <- cor(merged_weights[,2],merged_weights[,3])^2
  plot_list[[i]] <- ggplot(merged_weights,aes(x = RNA,y = prot)) + geom_point() +
  geom_smooth(method='lm', formula= y~x) + 
  theme_minimal() + 
  theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank())
}

r2_vec <- unlist(r2_list)
r2_vec
```

    ## [1] 0.041804090 0.411093174 0.038741311 0.421953345 0.002666596 0.291243355
    ## [7] 0.027585257 0.058177119 0.019215245

``` r
ggsave(filename = "results/mofa/RNA_prot_cor.pdf",plot = do.call("grid.arrange", c(plot_list, ncol = 3)), device = "pdf")
```

    ## Saving 7 x 7 in image

![](MOFA_to_COSMOS_files/figure-gfm/correlation%20RNA/prot-1.png)<!-- -->
Coherently, only the factors where there is variance explained for both
RNA and proteins are correlated.

Next, we perform filtering steps in order to identify top deregulated
features. Here, the individual filtering is based on absolute factor
weight or activity score. The first step involves the visualization of
the respective distribution to define an appropriate threshold.

First, the RNA factor weights are filtered based on a threshold defined
by their distribution as well a threshold defined by the distribution of
the protein factor weights.

``` r
##RNA
{plot(density(RNA_input))
abline(v = -0.2)
abline(v = 0.2)}
```

![](MOFA_to_COSMOS_files/figure-gfm/Plot:%20RNA%20weights%20&%20protein%20weights-1.png)<!-- -->

``` r
{plot(density(prot_input))
abline(v = -0.05)
abline(v = 0.05)}
```

![](MOFA_to_COSMOS_files/figure-gfm/Plot:%20RNA%20weights%20&%20protein%20weights-2.png)<!-- -->

Based on the plots and indicated by the straight lines, the RNA weights
threshold is set to -0.2 and 0.2, and the protein weights threshold to
-0.05 and 0.05. RNA factor weight values lying inside this threshold are
set to 0 and previously filtered out genes (before MOFA analysis) are
also included with value 0.

``` r
for(gene in names(RNA_input)) {
  if (RNA_input[gene] > -0.2 & RNA_input[gene] < 0.2)
  {
    RNA_input[gene] <- 0
  } else
  {
    RNA_input[gene] <- sign(RNA_input[gene]) * 10
  }
  if (gene %in% names(prot_input)) #will need to add a sign consistency check between RNA and prot as well
  {
    if (prot_input[gene] > -0.05 & prot_input[gene] < 0.05)
    {
      RNA_input[gene] <- 0
    } else
    {
      RNA_input[gene] <- sign(RNA_input[gene]) * 10
    }
  }
}

RNA_log2_FPKM <- RNA_log2_FPKM[!(RNA_log2_FPKM %in% names(RNA_input))]
genes <- RNA_log2_FPKM
RNA_log2_FPKM <- rep(0,length(RNA_log2_FPKM))
names(RNA_log2_FPKM) <- genes

RNA_input <- c(RNA_input, RNA_log2_FPKM)
```

The same procedure is repeated for the transcription factor activities.

``` r
## Transcription factors
{plot(density(signaling_input))
abline(v = -0.5)
abline(v = 3.5)}
```

![](MOFA_to_COSMOS_files/figure-gfm/Plot:%20TF%20activity%20estimates-1.png)<!-- -->

The threshold is set to -0.5 and 3.5 respectively. Further, the
TF_to_remove variable is later used to remove transcription factors with
activities below this threshold from the prior knowledge network.

``` r
TF_to_remove <- signaling_input[signaling_input > -0.5 & signaling_input < 3.5]
signaling_input <- signaling_input[signaling_input < -0.5 | signaling_input > 3.5]
```

The same procedure is repeated for the ligand-receptor activities.

``` r
## Ligand-receptor interactions
{plot(density(ligrec_input))
abline(v = -0.5)
abline(v = 2.5)}
```

![](MOFA_to_COSMOS_files/figure-gfm/Plot:%20Ligand-receptor%20activity%20estimates-1.png)<!-- -->

Here, only activities higher than 2.5 or lower than -0.5 are kept (see
straight lines in plot). This is an arbitrary theshold aimed at keeping
a relativelly high number of interesting fetures to connect in the rest
of the pipeline. Since at this point we are interested in the receptors,
the names are adjusted accordingly and if there are multiple mentions of
a receptor, its activity is calculated by the mean of the associated
ligand-receptor interactions.

``` r
ligrec_input <- ligrec_input[ligrec_input > 2.5 | ligrec_input < -0.5]

rec_inputs <- ligrec_input
names(rec_inputs) <- gsub(".+___","",names(rec_inputs))
rec_inputs <- tapply(rec_inputs, names(rec_inputs), mean)
receptors <- names(rec_inputs)
rec_inputs <- as.numeric(rec_inputs)
names(rec_inputs) <- receptors
```

Finally, the same procedure is repeated for the metabolite factor
weights.

``` r
## Metabolites
{plot(density(metab_inputs))
abline(v = -0.2)
abline(v = 0.2)}
```

![](MOFA_to_COSMOS_files/figure-gfm/Plot:%20Metabolite%20factor%20weights-1.png)<!-- -->

The threshold is set to 0.2 and -0.2 respectively. Further, the
metab_to_exclude variable is later used to remove metabolites with
activities below this threshold from the prior knowledge network.
Moreover, metabolite names are translated into HMDB format and metabolic
compartment codes (m = mitochondria, c = cytosol) as well as metab\_\_
prefix is added to HMDB IDs. More information regarding compartment
codes can be found [here](http://bigg.ucsd.edu/compartments).

``` r
metab_to_HMDB <- as.data.frame(
  read_csv("data/metabolomic/MetaboliteToHMDB.csv"))
metab_to_HMDB <- metab_to_HMDB[metab_to_HMDB$common %in% names(metab_inputs),]
metab_inputs <- metab_inputs[metab_to_HMDB$common]
names(metab_inputs) <- metab_to_HMDB$HMDB

metab_inputs <- cosmosR::prepare_metab_inputs(metab_inputs, compartment_codes = c("m","c"))

metab_to_exclude <- metab_inputs[abs(metab_inputs) < 0.2]
metab_inputs <- metab_inputs[abs(metab_inputs) > 0.2]
```

Next, the prior knowledge network (PKN) is loaded. To see how the full
meta PKN was assembled, see
[PKN](https://github.com/saezlab/meta_PKN_BIGG.git). Since we are not
interested in including transcription factors and metabolites that were
previously filtered out, we also remove the respective nodes from the
PKN.

``` r
## Load and filter meta network
# data("meta_network")
# save(meta_network, file = "support/meta_network.RData")
load("support/meta_network.RData")

meta_network <- meta_network_cleanup(meta_network)

meta_network <- meta_network[!(meta_network$source %in% names(TF_to_remove)) 
                             & !(meta_network$target %in% names(TF_to_remove)),]
meta_network <- meta_network[!(meta_network$source %in% names(metab_to_exclude)) 
                             & !(meta_network$target %in% names(metab_to_exclude)),]
```

Then the datasets are merged, entries that are not included in the PKN
are removed and the filtering steps are completed. The merging in this
case is based on the idea of deriving the forward network via the
receptor to transcription factor/metabolite direction. We will perform
first the first cosmos run between recptors and downstream TFs and
metabolites. The second run will focus on connecting TFs awith
downstream ligands.

Thus, we first define the upstream input as receptors, and the
downstream input as TFs and metabolites. Since the scoring of the
network is sensitive to the scale of the downstream measurments, we are
multiplying the metabolite weights by 10 so that they are on a similar
scale as the TF scores. Other scaling methods can be performed at the
user’S discretion, this one is just a rough scaling, but good enough in
this case.

``` r
# upstream_inputs <- c(rec_inputs)
# downstream_inputs <- c(metab_inputs*10, signaling_input)
# 
# upstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(upstream_inputs, meta_network)
# downstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(downstream_inputs, meta_network)
# 
# rec_to_TF_cosmos_input <- list(upstream_inputs_filtered, downstream_inputs_filtered)
# save(rec_to_TF_cosmos_input, file = "data/cosmos/rec_to_TF_cosmos_input.Rdata")
```

### MOON (meta footprint analysis)

In this part, we will run the network scoring to generate mechanistic
hypotheses connecting ligands, receptors, TFs and metabolites, while
pruning the network with information from both RNA and protein abundance
to constrain the flux of signal consistent with the molecular data
available.

In this approach, the network is compressed to avoid having redundant
paths in the network leading to inflating the importance of some input
measurements.

``` r
##filter expressed genes from PKN
# data("meta_network")
# save(meta_network, file = "support/meta_network.RData")
# load("support/meta_network.RData")

# meta_network <- meta_network_cleanup(meta_network)

expressed_genes <- as.data.frame(read_csv("data/RNA/RNA_log2_FPKM_clean.csv"))
```

    ## Rows: 11265 Columns: 61
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (1): Genes
    ## dbl (60): 786-0, A498, A549/ATCC, ACHN, BT-549, CAKI-1, CCRF-CEM, COLO 205, ...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
expressed_genes <- setNames(rep(1,length(expressed_genes[,1])), expressed_genes$Genes)
meta_network_filtered <- cosmosR:::filter_pkn_expressed_genes(names(expressed_genes), meta_pkn = meta_network)
```

    ## [1] "COSMOS: removing unexpressed nodes from PKN..."
    ## [1] "COSMOS: 14430 interactions removed"

``` r
##format metab inputs
metab_inputs <- as.numeric(scale(weights$metab[,4], center = F))
names(metab_inputs) <- row.names(weights$metab)

metab_to_HMDB <- as.data.frame(
  read_csv("data/metabolomic/MetaboliteToHMDB.csv"))
```

    ## Rows: 139 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (2): common, HMDB
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
metab_to_HMDB <- metab_to_HMDB[metab_to_HMDB$common %in% names(metab_inputs),]
metab_inputs <- metab_inputs[metab_to_HMDB$common]
names(metab_inputs) <- metab_to_HMDB$HMDB

metab_inputs <- cosmosR::prepare_metab_inputs(metab_inputs, compartment_codes = c("m","c"))
```

    ## [1] "Adding compartment codes."

``` r
#prepare upstream inputs
upstream_inputs <- c(rec_inputs) #the upstream input should be filtered for most significant

TF_inputs <- scale(ligrec_TF_moon_inputs$TF, center = F)
TF_inputs <- setNames(TF_inputs[,1], row.names(TF_inputs))

downstream_inputs <- c(metab_inputs, TF_inputs) #the downstream input should be complete

upstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(upstream_inputs, meta_network_filtered)
```

    ## [1] "COSMOS: 5 input/measured nodes are not in PKN any more: BCAM, CD63, LRP10, NOTCH1, RXRA and 0 more."

``` r
downstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(downstream_inputs, meta_network_filtered)
```

    ## [1] "COSMOS: 429 input/measured nodes are not in PKN any more: Metab__HMDB0011747_m, Metab__HMDB0000755_m, Metab__HMDB0001294_m, Metab__HMDB0001508_m, Metab__HMDB0000012_m, Metab__HMDB0001191_m and 423 more."

``` r
#Filter inputs and prune the meta_network to only keep nodes that can be found downstream of the inputs
#The number of step is quite flexible, 7 steps already covers most of the network

n_steps <- 6

# in this step we prune the network to keep only the relevant part between upstream and downstream nodes
meta_network_filtered <- cosmosR:::keep_controllable_neighbours(meta_network_filtered
                                                       , n_steps, 
                                                       names(upstream_inputs_filtered))
```

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 6 steps"
    ## [1] "COSMOS: 29459 from  42695 interactions are removed from the PKN"

``` r
downstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(downstream_inputs_filtered, meta_network_filtered)
```

    ## [1] "COSMOS: 28 input/measured nodes are not in PKN any more: Metab__HMDB0000168_m, Metab__HMDB0000056_m, Metab__HMDB0001051_m, Metab__HMDB0000139_m, Metab__HMDB0000687_m, Metab__HMDB0000696_m and 22 more."

``` r
meta_network_filtered <- cosmosR:::keep_observable_neighbours(meta_network_filtered, n_steps, names(downstream_inputs_filtered))
```

    ## [1] "COSMOS: removing nodes that are not observable by measurements within 6 steps"
    ## [1] "COSMOS: 6200 from  13236 interactions are removed from the PKN"

``` r
upstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(upstream_inputs_filtered, meta_network_filtered)
```

    ## [1] "COSMOS: 4 input/measured nodes are not in PKN any more: EPHA6, GPC1, SDC1, SIRPA and 0 more."

``` r
write_csv(meta_network_filtered, file = "results/cosmos/moon/meta_network_filtered.csv")

#compress the network to avoid redundant master controllers
meta_network_compressed_list <- compress_same_children(meta_network_filtered, sig_input = upstream_inputs_filtered, metab_input = downstream_inputs_filtered)

meta_network_compressed <- meta_network_compressed_list$compressed_network

meta_network_compressed <- meta_network_cleanup(meta_network_compressed)

load(file = "support/dorothea_df.RData")

# RNA_input <- as.numeric(weights$RNA[,4])
# names(RNA_input) <- gsub("_RNA$","",row.names(weights$RNA))
```

In this step, we format the MOFA weight and activity score information
so that it can be appended to the moon network result afterward.

``` r
data("HMDB_mapper_vec")

## MOFA weights
MOFA_weights <- get_weights(model, factors = 4, as.data.frame = T)
MOFA_weights <- MOFA_weights %>%
  arrange(feature, view) %>%
  filter(!duplicated(feature))
MOFA_weights <- MOFA_weights[,c(1,3)]
colnames(MOFA_weights) <- c("Nodes", "mofa_weights")

#remove RNA and proteo tags and integrate values
MOFA_weights$Nodes <- gsub("_RNA$","",MOFA_weights$Nodes)
MOFA_weights$Nodes <- gsub("_proteo$","",MOFA_weights$Nodes)
MOFA_weights$sign <- sign(MOFA_weights$mofa_weights)
MOFA_weights$mofa_weights <- abs(MOFA_weights$mofa_weights)

MOFA_weights <- MOFA_weights %>% group_by(Nodes) %>% summarise_each(funs(max(., na.rm = TRUE)))
```

    ## Warning: `summarise_each()` was deprecated in dplyr 0.7.0.
    ## ℹ Please use `across()` instead.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: `funs()` was deprecated in dplyr 0.8.0.
    ## ℹ Please use a list of either functions or lambdas:
    ## 
    ## # Simple named list: list(mean = mean, median = median)
    ## 
    ## # Auto named with `tibble::lst()`: tibble::lst(mean, median)
    ## 
    ## # Using lambdas list(~ mean(., trim = .2), ~ median(., na.rm = TRUE))
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
MOFA_weights <- as.data.frame(MOFA_weights)
MOFA_weights$mofa_weights <- MOFA_weights$mofa_weights * MOFA_weights$sign
MOFA_weights <- MOFA_weights[,-3]

#make metabolite names coherent between mofa and cosmos
#for now i just add the metab input that is already formatted
#I will miss some weights that were not input but will fix that at later time
MOFA_weights_metabaddon <- as.data.frame(metab_inputs)
names(MOFA_weights_metabaddon) <- "mofa_weights"
MOFA_weights_metabaddon$Nodes <- row.names(MOFA_weights_metabaddon)

MOFA_weights_metabaddon[, 2] <- sapply(MOFA_weights_metabaddon[, 2], function(x, HMDB_mapper_vec) {
        x <- gsub("Metab__", "", x)
        suffixe <- stringr::str_extract(x, "_[a-z]$")
        x <- gsub("_[a-z]$", "", x)
        if (x %in% names(HMDB_mapper_vec)) {
            x <- HMDB_mapper_vec[x]
            x <- paste("Metab__", x, sep = "")
        }
        if (!is.na(suffixe)) {
            x <- paste(x, suffixe, sep = "")
        }
        return(x)
    }, HMDB_mapper_vec = HMDB_mapper_vec)

MOFA_weights <- as.data.frame(rbind(MOFA_weights, MOFA_weights_metabaddon))


## Ligands
lig_weights <- ligrec_TF_moon_inputs$ligrec
names(lig_weights) <- gsub("___.+","",names(lig_weights))
lig_weights <- tapply(lig_weights, names(lig_weights), mean)
ligands <- names(lig_weights)
lig_weights <- as.numeric(lig_weights)
names(lig_weights) <- ligands
lig_weights <- data.frame(Nodes = names(lig_weights), feature_weights = lig_weights)

## Receptors
rec_weights <- ligrec_TF_moon_inputs$ligrec
names(rec_weights) <- gsub(".+___","",names(rec_weights))
rec_weights <- tapply(rec_weights, names(rec_weights), mean)
receptors <- names(rec_weights)
rec_weights <- as.numeric(rec_weights)
names(rec_weights) <- receptors
rec_weights <- data.frame(Nodes = names(rec_weights), feature_weights = rec_weights)

## TF
TF_weights <- ligrec_TF_moon_inputs$TF
names(TF_weights) <- gsub("_TF","",names(TF_weights))
TF_weights <- TF_weights[!(names(TF_weights) %in% c(rec_weights$Nodes,lig_weights$Nodes))]
TF_weights <- data.frame(Nodes = names(TF_weights), feature_weights = TF_weights)



## Combine
feature_weights <- as.data.frame(rbind(TF_weights, rec_weights, lig_weights))
```

To also identify whether a node is a ligand or receptor we also load the
liana consensus LR ressource.

``` r
# LIANA
load("support/ligrec_ressource.RData")
ligrec_df <- ligrec_ressource[,c("source_genesymbol","target_genesymbol")]
ligrec_df <- distinct(ligrec_df)
names(ligrec_df) <- c("Node1","Node2")

ligrec_df$Node1 <- gsub("-","_",ligrec_df$Node1)
ligrec_df$Node2 <- gsub("-","_",ligrec_df$Node2)
ligrec_df$Sign <- 1
ligrec_df$Weight <- 1
```

In this step. we run the first moon analysis, to connect receptors with
downstream TFs and ligands. MOON function is run iterativelly in a loop
until the resulting scored network does not include edges incoherent
with the RNA and proteomic data measurents. then, the network is
decompressed, IDs are translated in a human-friendly format and a
subnetowrk comparable to the classic cosmos solution network is
generated with the reduce_solution_network function.

``` r
meta_network_rec_to_TFmetab <- meta_network_compressed

#We run moon in a loop until TF-target coherence convergences
before <- 1
after <- 0
i <- 1
while (before != after & i < 10) {
  before <- length(meta_network_rec_to_TFmetab[,1])
  start_time <- Sys.time()
  moon_res <- cosmosR::moon(upstream_input = upstream_inputs_filtered, 
                                                 downstream_input = downstream_inputs_filtered, 
                                                 meta_network = meta_network_rec_to_TFmetab, 
                                                 n_layers = n_steps, 
                                                 statistic = "ulm") 
  
  meta_network_rec_to_TFmetab <- filter_incohrent_TF_target(moon_res, dorothea_df, meta_network_rec_to_TFmetab, RNA_input)
  end_time <- Sys.time()
  
  print(end_time - start_time)
  
  after <- length(meta_network_rec_to_TFmetab[,1])
  i <- i + 1
}
```

    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5
    ## [1] 6
    ## Time difference of 0.1666861 secs
    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5
    ## [1] 6
    ## Time difference of 0.2143788 secs

``` r
if(i < 10)
{
  print(paste("Converged after ",paste(i-1," iterations", sep = ""),sep = ""))
} else
{
  print(paste("Interupted after ",paste(i," iterations. Convergence uncertain.", sep = ""),sep = ""))
}
```

    ## [1] "Converged after 2 iterations"

``` r
moon_res <- decompress_moon_result(moon_res, meta_network_compressed_list, meta_network_rec_to_TFmetab)

#save the whole res for later
moon_res_rec_to_TFmet <- moon_res[,c(4,2,3)]
moon_res_rec_to_TFmet[,1] <- translate_column_HMDB(moon_res_rec_to_TFmet[,1], HMDB_mapper_vec)

levels <- moon_res[,c(4,3)]

moon_res <- moon_res[,c(4,2)]
names(moon_res)[1] <- "source"

plot(density(moon_res$score))
abline(v = 1)
abline(v = -1)
```

![](MOFA_to_COSMOS_files/figure-gfm/Run%20moon%20rec_to_TFmetab-1.png)<!-- -->

``` r
solution_network <- reduce_solution_network(decoupleRnival_res = moon_res, 
                                            meta_network = meta_network_filtered,
                                            cutoff = 1.5, 
                                            upstream_input = upstream_inputs_filtered, 
                                            RNA_input = RNA_input, 
                                            n_steps = n_steps)
```

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 6 steps"
    ## [1] "COSMOS: 385 from  901 interactions are removed from the PKN"

``` r
SIF_rec_to_TFmetab <- solution_network$SIF
names(SIF_rec_to_TFmetab)[3] <- "sign"
ATT_rec_to_TFmetab <- solution_network$ATT

translated_res <- translate_res(SIF_rec_to_TFmetab,ATT_rec_to_TFmetab,HMDB_mapper_vec)

levels_translated <- translate_res(SIF_rec_to_TFmetab,levels,HMDB_mapper_vec)[[2]]

SIF_rec_to_TFmetab <- translated_res[[1]]
ATT_rec_to_TFmetab <- translated_res[[2]]

## Add weights to data
ATT_rec_to_TFmetab <- merge(ATT_rec_to_TFmetab, MOFA_weights, all.x = T)
ATT_rec_to_TFmetab <- merge(ATT_rec_to_TFmetab, feature_weights, all.x = T)
names(ATT_rec_to_TFmetab)[2] <- "AvgAct"

ATT_rec_to_TFmetab$NodeType <- ifelse(ATT_rec_to_TFmetab$Nodes %in% levels_translated[levels_translated$level == 0,1],1,0)

ATT_rec_to_TFmetab$NodeType <- ifelse(ATT_rec_to_TFmetab$Nodes %in% ligrec_df$Node1, 2, ifelse(ATT_rec_to_TFmetab$Nodes %in% ligrec_df$Node2, 3, ifelse(ATT_rec_to_TFmetab$Nodes %in% dorothea_df$source, 4, ATT_rec_to_TFmetab$NodeType))) 

names(SIF_rec_to_TFmetab)[4] <- "Weight"

write_csv(SIF_rec_to_TFmetab,file = "results/cosmos/moon/SIF_rec_TFmetab.csv")
write_csv(ATT_rec_to_TFmetab,file = "results/cosmos/moon/ATT_rec_TFmetab.csv")
```

Next, we also score a network connecting the TF with downstream ligands.
The procedure is the same as the previous section, only the upstream and
downstream input definition is changing.

``` r
##filter expressed genes from PKN
dorothea_PKN <- dorothea_df[,c(1,3,2)]
names(dorothea_PKN)[2] <- "interaction"

dorothea_PKN_filtered <- cosmosR:::filter_pkn_expressed_genes(names(expressed_genes), meta_pkn = dorothea_PKN)
```

    ## [1] "COSMOS: removing unexpressed nodes from PKN..."
    ## [1] "COSMOS: 12053 interactions removed"

``` r
#prepare upstream inputs

#the upstream input should be filtered for most significant
upstream_inputs <- setNames(TF_weights$feature_weights, TF_weights$Nodes)
upstream_inputs <- upstream_inputs[abs(upstream_inputs) > 2]

lig_inputs <- ligrec_TF_moon_inputs$ligrec
names(lig_inputs) <- gsub("___.+","",names(lig_inputs))
lig_inputs <- tapply(lig_inputs, names(lig_inputs), mean)
ligands <- names(lig_inputs)
lig_inputs <- as.numeric(lig_inputs)
names(lig_inputs) <- ligands

lig_inputs <- scale(lig_inputs, center = F)
lig_inputs <- setNames(lig_inputs[,1], row.names(lig_inputs))

downstream_inputs <- lig_inputs #the downstream input should be complete

upstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(upstream_inputs, dorothea_PKN_filtered)
```

    ## [1] "COSMOS: 22 input/measured nodes are not in PKN any more: SPI1, NR0B2, ESR1, AR, NKX2-5, RARB and 16 more."

``` r
downstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(downstream_inputs, dorothea_PKN_filtered)
```

    ## [1] "COSMOS: 19 input/measured nodes are not in PKN any more: ACTR2, ADAM9, AGRN, ARPC5, EFNA3, EFNB1 and 13 more."

``` r
#Filter inputs and prune the meta_network to only keep nodes that can be found downstream of the inputs
#The number of step is quite flexible, 7 steps already covers most of the network

n_steps <- 1

# in this step we prune the network to keep only the relevant part between upstream and downstream nodes
dorothea_PKN_filtered <- cosmosR:::keep_controllable_neighbours(dorothea_PKN_filtered
                                                       , n_steps, 
                                                       names(upstream_inputs_filtered))
```

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 1 steps"
    ## [1] "COSMOS: 1896 from  17542 interactions are removed from the PKN"

``` r
downstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(downstream_inputs_filtered, dorothea_PKN_filtered)
```

    ## [1] "COSMOS: 9 input/measured nodes are not in PKN any more: BMP1, ETV5, GAS6, GDF11, ITGB3BP, MAML2 and 3 more."

``` r
dorothea_PKN_filtered <- cosmosR:::keep_observable_neighbours(dorothea_PKN_filtered, n_steps, names(downstream_inputs_filtered))
```

    ## [1] "COSMOS: removing nodes that are not observable by measurements within 1 steps"
    ## [1] "COSMOS: 12699 from  15646 interactions are removed from the PKN"

``` r
upstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(upstream_inputs_filtered, dorothea_PKN_filtered)
```

    ## [1] "COSMOS: 7 input/measured nodes are not in PKN any more: MTF1, NR2C2, RBPJ, SRF, NCOA2, E2F2 and 1 more."

``` r
# load(file = "support/dorothea_df.RData")

# RNA_input <- as.numeric(weights$RNA[,4])
# names(RNA_input) <- gsub("_RNA$","",row.names(weights$RNA))
```

``` r
meta_network_TF_lig <- dorothea_PKN_filtered

write_csv(meta_network_TF_lig, file = "results/cosmos/moon/meta_network_TF_lig.csv")

before <- 1
after <- 0
i <- 1
while (before != after & i < 10) {
  before <- length(meta_network_TF_lig[,1])
  moon_res <- cosmosR::moon(upstream_input = upstream_inputs_filtered, 
                                                 downstream_input = downstream_inputs_filtered, 
                                                 meta_network = meta_network_TF_lig, 
                                                 n_layers = n_steps, 
                                                 statistic = "ulm") 
  
  meta_network_TF_lig <- filter_incohrent_TF_target(moon_res, dorothea_df, meta_network_TF_lig, RNA_input)
  after <- length(meta_network_TF_lig[,1])
  i <- i + 1
}

if(i < 10)
{
  print(paste("Converged after ",paste(i-1," iterations", sep = ""),sep = ""))
} else
{
  print(paste("Interupted after ",paste(i," iterations. Convergence uncertain.", sep = ""),sep = ""))
}
```

    ## [1] "Converged after 1 iterations"

``` r
moon_res_TF_lig <- moon_res
moon_res_TF_lig[,1] <- translate_column_HMDB(moon_res_TF_lig[,1], HMDB_mapper_vec)

levels <- moon_res[,c(1,3)]

moon_res <- moon_res[,c(1,2)]
names(moon_res)[1] <- "source"

plot(density(moon_res$score))
abline(v = 1)
abline(v = -1)
```

![](MOFA_to_COSMOS_files/figure-gfm/run%20moon%20TF%20to%20lig-1.png)<!-- -->

``` r
solution_network <- reduce_solution_network(decoupleRnival_res = moon_res, 
                                            meta_network = as.data.frame(dorothea_PKN_filtered[,c(1,3,2)]),
                                            cutoff = 1.5, 
                                            upstream_input = upstream_inputs_filtered, 
                                            RNA_input = RNA_input, 
                                            n_steps = n_steps)
```

    ## [1] "COSMOS: removing nodes that are not reachable from inputs within 1 steps"
    ## [1] "COSMOS: 29 from  143 interactions are removed from the PKN"

``` r
SIF_TF_lig <- solution_network$SIF
names(SIF_TF_lig)[3] <- "sign"
ATT_TF_lig <- solution_network$ATT


translated_res <- translate_res(SIF_TF_lig,ATT_TF_lig,HMDB_mapper_vec)

levels_translated <- translate_res(SIF_TF_lig,levels,HMDB_mapper_vec)[[2]]

SIF_TF_lig <- translated_res[[1]]
ATT_TF_lig <- translated_res[[2]]

## Add weights to data
ATT_TF_lig <- merge(ATT_TF_lig, MOFA_weights, all.x = T)
ATT_TF_lig <- merge(ATT_TF_lig, feature_weights, all.x = T)
names(ATT_TF_lig)[2] <- "AvgAct"

ATT_TF_lig$NodeType <- ifelse(ATT_TF_lig$Nodes %in% levels_translated[levels_translated$level == 0,1],1,0)

ATT_TF_lig$NodeType <- ifelse(ATT_TF_lig$Nodes %in% ligrec_df$Node1, 2, ifelse(ATT_TF_lig$Nodes %in% ligrec_df$Node2, 3, ifelse(ATT_TF_lig$Nodes %in% dorothea_df$source, 4, ATT_TF_lig$NodeType))) 

names(SIF_TF_lig)[4] <- "Weight"

write_csv(SIF_TF_lig,file = "results/cosmos/moon/SIF_TF_lig.csv")
write_csv(ATT_TF_lig,file = "results/cosmos/moon/ATT_TF_lig.csv")
```

Both networks (Receptos -\> TF and metabs; TF -\> Ligands) are merged
together into a single compbined scored network. We also add the
conections between ligands and their corresponding receptors.

``` r
combined_SIF_moon <- as.data.frame(rbind(SIF_rec_to_TFmetab, SIF_TF_lig))
combined_SIF_moon <- unique(combined_SIF_moon)

combined_ATT_moon <- as.data.frame(rbind(ATT_rec_to_TFmetab, ATT_TF_lig))
combined_ATT_moon <- combined_ATT_moon %>% group_by(Nodes) %>% summarise_each(funs(mean(., na.rm = TRUE)))
```

    ## Warning: `summarise_each()` was deprecated in dplyr 0.7.0.
    ## ℹ Please use `across()` instead.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: `funs()` was deprecated in dplyr 0.8.0.
    ## ℹ Please use a list of either functions or lambdas:
    ## 
    ## # Simple named list: list(mean = mean, median = median)
    ## 
    ## # Auto named with `tibble::lst()`: tibble::lst(mean, median)
    ## 
    ## # Using lambdas list(~ mean(., trim = .2), ~ median(., na.rm = TRUE))
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
combined_ATT_moon <- as.data.frame(combined_ATT_moon)

ligrec_ressource_addon <- ligrec_ressource[
  ligrec_ressource$source_genesymbol %in% combined_SIF_moon$target &
    ligrec_ressource$target_genesymbol %in% combined_SIF_moon$source
, c(10,12)]
ligrec_ressource_addon$sign <- 1
ligrec_ressource_addon$Weight <- 1
names(ligrec_ressource_addon)[c(1,2)] <- c("source","target")
ligrec_ressource_addon <- unique(ligrec_ressource_addon)

combined_SIF_moon <- as.data.frame(rbind(combined_SIF_moon, ligrec_ressource_addon))

#It appears I may need to only consider direct TF regulations for TF to lig

write_csv(combined_SIF_moon,file = "results/cosmos/moon/combined_SIF_moon.csv")
write_csv(combined_ATT_moon,file = "results/cosmos/moon/combined_ATT_moon.csv")

write_csv(moon_res_rec_to_TFmet, file = "results/cosmos/moon/moon_res_rec_to_TFmet.csv")
write_csv(moon_res_TF_lig, file = "results/cosmos/moon/moon_res_TF_lig.csv")
```

### Further analyses

Further MOFA-COSMOS downstream analyses focusing on analyzing network
control structures and cell line specific MOFA transcriptomics are given
in ./pathway_control_analysis_moon.Rmd and ./Cell_line_MOFA_space.Rmd

### Session info

``` r
sessionInfo()
```

    ## R version 4.2.0 (2022-04-22)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Monterey 12.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] RColorBrewer_1.1-3   gridExtra_2.3        pheatmap_1.0.12     
    ##  [4] ggfortify_0.4.15     ggplot2_3.4.0        tidyr_1.3.0         
    ##  [7] GSEABase_1.58.0      graph_1.74.0         annotate_1.74.0     
    ## [10] XML_3.99-0.13        AnnotationDbi_1.58.0 IRanges_2.30.1      
    ## [13] S4Vectors_0.34.0     Biobase_2.56.0       BiocGenerics_0.42.0 
    ## [16] reshape2_1.4.4       dplyr_1.1.3          decoupleR_2.9.1     
    ## [19] liana_0.1.5          MOFA2_1.6.0          cosmosR_1.5.2       
    ## [22] readr_2.1.4         
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] utf8_1.2.3                  reticulate_1.28            
    ##   [3] tidyselect_1.2.0            RSQLite_2.2.20             
    ##   [5] htmlwidgets_1.6.1           grid_4.2.0                 
    ##   [7] BiocParallel_1.30.4         Rtsne_0.16                 
    ##   [9] devtools_2.4.5              munsell_0.5.0              
    ##  [11] ScaledMatrix_1.4.1          ragg_1.2.5                 
    ##  [13] codetools_0.2-18            statmod_1.5.0              
    ##  [15] scran_1.24.1                future_1.30.0              
    ##  [17] miniUI_0.1.1.1              withr_2.5.0                
    ##  [19] colorspace_2.1-0            progressr_0.13.0           
    ##  [21] filelock_1.0.2              highr_0.10                 
    ##  [23] logger_0.2.2                knitr_1.42                 
    ##  [25] rstudioapi_0.14             SingleCellExperiment_1.18.1
    ##  [27] listenv_0.9.0               labeling_0.4.2             
    ##  [29] MatrixGenerics_1.8.1        GenomeInfoDbData_1.2.8     
    ##  [31] farver_2.1.1                bit64_4.0.5                
    ##  [33] rhdf5_2.40.0                basilisk_1.8.1             
    ##  [35] parallelly_1.34.0           vctrs_0.6.1                
    ##  [37] generics_0.1.3              xfun_0.42                  
    ##  [39] R6_2.5.1                    doParallel_1.0.17          
    ##  [41] GenomeInfoDb_1.32.4         clue_0.3-63                
    ##  [43] rsvd_1.0.5                  locfit_1.5-9.7             
    ##  [45] bitops_1.0-7                rhdf5filters_1.8.0         
    ##  [47] cachem_1.0.7                DelayedArray_0.22.0        
    ##  [49] vroom_1.6.1                 promises_1.2.0.1           
    ##  [51] scales_1.2.1                gtable_0.3.1               
    ##  [53] beachmat_2.12.0             globals_0.16.2             
    ##  [55] processx_3.8.0              rlang_1.1.0                
    ##  [57] systemfonts_1.0.4           splines_4.2.0              
    ##  [59] GlobalOptions_0.1.2         checkmate_2.1.0            
    ##  [61] BiocManager_1.30.19         yaml_2.3.7                 
    ##  [63] backports_1.4.1             httpuv_1.6.8               
    ##  [65] tools_4.2.0                 usethis_2.1.6              
    ##  [67] ellipsis_0.3.2              sessioninfo_1.2.2          
    ##  [69] Rcpp_1.0.10                 plyr_1.8.8                 
    ##  [71] sparseMatrixStats_1.8.0     progress_1.2.2             
    ##  [73] zlibbioc_1.42.0             purrr_1.0.1                
    ##  [75] RCurl_1.98-1.10             ps_1.7.2                   
    ##  [77] basilisk.utils_1.8.0        prettyunits_1.1.1          
    ##  [79] GetoptLong_1.0.5            cowplot_1.1.1              
    ##  [81] urlchecker_1.0.1            SeuratObject_4.1.3         
    ##  [83] SummarizedExperiment_1.26.1 ggrepel_0.9.2              
    ##  [85] cluster_2.1.4               fs_1.6.1                   
    ##  [87] magrittr_2.0.3              circlize_0.4.15            
    ##  [89] matrixStats_0.63.0          pkgload_1.3.2              
    ##  [91] hms_1.1.3                   mime_0.12                  
    ##  [93] evaluate_0.20               xtable_1.8-4               
    ##  [95] readxl_1.4.2                shape_1.4.6                
    ##  [97] compiler_4.2.0              tibble_3.2.1               
    ##  [99] crayon_1.5.2                htmltools_0.5.5            
    ## [101] mgcv_1.8-41                 later_1.3.0                
    ## [103] tzdb_0.3.0                  DBI_1.1.3                  
    ## [105] corrplot_0.92               ComplexHeatmap_2.15.4      
    ## [107] rappdirs_0.3.3              Matrix_1.5-3               
    ## [109] cli_3.6.1                   parallel_4.2.0             
    ## [111] metapod_1.4.0               igraph_1.4.2               
    ## [113] GenomicRanges_1.48.0        forcats_1.0.0              
    ## [115] pkgconfig_2.0.3             dir.expiry_1.4.0           
    ## [117] OmnipathR_3.9.6             sp_1.6-0                   
    ## [119] scuttle_1.6.3               xml2_1.3.3                 
    ## [121] foreach_1.5.2               dqrng_0.3.0                
    ## [123] XVector_0.36.0              rvest_1.0.3                
    ## [125] stringr_1.5.0               callr_3.7.3                
    ## [127] digest_0.6.31               Biostrings_2.64.1          
    ## [129] rmarkdown_2.21              cellranger_1.1.0           
    ## [131] uwot_0.1.14                 edgeR_3.38.4               
    ## [133] DelayedMatrixStats_1.18.2   curl_5.0.0                 
    ## [135] shiny_1.7.4                 rjson_0.2.21               
    ## [137] nlme_3.1-161                lifecycle_1.0.3            
    ## [139] jsonlite_1.8.4              Rhdf5lib_1.18.2            
    ## [141] BiocNeighbors_1.14.0        limma_3.52.4               
    ## [143] fansi_1.0.4                 pillar_1.9.0               
    ## [145] lattice_0.20-45             KEGGREST_1.36.3            
    ## [147] fastmap_1.1.1               httr_1.4.5                 
    ## [149] pkgbuild_1.4.0              glue_1.6.2                 
    ## [151] remotes_2.4.2               png_0.1-8                  
    ## [153] iterators_1.0.14            bit_4.0.5                  
    ## [155] bluster_1.6.0               stringi_1.7.12             
    ## [157] profvis_0.3.7               HDF5Array_1.24.2           
    ## [159] blob_1.2.3                  textshaping_0.3.6          
    ## [161] BiocSingular_1.12.0         memoise_2.0.1              
    ## [163] irlba_2.3.5.1               future.apply_1.10.0
