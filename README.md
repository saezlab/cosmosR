# cosmosR <img src="inst/figures/logo.png" align="right" height="139">

<!-- badges: start -->
<!-- badges: end -->

## Overview

COSMOS (Causal Oriented Search of Multi-Omic Space) is a method that integrates phosphoproteomics, transcriptomics, and metabolomics data sets. COSMOS leverages extensive prior knowledge of signaling pathways, metabolic networks, and gene regulation  with computational methods to estimate activities of transcription factors and kinases as well as network-level causal reasoning. This pipeline can provide mechanistic explanations for experimental observations across multiple omic data sets. 


<img src="man/figures/intro_data.png" align="center" width="800">

COSMOS uses [CARNIVAL](https://saezlab.github.io/CARNIVAL/)’s Integer Linear Programming (ILP) optimization strategy to find the smallest coherent subnetwork causally connecting as many deregulated TFs, kinases/phosphatases and metabolites as possible. The subnetwork is extracted from a novel integrated PKN (available [here](http://metapkn.omnipathdb.org/)) spanning signaling, transcriptional regulation and metabolism.  Transcription factors activities are inferred from gene expression with [DoRothEA](https://saezlab.github.io/dorothea/), a meta resource of TF/target links. Kinase activities are inferred from phosphoproteomic with a kinase/substrate network of [Omnipath](http://omnipathdb.org/), a meta resource of protein-protein. [CARNIVAL](https://saezlab.github.io/CARNIVAL/) was adapted to find mechanistic hypotheses connecting the TF and kinase activities with metabolites from a signaling/metabolic prior knowledge network combining [Omnipath](http://omnipathdb.org/), [STITCHdb](http://stitch.embl.de/) and [Recon3D](https://www.vmh.life/). 


You can also use COSMOS if you don't have metabolomic data, to connect TF activities (from transcriptomic) with kinase activities (from phosphoproteomic) for exmaple !

<img src="man/figures/graphical_abstract.png" align="center" width="800">


## Installation

R >= 4.1 is required
```r
# install from bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("cosmosR")

# We advise to instal from github to get the latest version of the tool.
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("saezlab/cosmosR")
```

If you don't have R 4.1, you can also clone the github repository on your machine, create a new R project with R studio from the cosmosR folder, change the R version to your own R version in the DESCRIPTION file and then install it with devtools:install()

But 4.1 is advised in any case.

### Prerequisites

COSMOS is dependent on CARNIVAL for exhibiting the signalling pathway optimisation.
CARNIVAL requires the interactive version of IBM Cplex, Gurobi or CBC-COIN solver as the network optimiser. The IBM ILOG Cplex is freely available through Academic Initiative [here](https://www.ibm.com/products/ilog-cplex-optimization-studio). Gurobi license is also free for academics, request a license [here](https://www.gurobi.com/downloads/end-user-license-agreement-academic/). The [CBC](https://projects.coin-or.org/Cbc) solver is open source and freely available for any user, but has a significantly lower performance than CPLEX or Gurobi. Obtain CBC executable directly usable for cosmos [here](https://ampl.com/products/solvers/open-source/#cbc). Alternatively for small networks, users can rely on the freely available [lpSolve R-package](https://cran.r-project.org/web/packages/lpSolve/index.html), which is automatically installed with the package.


Small note to package dependencies:

 CARNIVAL is currently under active development. We try to ensure the compatibility of both packages. In case the current CARNIVAL installation non the less colides with COSMOS please report and consider installing an older stable version to run COSMSOS.
```r
# Install an older stable version of CARNIVAL (1.3.0)
remotes::install_github("saezlab/CARNIVAL@b3a84c6ba9706547caca02644566d75ee621f568")
```

## Tutorial (video)

We recorded a video guide for cosmosR tutorial in the context of a course organised by EBI-EMBL. 
You can access the recording at this link for a step by step introduction to cosmosR : 
https://embl-ebi.cloud.panopto.eu/Panopto/Pages/Viewer.aspx?id=318f7091-b6bf-44ee-939f-adb10121fc1b

## Tutorial (NCI60 playground)

We made a repository that contains pre-processed inputs and an example script to use cosmos with the NCI60 RNA+metabolomic datasets.
You can find the repository [here](https://github.com/saezlab/NCI60_cosmos).

## Access

The meta PKN used with the biorXiv version of COSMOS is available [here](http://metapkn.omnipathdb.org/).

An updated meta PKN is available with the package (using data(meta_network) in R)

## Citation
If you use cosmosR for your research please cite the [original publication](https://www.embopress.org/doi/full/10.15252/msb.20209730): 

> Dugourd A, Kuppe C, Sciacovelli M, Gjerga E, Gabor A, Emdal KB, Vieira V, Bekker-Jensen DB, Kranz J, Bindels EMJ, Jesper V Olsen, Christian Frezza, Rafael Kramann, Julio Saez-Rodriguez et al (2021) Causal integration of multi-omics data with prior knowledge to generate mechanistic hypotheses. Mol Syst Biol 17: e9730

## License

The code is distributed under the GNU General Public License v3.0. The meta PKN is distributed under the Attribution-NonCommercial 4.0 International (CC-BY-NC 4.0) License.
