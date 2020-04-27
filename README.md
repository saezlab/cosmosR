# cosmos <img src="man/figures/logo.png" align="right" height="139">

<!-- badges: start -->
<!-- badges: end -->

## Overview

COSMOS (Causal Oriented Search of Multi-Omic Space) is a method that integrates phosphoproteomics, transcriptomics, and metabolomics data sets. COSMOS leverages extensive prior knowledge of signaling pathways, metabolic networks, and gene regulation  with computational methods to estimate activities of transcription factors and kinases as well as network-level causal reasoning. This pipeline can provide mechanistic explanations for experimental observations across multiple omic data sets. 


<img src="man/figures/graphical_abstract.png" align="center" width="800">

COSMOS uses [CARNIVAL](https://saezlab.github.io/CARNIVAL/)’s Integer Linear Programming (ILP) optimization strategy to find the smallest coherent subnetwork causally connecting as many deregulated TFs, kinases/phosphatases and metabolites as possible. The subnetwork is extracted from a novel integrated PKN (available [here](http://metapkn.omnipathdb.org/)) spanning signaling, transcriptional regulation and metabolism.  Transcription factors activities are inferred from gene expression with [DoRothEA](https://saezlab.github.io/dorothea/), a meta resource of TF/target links. Kinase activities are inferred from phosphoproteomic with a kinase/substrate network of [Omnipath](http://omnipathdb.org/), a meta resource of protein-protein. [CARNIVAL](https://saezlab.github.io/CARNIVAL/) was adapted to find mechanistic hypotheses connecting the TF and kinase activities with metabolites from a signaling/metabolic prior knowledge network combining [Omnipath](http://omnipathdb.org/), [STITCHdb](http://stitch.embl.de/) and [Recon3D](https://www.vmh.life/). 

## Tutorial

Check the Tutorial.pdf (or click [here](https://github.com/saezlab/COSMOS/blob/master/Tutorial.pdf)) file in the root github folder to have an example of running COSMOS with from TF, kinase, phosphatase activities and metabolite abundances. You need to download the github cosmos folder to run the tutorial.

## Access

The integrated PKN used with COSMOS is available [here](http://metapkn.omnipathdb.org/).

## Citation
If you use cosmos for your research please cite the [original publication](https://www.biorxiv.org/content/10.1101/2020.04.23.057893v1): 

> Aurelien dugourd, Christoph Kuppe, Marco Sciacovelli, Enio Gjerga, Kristina Bennet Emdal, Dorte Breinholdt Bekker-Jensen, Jennifer Kranz, Eric J.M. Bindels, Sofia Costa, Jesper V Olsen, Christian Frezza, Rafael Kramann, Julio Saez-Rodriguez. Causal integration of multi-omics data with prior knowledge to generate mechanistic hypotheses. bioRxiv 2020.04.23.057893 (2020) doi:10.1101/2020.04.23.057893.

## License

The code is distributed under the GNU General Public License v3.0. The meta PKN is distributed under the Attribution-NonCommercial 4.0 International (CC-BY-NC 4.0) License.
