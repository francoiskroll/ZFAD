# Exploring single-cell RNAseq data from Raj et al., 2020 (Schier lab)

### Download/prepare the data

Start with **prepData.R**. It will download the RNAseq data (Seurat object) for each developmental stage (12 of them) and create one .rds object called _seuratAll.rds_.  

All the files in folder _seuratdata_ are collated into _seuralAll.rds_. Having one file is easier to handle/store. So once you ran **prepData.R** and you have _seuratAll.rds_, it is perhaps a good idea to delete folder _seuratdata_ all together as it is heavy (I put in Online only on Dropbox).  

### Look for zebrafish orthologues of Alzheimer's-related genes

This is happening in **ZFADexpression.R**.
