# single-cell Tissue Expression Values Aggregator in zebrafish - scTEVAz

## Loading raw data and sources

The following sections comment and describe the data sources and code needed to generate processed files for this project. We provide the information on the original sources of the data and describe the initial steps of preparing these datasets for usage in the later sections of the project.

For all the datasets we ran the following code to install packages:

```r
# install and load packages
install.packages("remotes")
remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
install.packages("paletteer")
remotes::install_github("mojaveazure/seurat-disk", force = TRUE)
```

### Daniocell

To obtain the complete Daniocell data, we visited this project website:

https://daniocell.nichd.nih.gov/

and downloaded the following files

* Daniocell2023_SeuratV4.rds (Seurat object)
* cluster_annotations.csv (Cluster annotations)
* daniocell_load.R (R script to load object and annotations)

We modified the original script slightly and here is the final code we used:

```r
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(paletteer)


# UPDATE THESE PATHS TO REFLECT WHERE YOU HAVE PUT THE FILES:
path.to.seurat.object <- "Daniocell2023_SeuratV4.rds"
path.to.annotations <- "cluster_annotations.csv"

# Load the Daniocell Seurat object and the latest annotatio]ns
daniocell <- readRDS(path.to.seurat.object)
daniocell.annot <- read.csv(path.to.annotations)

# To add the latest annotations to the Seurat metadata
rownames(daniocell.annot) <- daniocell.annot$clust
daniocell.annot <- daniocell.annot[,c("tissue", "identity.super", "identity.sub", "identity.super.short", "identity.sub.short", "zfin")]
daniocell@meta.data <- cbind(daniocell@meta.data, daniocell.annot[daniocell@meta.data$cluster,])
```
The code is also available here: [Daniocell_load_data.R](https://github.com/SergeyPry/scTEVAz/tree/main/01_Loading_raw_data_and_sources/Daniocell_load_data.R).
Later sections will discuss how we can use the imported Seurat object to extract the data for individual genes or groups of genes.


### Zebrafish Cell Landscape

To obtain the Zebrafish Cell Landscape (ZCL) data, we first visited its About page: 
https://bis.zju.edu.cn/ZCL/index.html 
and then clicked on its data repository page:
https://figshare.com/s/1ab3c6d7648d12247eb2

and downloaded the following files from the repository:
* ZCDL.rdata (Seurat object)
* ZCDL_cellinfo.csv (Cell annotations)


We used the following code to import the data:
```r
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(paletteer)

# UPDATE THESE PATHS TO REFLECT WHERE YOU HAVE PUT THE FILES:
path.to.seurat.object <- "ZCDL.rdata"
path.to.annotations <- "ZCDL_cellinfo.csv"

# Load the Daniocell Seurat object and the latest annotations
load(path.to.seurat.object)

ZCL_metadata <- read.csv(path.to.annotations)

rownames(ZCL_metadata) <- ZCL_metadata$barcodes

ZCL_metadata <- ZCL_metadata[,c("cluster", "stage", "cell_type", "cell_lineage", "tsne_x", "tsne_y" )]
```
The code is also available here: [ZCL_load_data.R](https://github.com/SergeyPry/scTEVAz/tree/main/01_Loading_raw_data_and_sources/ZCL_load_data.R).
Later sections will discuss how we can use the imported Seurat object to extract the data for individual genes or groups of genes.

### Zebrahub
To obtain the Zebrahub data, we first visited its Data page: 
https://zebrahub.sf.czbiohub.org/data

and then clicked on its data Google Drive page:
https://drive.google.com/file/d/1jwoNqdeRTsauf5WwRCBvRQ_zdtzFBON7/view?usp=drive_link

We downloaded the following file from the repository: **zf_atlas_full_v4_release.h5ad**

Here is the code to import the data:
```r
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(SeuratDisk)

# run the line below only once
#Convert("zf_atlas_full_v4_release.h5ad", dest="h5seurat", overwrite=TRUE)

#STEP2:loading the h5Seurat file
zhub_seurat <- LoadH5Seurat("zf_atlas_full_v4_release.h5seurat",assays = "RNA")

# store meta data in a separate object
zebrahub_meta <- zhub_seurat@meta.data
```
The code is also available here: [Zebrahub_load_data.R](https://github.com/SergeyPry/scTEVAz/tree/main/01_Loading_raw_data_and_sources/Zebrahub_load_data.R).
Later sections will discuss how we can use the imported Seurat object to extract the data for individual genes or groups of genes.
