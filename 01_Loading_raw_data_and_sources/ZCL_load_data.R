# install and load packages
#install.packages("remotes")

#remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
#remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

setwd("c:/Bioinformatics/00_Daniocell_data")

library(tidyverse)
library(dplyr)


# UPDATE THESE PATHS TO REFLECT WHERE YOU HAVE PUT THE FILES:
path.to.seurat.object <- "ZCDL.rdata"
path.to.annotations <- "ZCDL_cellinfo.csv"

# Load the Daniocell Seurat object and the latest annotations
load(path.to.seurat.object)

ZCL_metadata <- read.csv(path.to.annotations)

rownames(ZCL_metadata) <- ZCL_metadata$barcodes

ZCL_metadata <- ZCL_metadata[,c("cluster", "stage", "cell_type", "cell_lineage", "tsne_x", "tsne_y" )]