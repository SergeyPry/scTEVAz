# install and load packages
#install.packages("remotes")

#remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
#remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

setwd("c:/Bioinformatics/00_Daniocell_data")

library(dplyr)
library(tidyverse)
library(stringr)


# UPDATE THESE PATHS TO REFLECT WHERE YOU HAVE PUT THE FILES:
path.to.seurat.object <- "Daniocell2023_SeuratV4.rds"
path.to.annotations <- "cluster_annotations.csv"

# Load the Daniocell Seurat object and the latest annotations
daniocell <- readRDS(path.to.seurat.object)
daniocell.annot <- read.csv(path.to.annotations)

# To add the latest annotations to the Seurat metadata
rownames(daniocell.annot) <- daniocell.annot$clust
daniocell.annot <- daniocell.annot[,c("tissue", "identity.super", "identity.sub", "identity.super.short", "identity.sub.short", "zfin")]
daniocell@meta.data <- cbind(daniocell@meta.data, daniocell.annot[daniocell@meta.data$cluster,])