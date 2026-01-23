# install and load packages
#install.packages("remotes")

#remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
#remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

setwd("c:/Bioinformatics/00_Daniocell_data")

library(dplyr)
library(tidyverse)

#remotes::install_github("mojaveazure/seurat-disk", force = TRUE)

library(SeuratDisk)

# run the line below only once
#Convert("zf_atlas_full_v4_release.h5ad", dest="h5seurat", overwrite=TRUE)

#STEP2:loading the h5Seurat file
zhub_seurat <- LoadH5Seurat("zf_atlas_full_v4_release.h5seurat",assays = "RNA")

# store meta data in a separate object
zebrahub_meta <- zhub_seurat@meta.data

nrow(zebrahub_meta)

# collect the data
assay_data <- zhub_seurat@assays

