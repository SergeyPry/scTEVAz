# install and load packages
#install.packages("remotes")

#remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
#remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

#setwd("c:/Bioinformatics/00_Daniocell_data")

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

############################## Zebrahub dataset ###############################

############################# large slice aggregation for size factors #######
gene_names <- as.vector(assay_data$RNA@data@Dimnames[[1]])

# subset the data to the current slice
gene_data <- as.data.frame(assay_data$RNA@data[gene_names[1:1000],])

# exponentiate the data
gene_data[] <- lapply(gene_data, exp)

# transpose a data frame
gene_data <- as.data.frame(t(gene_data))

# add metadata
gene_data$stage <- zebrahub_meta$developmental_stage
gene_data$timepoint <- zebrahub_meta$timepoint
gene_data$cell_lineage<- zebrahub_meta$zebrafish_anatomy_ontology_class

# subtract 1 from everything
minus_one <- function(x) {
  x - 1
}

gene_data[1:1000] <- lapply(gene_data[1:1000], minus_one) 

# make summary of expression values

gene_data_long <- gene_data |> 
  pivot_longer(cols = -c( stage, timepoint, cell_lineage),
               names_to = "gene", 
               values_to = "counts")


# aggregate genes and stages

gene_data_long_summary <- gene_data_long |> 
  group_by( gene, stage) |> 
  summarise(counts_sum = mean(sum(counts))) |> 
  arrange(gene, stage, .by_group = TRUE)


################################## loop over all genes ######################

for(i in seq(1001, 31001, by=1000) ){
  print(i)
  
  # handle the end of the vector
  if(i == 31001){
    j = 32060
  } else{
    j = i + 999
  }
  
  
  # subset the data to the current slice
  gene_data <- as.data.frame(assay_data$RNA@data[gene_names[i : j],])
  
  # exponentiate the data
  gene_data[] <- lapply(gene_data, exp)
  
  gene_data[] <- lapply(gene_data, minus_one)
  
  # transpose a data frame
  gene_data <- as.data.frame(t(gene_data))
  
  
  # add metadata
  gene_data$stage <- zebrahub_meta$developmental_stage
  gene_data$timepoint <- zebrahub_meta$timepoint
  gene_data$cell_lineage<- zebrahub_meta$zebrafish_anatomy_ontology_class
  
  
  # make summary of expression values
  
  gene_data_long <- gene_data |> 
    pivot_longer(cols = -c( stage, timepoint, cell_lineage),
                 names_to = "gene", 
                 values_to = "counts")
  
  
  # aggregate genes and stages
  
  gene_data_long_sum <- gene_data_long |> 
    group_by( gene, stage) |> 
    summarise(counts_sum = mean(sum(counts))) |> 
    arrange(gene, stage, .by_group = TRUE)
  
  
  gene_data_long_summary <- rbind(gene_data_long_summary, gene_data_long_sum)
  
}


gene_wide_summary <- pivot_wider(gene_data_long_summary,     
                                 names_from = stage,
                                 values_from = counts_sum)

gene_wide_summary <- as.data.frame(gene_wide_summary)


rownames(gene_wide_summary) <- gene_wide_summary$gene

gene_wide_summary <- gene_wide_summary[,2:11]

write.csv(gene_wide_summary, "./size_factors_data/zebrahub_count_sums.csv", quote = FALSE)

