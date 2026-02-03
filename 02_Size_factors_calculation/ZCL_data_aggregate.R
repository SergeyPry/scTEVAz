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


############################## ZCL dataset ###############################

############################# large slice aggregation for size factors #######

# collect data

assay_data <- pbmc@assays
gene_names <- as.vector(assay_data$RNA@data@Dimnames[[1]])


# subset the data to the current slice
gene_data <- as.data.frame(as.matrix(assay_data$RNA@data[gene_names[1:100],]))

# transpose a data frame
gene_data <- as.data.frame(t(gene_data))

# add metadata
gene_data$cell_id <- rownames(gene_data)
gene_data$stage <- pbmc@meta.data[rownames(gene_data),]$stage


# make a longer version of the data frame

gene_data <- gene_data |> 
  pivot_longer(cols = -c(cell_id, stage),
               names_to = "gene", 
               values_to = "counts")


# aggregate genes and stages

gene_data_long_summary <- gene_data |> 
  group_by( gene, stage) |> 
  summarise(counts_sum = sum(counts)) |> 
  arrange(gene, stage, .by_group = TRUE)


################################## loop over all genes ######################

for(i in seq(101, 27501, by=100) ){
  print(i)
  
  # handle the end of the sequence
  if(i == 27501){
    j = 27538
  } else{
    j = i + 99
  }
  
  
  # subset the data to the current slice
  gene_data <- as.data.frame(as.matrix(assay_data$RNA@data[gene_names[i: j],]))
  
  # transpose a data frame
  gene_data <- as.data.frame(t(gene_data))
  
  # add metadata
  gene_data$cell_id <- rownames(gene_data)
  gene_data$stage <- pbmc@meta.data[rownames(gene_data),]$stage
  
  # make a tidy version of the data
  gene_data <- gene_data |> 
    pivot_longer(cols = -c(cell_id, stage),
                 names_to = "gene", 
                 values_to = "counts")
  
  # aggregate genes and stages
  
  gene_data_long_sum <- gene_data|> 
    group_by( gene, stage) |> 
    summarise(counts_sum = sum(counts)) |> 
    arrange(gene, stage, .by_group = TRUE)
  
  
  gene_data_long_summary <- rbind(gene_data_long_summary, gene_data_long_sum)
  
}

# reshape to generate count sums for all genes and stages
gene_wide_summary <- pivot_wider(gene_data_long_summary,     
                                 names_from = stage,
                                 values_from = counts_sum)

gene_wide_summary <- as.data.frame(gene_wide_summary)


rownames(gene_wide_summary) <- gene_wide_summary$gene

gene_wide_summary <- gene_wide_summary[,2:6]

write.csv(gene_wide_summary, "./size_factors_data/ZCL_count_sums.csv", quote = FALSE)



