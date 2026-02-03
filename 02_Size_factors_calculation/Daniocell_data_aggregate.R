# install and load packages
#install.packages("remotes")

#remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
#remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

#setwd("c:/Bioinformatics/00_Daniocell_data")

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


############################## Daniocell dataset ###############################


############################# large slice aggregation for size factors #######

# collect the data
assay_data <- daniocell@assays
gene_names <- as.vector(assay_data$RNA@data@Dimnames[[1]])

# subset the data to the current slice
gene_data <- as.data.frame(assay_data$RNA@data[gene_names[1:1000],])

# exponentiate the data
gene_data[] <- lapply(gene_data, exp)

# transpose a data frame
gene_data <- as.data.frame(t(gene_data))

# add metadata
gene_data$cell_id <- rownames(gene_data)
gene_data$stage <- daniocell@meta.data[rownames(gene_data),]$stage.group

# subtract 1 from everything
minus_one <- function(x) {
  x - 1
}

gene_data[1:1000] <- lapply(gene_data[1:1000], minus_one) 

# make a tidy data frame

gene_data_long <- gene_data |> 
  pivot_longer(cols = -c(cell_id, stage),
               names_to = "gene", 
               values_to = "counts")


# aggregate genes and stages

gene_data_long_summary <- gene_data_long |> 
  group_by( gene, stage) |> 
  summarise(counts_sum = sum(counts)) |> 
  arrange(gene, stage, .by_group = TRUE)


################################## loop over all genes ######################

for(i in seq(1001, 36001, by=1000) ){
  print(i)
  
  # handle the 
  if(i == 36001){
    j = 36250
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
  gene_data$cell_id <- rownames(gene_data)
  gene_data$stage <- daniocell@meta.data[rownames(gene_data),]$stage.group
      
  # make a tidy version of the data
  gene_data_long <- gene_data |> 
        pivot_longer(cols = -c(cell_id, stage),
                     names_to = "gene", 
                     values_to = "counts")
      
  # aggregate genes and stages
      
  gene_data_long_sum <- gene_data_long |> 
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

gene_wide_summary <- gene_wide_summary[,2:15]

write.csv(gene_wide_summary, "./size_factors_data/daniocell_count_sums.csv", quote = FALSE)

