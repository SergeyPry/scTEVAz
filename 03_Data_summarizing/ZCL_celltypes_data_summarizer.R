# install and load packages
#install.packages("remotes")

#remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
#remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

#setwd("c:/Bioinformatics/00_Daniocell_data")

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

nrow(pbmc@meta.data)
#[1] 1088106

nrow(ZCL_metadata)
#[1] 1082680

################### Store data and gene names ###############################
assay_data <- pbmc@assays
gene_names <- as.vector(assay_data$RNA@data@Dimnames[[1]])


########## Algorithm

# Pre-compute:
# 1. Subset the large dataset by gene name and filter.
# 2. Compute summaries.
# 3. Store list of summaries for future retrieval.


# Retrieval:
# 1. Retrieve the stored list
# 2. Subset to the input gene.
# 3. Make plots.

# prevent scientific notation in R
options(scipen=999)

size_factors_df <- data.frame(stage = c("21Day" ,  "22Month", "24hpf",  "3Month", "72hpf"),
                              size_factors = c(0.7470865, 1.8293085, 0.2447139, 3.1254952, 0.9347619))

size_factors_df$stage <- factor(size_factors_df$stage, levels = c("24hpf", "72hpf", "21Day" , "3Month",  "22Month"   ))

slash_genes <- gene_names[str_detect(gene_names, '/')]
colon_dash_genes <- gene_names[str_detect(gene_names, ':')]

# the idea is that if the user requests such a gene
# you convert the id by removing : and - and then retrieve its data.
# to obtain the proper name, it will make sense to have a dictionary to map the derived name
# to the original one.

# dictionary/list generation:
# if a gene contains a '/' or ':', 
#     keep the original and make a derived version without special symbols
# add the derived version as key and the original one as the value and
# store the resulting list as a data frame and retrieve for later use 
# when handling user inputs


# Process and save data frame with tissue and stage count sums
for(gene in gene_names){
  
  
  # keep track of the 
  print(gene)
  
  
  # 1.1 Subset the data 
  
  # subset the data to the input gene
  gene_data <- as.data.frame(assay_data$RNA@data[gene,])
  
  
  # the gene name has to work in the process of saving the data
  if(gene %in% slash_genes){
    
    gene <- str_replace(gene, '/', '')
    
  }
  
  #also handle the genes with ':' and '-'
  
  if(gene %in% colon_dash_genes){
    
    gene <- str_replace(gene, ':', '')
    gene <- str_replace(gene, '-', '')
    
  }
  
  
  colnames(gene_data) <- c("counts")
  
  gene_data$cell_id <- rownames(gene_data)
  
  
  # filter out the empty cells
  gene_data <- gene_data[gene_data$counts >= 1,]  
  
  gene_data$stage <- ZCL_metadata[rownames(gene_data),]$stage
  gene_data$cell_type <- ZCL_metadata[rownames(gene_data),]$cell_type
  
  
  # remove the NA data
  gene_data <- gene_data[!is.na(gene_data$cell_type),]
  
  
  ############################### SUm of counts expression distribution #########
  ###############################################################################
  
  gene_data <- inner_join(gene_data, size_factors_df, by = "stage")
  
  # adjust the counts by size factors
  gene_data$counts_norm <- gene_data$counts/gene_data$size_factors
  
  # order the stages
  gene_data$stage <- factor(gene_data$stage, levels = c("24hpf", "72hpf","21Day", "3Month","22Month"))
  
  
  gene_expression_summary <- gene_data |> 
    group_by(cell_type, stage) |> 
    summarise(counts_sum = sum(counts_norm)) |> 
    arrange(cell_type, stage, .by_group = TRUE)
  
  
  gene_wide_summary <- pivot_wider(gene_expression_summary,     
                                   names_from = cell_type,
                                   values_from = counts_sum) |> 
    arrange(stage)
  
  gene_wide_summary[is.na(gene_wide_summary)] <- 0
  
  write.csv(gene_wide_summary, paste0("./ZCL_dataset_celltypes/single_genes_dfs/ZCL_", gene, "_count_sums.csv"), quote = FALSE, row.names = FALSE)
  
}

##################################### normalized cell numbers output ###########

############################ Cell counts distribution ########################

# generate cell number factors to account for differences in cell abundance 
# between stages

ZCL_metadata$stage <- factor(ZCL_metadata$stage, levels = c("24hpf", "72hpf","21Day", "3Month","22Month"))
stage_cell_nums <- table(ZCL_metadata$stage)
stage_cell_nums

# 24hpf   72hpf   21Day  3Month 22Month 
# 44771  158848  121007  430420  327634 


# normalize by the 24 hpf
stage_cell_factors <- as.vector(stage_cell_nums/44771)

stage_cell_factors_df <- data.frame(stage = c("24hpf", "72hpf", "21Day", "3Month", "22Month"),
                                    size_factors = stage_cell_factors)


for(gene in gene_names[20347:length(gene_names)]){
  
  
  # check if the gene is already processed
  print(gene)
  
  # 1.1 Subset the data 
  
  # subset the data to the input gene
  gene_data <- as.data.frame(assay_data$RNA@data[gene,])
  
  
  # the gene name has to work in the process of saving the data
  if(gene %in% slash_genes){
    
    gene <- str_replace(gene, '/', '')
    
  }
  
  #also handle the genes with ':' and '-'
  
  if(gene %in% colon_dash_genes){
    
    gene <- str_replace(gene, ':', '')
    gene <- str_replace(gene, '-', '')
    
  }
  
  colnames(gene_data) <- c("counts")
  
  gene_data$cell_id <- rownames(gene_data)
  
  # filter out the empty cells
  gene_data <- gene_data[gene_data$counts >= 1,] 
  
  gene_data$stage <- ZCL_metadata[rownames(gene_data),]$stage
  gene_data$cell_type <- ZCL_metadata[rownames(gene_data),]$cell_type
  
  
  # remove the NA data
  gene_data <- gene_data[!is.na(gene_data$cell_type),]
  
  ##################### cell counts per tissue and stage
  
  # cell-type level aggregation
  gene_counts <- gene_data |> 
    group_by(cell_type, stage) |> 
    summarise(cell_count = n()) |> 
    arrange(cell_type, stage, .by_group = TRUE)
  
  
  gene_counts <- inner_join(gene_counts , stage_cell_factors_df, by = "stage")
  
  # order the stages
  gene_counts$stage <- factor(gene_counts$stage, levels = c("24hpf", "72hpf","21Day", "3Month","22Month"))
  
  
  # adjust the counts by size factors
  gene_counts$counts_norm <- gene_counts$cell_count/gene_counts$size_factors
  
  
  
  gene_counts_summary <- gene_counts[,c('cell_type', 'stage', 'counts_norm')]
  
  gene_wide_summary <- pivot_wider(gene_counts_summary,     
                                   names_from = cell_type,
                                   values_from = counts_norm) |> 
                                   arrange(stage)
  
  gene_wide_summary[is.na(gene_wide_summary)] <- 0
  
  write.csv(gene_wide_summary, paste0("./ZCL_dataset_celltypes/single_genes_cell_counts/ZCL_", gene, "_cell_counts.csv"), quote = FALSE, row.names = FALSE)
  
}

