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

# collect the data
assay_data <- daniocell@assays

########################## generating stage and tissue summary ################

gene_names <- rownames(assay_data$RNA@data)

########## Algorithm

# Pre-compute:
# 1. Subset the large dataset by gene name and filter.
# 2. Compute summaries.
# 3. Store list of summaries for future retrieval.


# Retrieval:
# 1. Retrieve the stored list
# 2. Subset to the input gene.
# 3. Make plots.


size_factors_df <- data.frame(stage = c("  3-4","  5-6","  7-9"," 10-12", " 14-21", " 24-34", " 36-46", " 48-58", " 60-70", " 72-82", " 84-94", " 96-106", "108-118", "120"),
                              size_factors = c(0.1072270, 0.1653038, 0.3126715, 0.3333023, 1.2474709, 1.6916442, 1.8080851, 2.7105623, 1.7793399, 2.0617826, 1.8196684, 1.7966809, 2.4282912, 1.0000000))

slash_genes <- gene_names[str_detect(gene_names, '/')]

colon_dash_genes <- gene_names[str_detect(gene_names, ':')]

# the idea is that if the user requests such a gene
# you convert the id by removing : and - and then retrieve its data.
# to obtain the proper name, it will make sense to have a dictionary to map the derived name
# to the original one.

# dictionary/list generation:
# if a gene contains a '/' or ':', 
#     keep the original and make a derived version without special symbols
# add the derived version as key and the original one as the value
# store the resulting list as a data frame and retrieve for later use 
# when handling user inputs

setwd("c:/Bioinformatics/00_Daniocell_data/scTEVAz/Daniocell_dataset")

# Process and save data frame with tissue and stage count sums
for(gene in gene_names){
  
  # 1.1 Subset the data 
  
  # default value
  gene_proc = ''
  
  # the gene name has to work in the process of saving the data
  if(gene %in% slash_genes){
    
    gene_proc <- str_replace(gene, '/', '')
    
  }
  
  #also handle the genes with ':' and '-'
  
  if(gene %in% colon_dash_genes){
    
    gene_proc <- str_replace(gene, ':', '')
    gene_proc <- str_replace(gene_proc, '-', '')
    
  }
  
  # handle the default case
  if(gene_proc == ''){
    gene_proc = gene
  }
  
  
  if(!file.exists(paste0("./single_genes_dfs/daniocell_", gene_proc, "_count_sums.csv"))){
    
    # check if the gene is already processed
    print(gene)
    
    # subset the data to the input gene
    gene_data <- as.data.frame(assay_data$RNA@data[gene,])
    
    # complete the data frame
    colnames(gene_data) <- c("expression")
    gene_data$cell_id <- rownames(gene_data)
    
    ############## convert the log-normalized values to count-like
    gene_data$counts <- exp(gene_data$expression) - 1
    gene_data <- gene_data[gene_data$counts >= 1,] 
    
    gene_data$tissue <- daniocell@meta.data[rownames(gene_data),]$tissue
    gene_data$stage <- daniocell@meta.data[rownames(gene_data),]$stage.group
    
    
    # explore how many cells do not have tissue information
    sum(gene_data$tissue == "", na.rm = TRUE)
    
    # 5123 out of 77909 selected cells have empty tissue value
    
    # clean the tissue factor 
    gene_data <- gene_data[gene_data$tissue != '',]
    gene_data <- gene_data[!is.na(gene_data$tissue),]
    
    ############################### Sum of counts expression distribution #########
    ###############################################################################
    
    # adjust the count values by size factors
    gene_data <- inner_join(gene_data, size_factors_df, by = "stage")
    
    # adjust the counts by size factors
    gene_data$counts_norm <- gene_data$counts/gene_data$size_factors
    
    # tissue level aggregation
    gene_expression_summary <- gene_data |> 
      group_by(tissue, stage) |> 
      summarise(counts_sum = round(sum(counts_norm), 2) ) |> 
      arrange(tissue, stage, .by_group = TRUE)
    
    gene_wide_summary <- pivot_wider(gene_expression_summary,     
                                     names_from = tissue,
                                     values_from = counts_sum) |> 
      arrange(stage)
    
    gene_wide_summary[is.na(gene_wide_summary)] <- 0
    
    
    write.csv(gene_wide_summary, paste0("./new/daniocell_", gene_proc, "_count_sums.csv"), quote = FALSE, row.names = FALSE)
  }
  
}

##################################### normalized cell numbers output ###########

############################ Cell counts distribution ########################

# generate cell number factors to account for differences in cell abundance 
# between stages
stage_cell_nums <- table(daniocell@meta.data$stage.group)
stage_cell_factors <- as.vector(stage_cell_nums/41326)

stage_cell_factors_df <- data.frame(stage = c("  3-4","  5-6","  7-9"," 10-12", " 14-21", " 24-34", " 36-46", " 48-58", " 60-70", " 72-82", " 84-94", " 96-106", "108-118", "120"),
                                    size_factors = stage_cell_factors)


for(gene in gene_names){
  

  # default value
  gene_proc = ''
  
  # the gene name has to work in the process of saving the data
  if(gene %in% slash_genes){
    
    gene_proc <- str_replace(gene, '/', '')
    
  }
  
  #also handle the genes with ':' and '-'
  
  if(gene %in% colon_dash_genes){
    
    gene_proc <- str_replace(gene, ':', '')
    gene_proc <- str_replace(gene_proc, '-', '')
    
  }
  
  # handle the default case
  if(gene_proc == ''){
    gene_proc = gene
  }


  if(!file.exists(paste0("./single_genes_cell_counts/daniocell_", gene_proc, "_cell_counts.csv"))){
    
    # check if the gene is already processed
    print(gene)
    
  # 1.1 Subset the data 
  
  # subset the data to the input gene
  gene_data <- as.data.frame(assay_data$RNA@data[gene,])
  
  # complete the data frame
  colnames(gene_data) <- c("expression")
  
  ############## convert the log-normalized values to count-like
  gene_data$counts <- exp(gene_data$expression) - 1
  gene_data <- gene_data[gene_data$counts >=1,] 
  
  
  gene_data$cell_id <- rownames(gene_data)
  gene_data$tissue <- daniocell@meta.data[rownames(gene_data),]$tissue
  gene_data$stage <- daniocell@meta.data[rownames(gene_data),]$stage.group
  
  
  # explore how many cells do not have tissue information
  sum(gene_data$tissue == "", na.rm = TRUE)
  
  # 5123 out of 77909 selected cells have empty tissue value
  
  # clean the tissue factor 
  gene_data <- gene_data[gene_data$tissue != '',]
  gene_data <- gene_data[!is.na(gene_data$tissue),]
  
  ##################### cell counts per tissue and stage
  
  # tissue level aggregation
  gene_counts <- gene_data |> 
    group_by(tissue, stage) |> 
    summarise(cell_count = n()) |> 
    arrange(tissue, stage, .by_group = TRUE)
  
  
  gene_counts <- inner_join(gene_counts , stage_cell_factors_df, by = "stage")
  
  # adjust the counts by size factors
  gene_counts$counts_norm <- round(gene_counts$cell_count/gene_counts$size_factors,2)
  
  
  gene_counts_summary <- gene_counts[,c('tissue', 'stage', 'counts_norm')]
  
  gene_wide_summary <- pivot_wider(gene_counts_summary,     
                                   names_from = tissue,
                                   values_from = counts_norm) |> 
    arrange(stage)
  
  gene_wide_summary[is.na(gene_wide_summary)] <- 0
  
  write.csv(gene_wide_summary, paste0("./new/daniocell_", gene_proc, "_cell_counts.csv"), quote = FALSE, row.names = FALSE)
  }
}





