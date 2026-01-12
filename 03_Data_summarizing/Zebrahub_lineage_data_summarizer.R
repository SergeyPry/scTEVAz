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

########################## generating stage and tissue summary ################

gene_names <- rownames(assay_data$RNA@data)

write.csv(gene_names, "zhub_gene_names.csv")


########## Algorithm

# Pre-compute:
# 1. Subset the large dataset by gene name and filter.
# 2. Compute summaries.
# 3. Store list of summaries for future retrieval.


# Retrieval:
# 1. Retrieve the stored list
# 2. Subset to the input gene.
# 3. Make plots.



# collect the data
assay_data <- zhub_seurat@assays


size_factors_df <- data.frame(stage = c("10hpf", "12hpf", "14hpf", "16hpf", "19hpf", "24hpf", "2dpf", "3dpf", "5dpf", "10dpf"),
                              size_factors = c(0.1208711, 0.5704448, 0.4491545, 0.6201209, 0.6681976, 1.2310201, 2.1067936, 3.5070977, 3.5066438, 2.7098503))

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


# Consider other potentially problematic symbols.  

# UPDATE: these missing genes are in reality duplicates

# Process and save data frame with tissue and stage count sums
for(gene in gene_names){
  
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
  
  colnames(gene_data) <- c("expression")
  
  gene_data$stage <- zebrahub_meta$timepoint
  gene_data$cell_lineage<- zebrahub_meta$zebrafish_anatomy_ontology_class
  
  
  ############## convert the log-normalized values to count-like
  gene_data$counts <- exp(gene_data$expression) - 1
  gene_data <- gene_data[gene_data$counts >=1,] 
  
  
  # clean the cell_lineage factor 
  gene_data <- gene_data[gene_data$cell_lineage != '',]
  gene_data <- gene_data[!is.na(gene_data$cell_lineage),]
  
  ############################### Sum of counts expression distribution #########
  ###############################################################################
  
  # adjust the count values by size factors
  gene_data <- inner_join(gene_data, size_factors_df, by = "stage")
  
  # adjust the counts by size factors
  gene_data$counts_norm <- gene_data$counts/gene_data$size_factors
  
  # tissue level aggregation
  gene_expression_summary <- gene_data |> 
    group_by(cell_lineage, stage) |> 
    summarise(counts_sum = round(sum(counts_norm), 2) ) |> 
    arrange(cell_lineage, stage, .by_group = TRUE)
  
  gene_wide_summary <- pivot_wider(gene_expression_summary,     
                                   names_from = cell_lineage,
                                   values_from = counts_sum) |> 
                                   mutate(stage = factor(stage, levels = c("10hpf", "12hpf", "14hpf", "16hpf", "19hpf", "24hpf", "2dpf", "3dpf", "5dpf", "10dpf")))
  
  gene_wide_summary$stage <- factor(gene_wide_summary$stage,levels = c("10hpf", "12hpf", "14hpf", "16hpf", "19hpf", "24hpf", "2dpf", "3dpf", "5dpf", "10dpf"))
  
  gene_wide_summary <- gene_wide_summary |> arrange(stage)
  
  gene_wide_summary[is.na(gene_wide_summary)] <- 0
  
  if(file.exists(paste0("./new/zebrahub_", gene, "_count_sums.csv"))){
    next
  }else{
    
    write.csv(gene_wide_summary, paste0("./new/zebrahub_", gene, "_count_sums.csv"), quote = FALSE, row.names = FALSE)  
  }
  
  
}



##################################### normalized cell numbers output ###########

############################ Cell counts distribution ########################

# generate cell number factors to account for differences in cell abundance 
# between stages
stage_cell_nums <- table(zebrahub_meta$timepoint)


stage_cell_factors <- as.vector(stage_cell_nums/12914)

stage_cell_factors_df <- data.frame(stage = c("10hpf", "12hpf", "14hpf", "16hpf", "19hpf", "24hpf", "2dpf", "3dpf", "5dpf", "10dpf"),
                                    size_factors = stage_cell_factors)


for(gene in gene_names){
  
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
  
  colnames(gene_data) <- c("expression")
  
  gene_data$stage <- zebrahub_meta$timepoint
  gene_data$cell_lineage<- zebrahub_meta$zebrafish_anatomy_ontology_class
  
  
  ############## convert the log-normalized values to count-like
  gene_data$counts <- exp(gene_data$expression) - 1
  gene_data <- gene_data[gene_data$counts >=1,] 
  
  
  # clean the cell_lineage factor 
  gene_data <- gene_data[gene_data$cell_lineage != '',]
  gene_data <- gene_data[!is.na(gene_data$cell_lineage),]
  
  ##################### cell counts per tissue and stage
  
  # tissue level aggregation
  gene_counts <- gene_data |> 
    group_by(cell_lineage, stage) |> 
    summarise(cell_count = n()) |> 
    arrange(cell_lineage, stage, .by_group = TRUE)
  
  
  gene_counts <- inner_join(gene_counts , stage_cell_factors_df, by = "stage")
  
  # adjust the counts by size factors
  gene_counts$counts_norm <- round(gene_counts$cell_count/gene_counts$size_factors,2)
  
  
  gene_counts_summary <- gene_counts[,c('cell_lineage', 'stage', 'counts_norm')]
  
  gene_wide_summary <- pivot_wider(gene_counts_summary,     
                                   names_from = cell_lineage,
                                   values_from = counts_norm) |> 
                                   mutate(stage = factor(stage, levels = c("10hpf", "12hpf", "14hpf", "16hpf", "19hpf", "24hpf", "2dpf", "3dpf", "5dpf", "10dpf")))
  
  gene_wide_summary$stage <- factor(gene_wide_summary$stage,levels = c("10hpf", "12hpf", "14hpf", "16hpf", "19hpf", "24hpf", "2dpf", "3dpf", "5dpf", "10dpf"))
  gene_wide_summary <- gene_wide_summary |> arrange(stage)
  
  gene_wide_summary[is.na(gene_wide_summary)] <- 0
  
  write.csv(gene_wide_summary, paste0("./new/zebrahub_", gene, "_cell_counts.csv"), quote = FALSE, row.names = FALSE)
  
}


