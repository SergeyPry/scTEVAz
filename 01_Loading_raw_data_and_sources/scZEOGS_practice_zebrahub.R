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

write.csv(gene_data_long_summary, "./size_factors_data/zebrahub_stage_gene_counts.csv")


gene_wide_summary <- pivot_wider(gene_data_long_summary,     
                                 names_from = stage,
                                 values_from = counts_sum)

gene_wide_summary <- as.data.frame(gene_wide_summary)


rownames(gene_wide_summary) <- gene_wide_summary$gene

gene_wide_summary <- gene_wide_summary[,2:11]

write.csv(gene_wide_summary, "./size_factors_data/zebrahub_count_sums.csv", quote = FALSE)




################################## isolation of gene-specific data #############


# collect the data
assay_data <- zhub_seurat@assays


# target gene
# target gene

plex9_1 <- "si:ch1073-385f13.3"
plex9_2 <- "si:ch1073-296i8.2"

gene <- plex9_2

gene_names <- as.vector(assay_data$RNA@data@Dimnames[[1]])
gene_data <- as.data.frame(assay_data$RNA@data[gene,])

# complete the data frame
colnames(gene_data) <- c("expression")

sum(gene_data$expression > 0)
# [1] 8935

gene_data$stage <- zebrahub_meta$developmental_stage
gene_data$timepoint <- zebrahub_meta$timepoint
gene_data$cell_lineage<- zebrahub_meta$zebrafish_anatomy_ontology_class

############## convert the log-normalized values to count-like values
gene_data$counts <- exp(gene_data$expression)




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


size_factors_df <- data.frame(stage = c("0 somites", "05 somites", "10 somites", "15 somites", "20 somites", "30 somites", "larval-2dpf", "larval-3dpf", "larval-5dpf", "larval-10dpf"),
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
for(gene in colon_dash_genes){
  
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
  
  gene_data$stage <- zebrahub_meta$developmental_stage
  gene_data$timepoint <- zebrahub_meta$timepoint
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
                                   arrange(stage)
  
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
stage_cell_nums <- table(zebrahub_meta$developmental_stage)


stage_cell_factors <- as.vector(stage_cell_nums/12914)

stage_cell_factors_df <- data.frame(stage = c("0 somites", "05 somites", "10 somites", "15 somites", "20 somites", "30 somites", "larval-2dpf", "larval-3dpf", "larval-5dpf", "larval-10dpf"),
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
  
  gene_data$stage <- zebrahub_meta$developmental_stage
  gene_data$timepoint <- zebrahub_meta$timepoint
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
                                   arrange(stage)
  
  gene_wide_summary[is.na(gene_wide_summary)] <- 0
  
  write.csv(gene_wide_summary, paste0("./new/zebrahub_", gene, "_cell_counts.csv"), quote = FALSE, row.names = FALSE)
  
}



