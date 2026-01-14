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

###############################################################################
# Part 1: Obtain the tissue distribution for a specific gene and visualize it across stages

############################## Daniocell dataset ###############################



################################## isolation of gene-specific data #############

# collect the data
assay_data <- daniocell@assays

# target gene
plex9_1 <- "si:ch1073-385f13.3"
plex9_2 <- "si:ch1073-296i8.2"

gene <- plex9_1

# subset the data to the input gene
gene_data <- as.data.frame(assay_data$RNA@data[gene,])

# complete the data frame
colnames(gene_data) <- c("expression")
gene_data$cell_id <- rownames(gene_data)
gene_data$cluster <- daniocell@meta.data[rownames(gene_data),]$cluster
gene_data$tissue <- daniocell@meta.data[rownames(gene_data),]$tissue
gene_data$stage <- daniocell@meta.data[rownames(gene_data),]$stage.group


############################# large slice aggregation for size factors #######
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

# make 

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

write.csv(gene_data_long_summary, "./size_factors_data/daniocell_stage_gene_counts.csv")


gene_wide_summary <- pivot_wider(gene_data_long_summary,     
                                 names_from = stage,
                                 values_from = counts_sum)

gene_wide_summary <- as.data.frame(gene_wide_summary)


rownames(gene_wide_summary) <- gene_wide_summary$gene

gene_wide_summary <- gene_wide_summary[,2:15]

write.csv(gene_wide_summary, "./size_factors_data/daniocell_count_sums.csv", quote = FALSE)

##############################################################################
# complete the data frame
colnames(gene_data) <- c("expression")
gene_data$cell_id <- rownames(gene_data)
gene_data$cluster <- daniocell@meta.data[rownames(gene_data),]$cluster
gene_data$tissue <- daniocell@meta.data[rownames(gene_data),]$tissue
gene_data$stage <- daniocell@meta.data[rownames(gene_data),]$stage.group

############## convert the log-normalized values to count-like
gene_data$counts <- exp(gene_data$expression) - 1

# subset by expression
nrow(gene_data)
# [1] 489686

gene_data <- gene_data[gene_data$counts >= 1,] # 1 corresponds to 0 and 
# > 2 would correspond to some real reads

nrow(gene_data)
# [1] 77909

# explore how many cells do not have tissue information
sum(gene_data$tissue == "", na.rm = TRUE)

# 5123 out of 77909 selected cells have empty tissue value

# clean the tissue factor 
gene_data <- gene_data[gene_data$tissue != '',]
gene_data <- gene_data[!is.na(gene_data$tissue),]

nrow(gene_data)
# [1] 72786

library(dplyr)
options(scipen=999)


# tissue level aggregation
gene_expression_summary <- gene_data |> 
  group_by(tissue, stage) |> 
  summarise(counts_sum = sum(counts)) |> 
  arrange(tissue, stage, .by_group = TRUE)

################## normalize the total counts by the size factors ############ 


write.csv(gene_expression_summary, "akna_expression_summary.csv")

############################### Visualization of the data ####################

library(ggplot2)

gene <- "plex9.1"

ggplot(gene_expression_summary, aes(x = stage, y = counts_sum, fill = stage)) +
  geom_col() +
  facet_wrap(~tissue, scales = "free_x") +
  ylab("Cell expression value sum")+
  ggtitle(paste("Summarised expression plot for", gene , "in zebrafish"))+
  theme(strip.text.x = element_text(size = 8, margin = margin(0.12,0,0.12,0, "cm")),
        axis.text.x = element_text(angle = 70, vjust=0.7, colour="grey20", size= 7.5, face="plain"),
        legend.title = element_text( size = 10, face = "bold"),
        legend.text = element_text( size = 10, face = "plain"),
        legend.key.size = unit(0.5, "cm"),
        panel.spacing = unit(0.1, "lines")
  )


##################### averaged expression

# tissue level aggregation
gene_mean_expr_summary <- gene_data |> 
  group_by(tissue, stage) |> 
  summarise(expr_mean = mean(expression)) |> 
  arrange(tissue, stage, .by_group = TRUE)

#write.csv(gene_expression_summary, "akna_expression_summary.csv")

############################### Visualization of the data ####################

library(ggplot2)

ggplot(gene_mean_expr_summary, aes(x = stage, y = expr_mean, fill = stage)) +
  geom_col() +
  facet_wrap(~tissue, scales = "free_x") +
  theme(strip.text.x = element_text(size = 8),
        axis.text.x = element_text(angle = 70, vjust=0.7, colour="grey20", size= 8, face="plain"),
        legend.title = element_text( size = 12, face = "bold"),
        legend.text = element_text( size = 12, face = "plain"),
        legend.key.size = unit(0.5, "cm"),
        panel.spacing = unit(0.5, "lines")
  )

########################## generating stage and tissue summary ################

gene_names <- rownames(assay_data$RNA@data)

write.csv(gene_names, "gene_names.csv")


########## Algorithm

# Pre-compute:
# 1. Subset the large dataset by gene name and filter.
# 2. Compute summaries.
# 3. Store list of summaries for future retrieval.


# Retrieval:
# 1. Retrieve the stored list
# 2. Subset to the input gene.
# 3. Make plots.



assay_data <- daniocell@assays


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


# Consider other potentially problematic symbols.  

# UPDATE: these missing genes are in reality duplicates

# Process and save data frame with tissue and stage count sums
for(gene in gene_names){

  
  # check if the gene is already processed
  print(gene)
  
  if(gene %in% colon_dash_genes){
    
    next
    
  }
    

  # 1.1 Subset the data 
  
  # subset the data to the input gene
  gene_data <- as.data.frame(assay_data$RNA@data[gene,])

  
  # the gene name has to work in the process of saving the data
  if(gene %in% slash_genes){
    
    gene <- str_replace(gene, '/', '')
    
  }
  
  #also handle the genes with ':' and '-'
  
  # if(gene %in% colon_dash_genes){
  #   
  #   gene <- str_replace(gene, ':', '')
  #   gene <- str_replace(gene, '-', '')
  #   
  # }


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

  print(gene_wide_summary)
  
  write.csv(gene_wide_summary, paste0("./new/daniocell_", gene, "_count_sums.csv"), quote = FALSE, row.names = FALSE)

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
  
  write.csv(gene_wide_summary, paste0("./single_genes_cell_counts/daniocell_", gene, "_cell_counts.csv"), quote = FALSE, row.names = FALSE)
  
}


################################## Figure out which genes were not saved ######
# 1. Process gene names as was done when saving the data frames.

slash_genes <- gene_names[str_detect(gene_names, '/')]

colon_dash_genes <- gene_names[str_detect(gene_names, ':')]

processed_gene_names <- c()

for(gene in gene_names){

  
  if(gene %in% slash_genes){
    
    gene <- str_replace(gene, '/', '')
    
  }
  
  
  if(gene %in% colon_dash_genes){
    
    gene <- str_replace(gene, ':', '')
    gene <- str_replace(gene, '-', '')
    
  }
  
  processed_gene_names <- c(processed_gene_names, gene)
  
  
}

# 2. Obtain all file names and extract the processed file names.

all_file_names <- list.files("./single_genes_dfs")

names_in_files <- c()

for(fn in all_file_names){
 items <- as.vector(str_split(fn, '_'))

 names_in_files <- c(names_in_files, items[[1]][2])

}

# 3. Determine the set difference between the original processed gene names and
# the gene names extracted from the file names.

missing_genes <- setdiff(processed_gene_names, names_in_files)
missing_genes

# [1] "ABCD2"    "ADAM12"   "ARF5"     "ARHGAP22" "CACNA2D3" "CALHM2"   "CEP170B"  "CMKLR1"   "COLGALT1" "FAM163B"  "GAA"     
# [12] "GABRR1"   "NAT16"    "NPC2"     "NUP85"    "PAMR1"    "PCMTD2"   "PLPP7"    "RFESD"    "RGS13"    "RNF14"    "RNF180"  
# [23] "SLC29A4"  "SLC45A4"  "SLC46A3"  "SLC7A1"   "SLCO3A1"  "TMEM235"  "TSTA3"    "ZNF521"   "march4"   "ust" 





