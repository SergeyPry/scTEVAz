# install and load packages
#install.packages("remotes")

#remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
#remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

#install.packages("paletteer")

#setwd("c:/Bioinformatics/00_Daniocell_data")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(paletteer)
library(SeuratObject)
library(Seurat)


# run the line below only once
#Convert("zf_atlas_full_v4_release.h5ad", dest="h5seurat", overwrite=TRUE)

#STEP2:loading the h5Seurat file
zhub_seurat <- LoadH5Seurat("zf_atlas_full_v4_release.h5seurat",assays = "RNA")

# store meta data in a separate object
zebrahub_meta <- zhub_seurat@meta.data

nrow(zebrahub_meta)
#[1] 120444

write.csv(zebrahub_meta, "zebrahub_metadata.csv")

###############################################################################
############################## Zebrahub dataset ############################## 

################################## isolation of gene-specific data ############

# collect the data
assay_data <- zhub_seurat@assays


# target gene
gene <- "akna"

# subset the data to one gene
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


# filter out the empty cells
gene_data <- gene_data[gene_data$counts > 1,]

nrow(gene_data)
# [1] 8935

# prevent scientific notation in R
options(scipen=999)


sum(is.na(gene_data$cell_lineage))
# [1] 0

sum(is.na(gene_data$stage))
# [1] 0


############################### Sum of counts expression distribution #########
###############################################################################


# adjust the count values by size factors

# copy size factors from the size_factors script for this dataset
size_factors_df <- data.frame(stage = c("0 somites", "05 somites", "10 somites", "15 somites", "20 somites", "30 somites", "larval-2dpf", "larval-3dpf", "larval-5dpf", "larval-10dpf"),
                              size_factors = c(0.1208711, 0.5704448, 0.4491545, 0.6201209, 0.6681976, 1.2310201, 2.1067936, 3.5070977, 3.5066438, 2.7098503))
# 
#
size_factors_df$size_factors <- size_factors_df$size_factors/1.2310201

size_factors_df$size_factors

gene_data <- inner_join(gene_data, size_factors_df, by = "stage")

# adjust the counts by size factors
gene_data$counts_norm <- gene_data$counts/gene_data$size_factors

# order the stages
gene_data$stage <- factor(gene_data$stage, levels = c("0 somites", "05 somites", "10 somites", "15 somites", "20 somites", "30 somites", "larval-2dpf", "larval-3dpf", "larval-5dpf", "larval-10dpf"))

# tissue level aggregation
gene_expression_summary <- gene_data |> 
  group_by(cell_lineage, stage) |> 
  summarise(counts_sum = sum(counts_norm)) |> 
  arrange(cell_lineage, stage, .by_group = TRUE)


# sorting by the sum of expression values
gene_expression_summary <- gene_expression_summary %>%
  group_by(tissue) %>%
  mutate(Total_sum = sum(counts_sum)) %>%
  arrange(desc(Total_sum))

# adjusting tissue levels
gene_expression_summary$tissue <- factor(gene_expression_summary$tissue, 
                                         levels = unique(gene_expression_summary.$tissue))


############################### Visualization of the data ####################


# Uncomment if another gene label is required
#gene <- "plex9.2"

# sqrt-scaled plot
ggplot(gene_expression_summary, aes(x = stage, y = counts_sum, fill = stage)) +
  geom_col() +
  scale_y_sqrt(breaks = c(0, 100,250, 500, 1000, 1500, 2000)) +
  facet_wrap(~cell_lineage, scales = "free_x") +
  ylab("Cell expression value sum")+
  ggtitle(paste("Summarised expression plot for", gene , "in zebrafish"))+
  theme_bw() + paletteer::scale_fill_paletteer_d("colorBlindness::paletteMartin") +
  theme(plot.title = element_text(size = 14),
        strip.text.x = element_text(size = 14, margin = margin(0.12,0,0.12,0, "cm")),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 11, face="plain"),
        legend.title = element_text( size = 10, face = "bold"),
        legend.text = element_text( size = 10, face = "plain"),
        legend.key.size = unit(0.5, "cm"),
        panel.spacing = unit(0.1, "lines")
  )

ggsave(paste0(gene, "_counts-sum-normalized_plot_sqrt_zhub.png"), height = 12, width = 15)

# linear plot
ggplot(gene_expression_summary, aes(x = stage, y = counts_sum, fill = stage)) +
  geom_col() +
  scale_y_continuous(breaks = c(0, 250, 500, 1000, 1500, 2000)) +
  facet_wrap(~cell_lineage, scales = "free_x") +
  ylab("Cell expression value sum")+
  ggtitle(paste("Summarised expression plot for", gene , "in zebrafish"))+
  theme_bw() + paletteer::scale_fill_paletteer_d("colorBlindness::paletteMartin") +
  theme(plot.title = element_text(size = 14),
        strip.text.x = element_text(size = 14, margin = margin(0.12,0,0.12,0, "cm")),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 11, face="plain"),
        legend.title = element_text( size = 10, face = "bold"),
        legend.text = element_text( size = 10, face = "plain"),
        legend.key.size = unit(0.5, "cm"),
        panel.spacing = unit(0.1, "lines")
  )

ggsave(paste0(gene, "_counts-sum-normalized_plot_linear_zhub.png"), height = 12, width = 15)

##############################################################################



############################ Cell counts distribution ########################

# generate cell number factors to account for differences in cell abundance 
# between stages

stage_cell_nums <- table(zebrahub_meta$developmental_stage)
stage_cell_nums

# 0 somites   05 somites   10 somites   15 somites   20 somites   30 somites  larval-2dpf 
# 1251         6699         3862         6297         5899        12914        15483 
# larval-3dpf  larval-5dpf larval-10dpf 
# 22473        24987        20579 

# normalize by the 24 hpf
stage_cell_factors <- as.vector(stage_cell_nums/12914)

stage_cell_factors_df <- data.frame(stage = c("0 somites", "05 somites", "10 somites", "15 somites", "20 somites", "30 somites", "larval-2dpf", "larval-3dpf", "larval-5dpf", "larval-10dpf"),
                                    size_factors = stage_cell_factors)

##################### cell counts per tissue and stage

# cell-type level aggregation
gene_counts <- gene_data |> 
  group_by(cell_lineage, stage) |> 
  summarise(cell_count = n()) |> 
  arrange(cell_lineage, stage, .by_group = TRUE)


gene_counts <- inner_join(gene_counts , stage_cell_factors_df, by = "stage")

# order the stages
gene_counts$stage <- factor(gene_counts$stage, levels = c("0 somites", "05 somites", "10 somites", "15 somites", "20 somites", "30 somites", "larval-2dpf", "larval-3dpf", "larval-5dpf", "larval-10dpf"))


# adjust the counts by size factors
gene_counts$counts_norm <- gene_counts$cell_count/gene_counts$size_factors

# sorting by the sum of expression values
gene_counts <- gene_counts %>%
  group_by(tissue) %>%
  mutate(Total_sum = sum(counts_norm)) %>%
  arrange(desc(Total_sum))


############################### Visualization of the data ####################

#linear plot
ggplot(gene_counts, aes(x = stage, y = counts_norm, fill = stage)) +
  geom_col() +
  scale_y_continuous(breaks = c(0,100, 200, 500, 1000, 1250))+
  facet_wrap(~cell_lineage, scales = "free_x") +
  ylab("Cell count")+
  ggtitle(paste("Summarised cell counts for", gene , "in zebrafish"))+
  theme_bw() + paletteer::scale_fill_paletteer_d("colorBlindness::paletteMartin") +
  theme(plot.title = element_text(size = 14),
        strip.text.x = element_text(size = 14, margin = margin(0.12,0,0.12,0, "cm")),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(angle = 70, vjust=0.6, colour="grey20", size= 11, face="plain"),
        legend.title = element_text( size = 10, face = "bold"),
        legend.text = element_text( size = 10, face = "plain"),
        legend.key.size = unit(0.5, "cm"),
        panel.spacing = unit(0.1, "lines")
  )

ggsave(paste0(gene, "_cell_counts-normalized_plot_zhub.png"), height = 12, width = 15)
