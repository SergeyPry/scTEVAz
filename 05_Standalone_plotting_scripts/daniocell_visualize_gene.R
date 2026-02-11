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


# UPDATE THESE PATHS TO REFLECT WHERE YOU HAVE PUT THE FILES:
path.to.seurat.object <- "Daniocell2023_SeuratV4.rds"
path.to.annotations <- "cluster_annotations.csv"

# Load the Daniocell Seurat object and the latest annotatio]ns
daniocell <- readRDS(path.to.seurat.object)
daniocell.annot <- read.csv(path.to.annotations)

# To add the latest annotations to the Seurat metadata
rownames(daniocell.annot) <- daniocell.annot$clust
daniocell.annot <- daniocell.annot[,c("tissue", "identity.super", "identity.sub", "identity.super.short", "identity.sub.short", "zfin")]
daniocell@meta.data <- cbind(daniocell@meta.data, daniocell.annot[daniocell@meta.data$cluster,])

###############################################################################
############################## Daniocell dataset ############################## 


################################## isolation of gene-specific data ############

# collect the data
assay_data <- daniocell@assays

gene <- "pdgfaa"

gene_names <- as.vector(assay_data$RNA@data@Dimnames[[1]])

# subset the data to the current gene
gene_data <- as.data.frame(assay_data$RNA@data[gene,])

# complete the data frame
colnames(gene_data) <- c("expression")
gene_data$cell_id <- rownames(gene_data)
gene_data$cluster <- daniocell@meta.data[rownames(gene_data),]$cluster
gene_data$tissue <- daniocell@meta.data[rownames(gene_data),]$tissue
gene_data$stage <- daniocell@meta.data[rownames(gene_data),]$stage.group


############## convert the log-normalized values to count-like values
gene_data$counts <- exp(gene_data$expression) - 1

# subset by expression
nrow(gene_data)
# [1] 489686

gene_data <- gene_data[gene_data$counts >= 1,] # 1 corresponds to 0 and 
# >= 2 would correspond to some real reads

nrow(gene_data)

# explore how many cells do not have tissue information
sum(gene_data$tissue == "", na.rm = TRUE)


# clean the tissue factor 
gene_data <- gene_data[gene_data$tissue != '',]
gene_data <- gene_data[!is.na(gene_data$tissue),]

nrow(gene_data)
# [1] 8343

# suppress the scientific notation
options(scipen=999)



############################### Sum of counts expression distribution #########
###############################################################################


# adjust the count values by size factors

size_factors_df <- data.frame(stage = c("  3-4","  5-6","  7-9"," 10-12", " 14-21", " 24-34", " 36-46", " 48-58", " 60-70", " 72-82", " 84-94", " 96-106", "108-118", "120"),
                              size_factors = c(0.1072270, 0.1653038, 0.3126715, 0.3333023, 1.2474709, 1.6916442, 1.8080851, 2.7105623, 1.7793399, 2.0617826, 1.8196684, 1.7966809, 2.4282912, 1.0000000))


gene_data <- inner_join(gene_data, size_factors_df, by = "stage")

# adjust the counts by size factors
gene_data$counts_norm <- gene_data$counts/gene_data$size_factors


# tissue level aggregation
gene_expression_summary <- gene_data |> 
  group_by(tissue, stage) |> 
  summarise(counts_sum = round(sum(counts_norm), 2) ) |> 
  arrange(tissue, stage, .by_group = TRUE)

# sorting by the sum of expression values
gene_expression_summary <- gene_expression_summary %>%
  group_by(tissue) %>%
  mutate(Total_sum = sum(counts_sum)) %>%
  arrange(desc(Total_sum))

# adjusting tissue levels
gene_expression_summary$tissue <- factor(gene_expression_summary$tissue, 
                              levels = unique(gene_expression_summary.$tissue))



############################### Visualization of the data ####################


# redefine gene name if the official symbol is too weird
# gene <- "kmt2ca"

# # subset tissues as required
# gene_expression_summary <- gene_expression_summary[gene_expression_summary$tissue %in% c("endoderm", "eye", "glial", "mesenchyme", 
#                                                                                   "neural", "olfactory", "otic", "paraxial mesoderm", 
#                                                                                          "periderm", "pgc", "pronephros", "spinal cord"),]

# sqrt-scaled plot
ggplot(gene_expression_summary, aes(x = stage, y = counts_sum, fill = stage)) +
  geom_col() +
  scale_y_sqrt(breaks = c(0, 50, 100,200, 500, 1000, 1200)) +
  facet_wrap(~tissue, scales = "free_x") +
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

ggsave(paste0(gene, "_counts-sum-normalized_plot_sqrt.png"), height = 12, width = 15)


# linear plot
ggplot(gene_expression_summary, aes(x = stage, y = counts_sum, fill = stage)) +
  geom_col() +
  scale_y_continuous(breaks = c(0, 100,200, 500, 1000, 1200)) +
  facet_wrap(~tissue, scales = "free_x") +
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

ggsave(paste0(gene, "_counts-sum-normalized_plot_linear.png"), height = 12, width = 15)


##############################################################################



############################ Cell counts distribution ########################

# generate cell number factors to account for differences in cell abundance 
# between stages
stage_cell_nums <- table(daniocell@meta.data$stage.group)
stage_cell_factors <- as.vector(stage_cell_nums/41326)

stage_cell_factors_df <- data.frame(stage = c("  3-4","  5-6","  7-9"," 10-12", " 14-21", " 24-34", 
                " 36-46", " 48-58", " 60-70", " 72-82", " 84-94", " 96-106", "108-118", "120"), size_factors = stage_cell_factors)

##################### cell counts per tissue and stage #######################

# tissue level aggregation
gene_counts <- gene_data |> 
  group_by(tissue, stage) |> 
  summarise(cell_count = n()) |> 
  arrange(tissue, stage, .by_group = TRUE)


gene_counts <- inner_join(gene_counts , stage_cell_factors_df, by = "stage")

# adjust the counts by size factors
gene_counts$counts_norm <- gene_counts$cell_count/gene_counts$size_factors


# sorting by the sum of expression values
gene_counts <- gene_counts %>%
  group_by(tissue) %>%
  mutate(Total_sum = sum(counts_norm)) %>%
  arrange(desc(Total_sum))


# subset the tissues or lineages as required
# gene_counts <- gene_counts[gene_counts$tissue %in% c("endoderm", "eye", "glial", "mesenchyme", 
#                                                     "neural", "olfactory", "otic", "paraxial mesoderm", 
#                                                     "periderm", "pgc", "pronephros", "spinal cord"),]


############################### Visualization of the data ####################

# linear data plot
ggplot(gene_counts, aes(x = stage, y = counts_norm, fill = stage)) +
  geom_col() +
  scale_y_continuous(breaks = c(0,100, 200, 400, 600, 700))+
  facet_wrap(~tissue, scales = "free_x") +
  ggtitle(paste("Summarised cell counts for", gene , "in zebrafish"))+
  ylab("Cell count")+
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

ggsave(paste0(gene, "_cell_counts-normalized_plot.png"), height = 12, width = 15)
