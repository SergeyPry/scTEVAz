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

# To add the latest annotations to the Seurat metadata
#pbmc@meta.data <- cbind(pbmc@meta.data, ZCL_metadata)

nrow(pbmc@meta.data)
#[1] 1088106

nrow(ZCL_metadata)
#[1] 1082680
# the metadata file is a bit smaller, which means it can be used for subsetting the larger file


###############################################################################
############################## ZCL dataset ############################## ###############################


################################## isolation of gene-specific data ############

# collect the data
assay_data <- pbmc@assays

gene <- "akna"

gene_names <- as.vector(assay_data$RNA@data@Dimnames[[1]])
gene_data <- as.data.frame(assay_data$RNA@data[gene,])


# complete the data frame
colnames(gene_data) <- c("counts")

sum(gene_data$counts > 0)
# [1] 3058

gene_data$cell_id <- rownames(gene_data)


gene_data$cluster <- ZCL_metadata[rownames(gene_data),]$cluster
gene_data$stage <- ZCL_metadata[rownames(gene_data),]$stage
gene_data$cell_type <- ZCL_metadata[rownames(gene_data),]$cell_type
gene_data$cell_lineage <- ZCL_metadata[rownames(gene_data),]$cell_lineage


# filter out the empty cells
gene_data <- gene_data[gene_data$counts >= 1,] 

nrow(gene_data)
# [1] 3058

# prevent scientific notation in R
options(scipen=999)

sum(is.na(gene_data$cell_type))
# [1] 2

sum(is.na(gene_data$cell_lineage))
# [1] 2

sum(is.na(gene_data$stage))
# [1] 2


# remove the NA data
gene_data <- gene_data[!is.na(gene_data$cell_type),]



############################### Sum of counts expression distribution #########
###############################################################################


# adjust the count values by size factors

size_factors_df <- data.frame(stage = c("21Day" ,  "22Month", "24hpf",  "3Month", "72hpf"),
                              size_factors = c(0.7470865, 1.8293085, 0.2447139, 3.1254952, 0.9347619))


gene_data <- inner_join(gene_data, size_factors_df, by = "stage")

# adjust the counts by size factors
gene_data$counts_norm <- gene_data$counts/gene_data$size_factors

# order the stages
gene_data$stage <- factor(gene_data$stage, levels = c("24hpf", "72hpf","21Day", "3Month","22Month"))


gene_expression_summary <- gene_data |> 
  group_by(cell_type, stage) |> 
  summarise(counts_sum = sum(counts_norm)) |> 
  arrange(cell_type, stage, .by_group = TRUE)


############################### Visualization of the data ####################
# 
library(ggplot2)


# redefine gene name if the official symbol is too weird
# gene <- "kmt2ca"

# subset tissues as required
gene_expression_summary <- gene_expression_summary[gene_expression_summary$cell_type %in% c("Epithelial cell", "Erythrocyte", "Hepatocyte", "Mesenchymal cell", 
                                                                                         "Neural cell", "Neural progenitor cell", "Osteoblast", "Primordial germ cell", 
                                                                                         "Radial glia", "Retinal cell", "Retinal cone cell", "Spermatocyte", "T cell"),]

ggplot(gene_expression_summary, aes(x = stage, y = counts_sum, fill = stage)) +
  geom_col() +
  scale_y_sqrt(breaks = c(0, 20, 50, 100,200, 300, 400)) +
  facet_wrap(~cell_type, scales = "free_x") +
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


ggsave(paste0(gene, "_counts-sum-normalized_plot_sqrt_ZCL.png"), height = 12, width = 15)


ggplot(gene_expression_summary, aes(x = stage, y = counts_sum, fill = stage)) +
  geom_col() +
  scale_y_continuous(breaks = c(0, 20, 50, 100,200, 300, 400)) +
  facet_wrap(~cell_type, scales = "free_x") +
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


ggsave(paste0(gene, "_counts-sum-normalized_plot_linear_ZCL.png"), height = 12, width = 15)



##############################################################################



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


# subset tissues as required
gene_counts <- gene_counts[gene_counts$cell_type %in% c("Epithelial cell", "Erythrocyte", "Hepatocyte", "Immune cell", "Keratinocyte", "Mesenchymal cell", 
                                                                                            "Neural cell", "Neural progenitor cell", "Osteoblast", "Primordial germ cell", 
                                                                                            "Radial glia", "Retinal cell", "Retinal cone cell", "Spermatocyte", "T cell"),]



############################### Visualization of the data ####################
ggplot(gene_counts, aes(x = stage, y = counts_norm, fill = stage)) +
  geom_col() +
  scale_y_continuous()+
  facet_wrap(~cell_type, scales = "free_x") +
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


ggsave(paste0(gene, "_cell-counts-normalized_plot_linear_ZCL.png"), height = 12, width = 15)

