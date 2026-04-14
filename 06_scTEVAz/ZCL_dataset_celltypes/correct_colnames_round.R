# Load the necessary library
library(dplyr)
library(stringr)

setwd("c:/Bioinformatics/00_Daniocell_data/scTEVAz/ZCL_dataset_celltypes")


all_countsum_files <- list.files("./single_genes_dfs/")

for(filename in all_countsum_files){
  print(filename)
  
  filepath <- paste0("./single_genes_dfs/", filename)
  
  df_main <- read.csv(filepath)
  
  colnames_new <- str_replace_all(colnames(df_main), '\\.', '_')
  
  colnames(df_main) <- colnames_new
  
  df_rounded <- df_main %>%
    mutate(across(-stage, ~round(.x, digits = 2)))

  write.csv(df_rounded, filepath, quote = FALSE, row.names = FALSE)

}


#################################

all_cellcount_files <- list.files("./single_genes_cell_counts/")

for(filename in all_cellcount_files){
  print(filename)
  
  filepath <- paste0("./single_genes_cell_counts/", filename)

  df_main <- read.csv(filepath)
  
  colnames_new <- str_replace_all(colnames(df_main), '\\.', '_')
  
  colnames(df_main) <- colnames_new
  
  df_rounded <- df_main %>%
    mutate(across(-stage, ~round(.x, digits = 2)))
  
  write.csv(df_rounded, filepath, quote = FALSE, row.names = FALSE)
  
}

