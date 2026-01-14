# Load the necessary library
library(dplyr)

setwd("c:/Bioinformatics/00_Daniocell_data/scTEVAz/Zebrahub_dataset")


# Mapping data frame 
df_map <- read.csv("stages_map.txt", sep = "\t")


# 1. Create Example Data
# Main data frame (unsorted)

all_countsum_files <- list.files("./single_genes_dfs/")

for(filename in all_countsum_files){
  print(filename)
  
  filepath <- paste0("./single_genes_dfs/", filename)
  
  skip_to_next <- FALSE # Initialize flag for each iteration
  
  
  # reading the file
  problemPlace <- tryCatch({
    
    df_main <- read.csv(filepath)
    
  }, error = function(e) {
    skip_to_next <<- TRUE # Set the flag using super-assignment
  })
  
  # if reading the file failed
  if (skip_to_next) {
    next # Skip the rest of the loop body
  }
 

  # 2. Map and Order
  result <- df_main %>%
    # Join the data frames to map the new column
    left_join(df_map, by = "stage") %>%
    # Convert 'code' to a factor using the order found in df_map
    mutate(stage = factor(stage, levels = df_map$stage)) %>%
    mutate(timepoint = factor(timepoint, levels = df_map$timepoint)) %>%
    # Sort the data frame based on that factor order
    arrange(timepoint)
  
  # View the result
  result$stage <- result$timepoint
  
  result <- result[, 1: ncol(df_main)]
  
  write.csv(result, filepath, quote = FALSE, row.names = FALSE)
  
}


#################################

all_cellcount_files <- list.files("./single_genes_cell_counts/")

for(filename in all_cellcount_files){
  print(filename)
  
  filepath <- paste0("./single_genes_cell_counts/", filename)
  
  skip_to_next <- FALSE # Initialize flag for each iteration
  
  
  # reading the file
  problemPlace <- tryCatch({
    
    df_main <- read.csv(filepath)
    
  }, error = function(e) {
    skip_to_next <<- TRUE # Set the flag using super-assignment
  })
  
  # if reading the file failed
  if (skip_to_next) {
    next # Skip the rest of the loop body
  }
  
  # 2. Map and Order
  result <- df_main %>%
    # Join the data frames to map the new column
    left_join(df_map, by = "stage") %>%
    # Convert 'code' to a factor using the order found in df_map
    mutate(stage = factor(stage, levels = df_map$stage)) %>%
    mutate(timepoint = factor(timepoint, levels = df_map$timepoint)) %>%
    # Sort the data frame based on that factor order
    arrange(timepoint)
  
  # View the result
  result$stage <- result$timepoint
  
  result <- result[, 1: ncol(df_main)]
  
  write.csv(result, filepath, quote = FALSE, row.names = FALSE)
  
}







