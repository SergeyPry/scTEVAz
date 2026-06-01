# Load the necessary library
library(dplyr)

setwd("c:/Bioinformatics/00_Daniocell_data/scTEVAz/Daniocell_dataset")


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
 
  
  num_table_cells <- (ncol(df_main) - 1) * nrow(df_main)
  
  newpath <- paste0("./filtered/", filename)
  
  if(num_table_cells > 1){
    write.csv(df_main, newpath, quote = FALSE, row.names = FALSE)
  }
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
  
  num_table_cells <- (ncol(df_main) - 1) * (nrow(df_main) - 1) 
  
  newpath <- paste0("./filtered/", filename)
  
  if(num_table_cells > 1){
    write.csv(df_main, newpath, quote = FALSE, row.names = FALSE)
  }
  
}

