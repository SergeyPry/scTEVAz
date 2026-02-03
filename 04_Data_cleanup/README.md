# single-cell Tissue Expression Values Aggregator in zebrafish - scTEVAz

## Data cleanup

**Detection of unreadable files and removal of very small files**
To perform quality control on the summarized data, we first tested whether each data file can be read and if the resulting data frame contains at least 2 data points. The script [check_reading_and_dims.R](https://github.com/SergeyPry/scTEVAz/tree/main/04_Data_cleanup/check_reading_and_dims.R) is also shown below:

```r
# Load the necessary library
library(dplyr)

# uncomment and choose which folder needs this processing step

#setwd("c:/Bioinformatics/00_Daniocell_data/scTEVAz/Daniocell_dataset")


################### count sums data ##########################################

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

################### cell counts data ##########################################

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
  
  num_table_cells <- (ncol(df_main) - 1) * nrow(df_main)  
  
  newpath <- paste0("./filtered/", filename)
  
  if(num_table_cells > 1){
    write.csv(df_main, newpath, quote = FALSE, row.names = FALSE)
  } 
}

```

After this script has finished running, we checked the warnings in RStudio and the number of files which passed the filtering step. In case of significant problems, the summarized data files were regenerated and the above script run again. The successfully filtered files were then copied to their respective folders after removing the original files.

**Updating of the available gene catalog for each dataset**
Since the above procedure reduces the number of genes whose pattern of tissue expression can be visualized by this tool, it is necessary to update the file, from which the set of available genes is derived. We use the following algorithm to achieve this goal:

1. Read filenames from both count-sums and cell-counts folders.
2. Extract processed gene names from both and make a union from both of them.
3. Subset the gene mapping data frame to this union of processed genes.
4. Store this gene mapping data frame for future use.

The script [extract_subset_names.R](https://github.com/SergeyPry/scTEVAz/tree/main/04_Data_cleanup/extract_subset_names.R) implements this algorithm for one of the ZCL datasets. For all other datasets, the code is identical except for the filenames of the gene tables.

```r
library(stringr)

# Algorithm:

# 1. Read filenames from both count-sums and cell-counts folders.
# 2. Extract processed gene names from both and make a union from both of them.
# 3. subset the gene mapping data frame to this union of processed genes.
# 4. Store this gene mapping data frame for future use.


# original data frame with both original and processed gene names
gene_df_filename <- "ZCL_map_df.csv"  # others: daniocell_map_df.csv, zhub_map_df.csv 

updated_df_filename <- "ZCL_map_celltypes_df.csv" # others: daniocell_map_df_updated.csv,
# ZCL_map_lineages_df.csv, zhub_map_df_updated.csv


# 1. Read filenames from both count-sums and cell-counts folders.

countsum_file_names <- list.files("./single_genes_dfs/")
cellcount_file_names <- list.files("./single_genes_cell_counts/")

all_file_names <- c(countsum_file_names, cellcount_file_names)

# 2. Extract processed gene names from both and make a union from both of them.

# set up and collect all of the available processed file names 
processed_names <- c()


for(filename in all_file_names){
  
  items <- str_split(filename, '_')  
  
  processed_names <- c(processed_names, items[[1]][2]) 
  
}

# take care of the processed name duplication
processed_names <- unique(processed_names)

# 3. subset the gene mapping data frame to this union of processed genes.

# read the data frame with all gene names
map_df <- read.csv(gene_df_filename, header = TRUE)

# subset 
map_df_updated <- map_df[map_df$proc_gene %in% processed_names,] 
  

# 4. Store this gene mapping data frame for future use.
write.csv(map_df_updated, file = updated_df_filename,  quote = FALSE, row.names = FALSE)

```
The resulting updated data frames were used in the app.