# Load DESeq2 package
library(DESeq2)

setwd("c:/Bioinformatics/00_Daniocell_data/size_factors_data/")

readCountData <- function(filename, sampTab){
  count_data <- read.csv(filename, sep = "\t", header = TRUE)
  rownames(count_data) <- count_data$gene
  count_data <- count_data[, 2:ncol(count_data)]
  colnames(count_data) <- rownames(sampTab)
  
  return(count_data)
}

# Assume 'countData' is your raw count matrix and 'colData' contains metadata
# Example:
# countData <- matrix(rnbinom(n=1000, mu=100, size=100), ncol=10)
# rownames(countData) <- paste0("gene", 1:100)
# colData <- data.frame(condition = factor(rep(c("A", "B"), each=5)))
# rownames(colData) <- paste0("sample", 1:10)


# # load and process the sample table - custom for each table
sampleTable <- read.csv("sample_table_zebrahub.txt", sep = "\t")
rownames(sampleTable) <- sampleTable$sampleName

# # load the raw count data
count_data <- readCountData("zebrahub_count_sums.csv", sampleTable)

saved_rn <- rownames(count_data)
count_data <- as.data.frame(sapply(count_data, ceiling))
rownames(count_data) <- saved_rn

head(count_data)

# # 1. Read in the count data into a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sampleTable,
                              design = ~condition)
nrow(dds)
# 32060

dds <- dds[rowSums(counts(dds) == 0) < 5, ]
dds <- dds[rowSums(counts(dds)) > 20, ]

nrow(dds)
#[1] 26091

dds <- estimateSizeFactors(dds)

size_factors <- sizeFactors(dds)
print(size_factors)

# 0 somites   05 somites   10 somites   15 somites   20 somites   30 somites  larval-2dpf 
# 0.1208711    0.5704448    0.4491545    0.6201209    0.6681976    1.2310201    2.1067936 
# larval-3dpf  larval-5dpf larval-10dpf 
# 3.5070977    3.5066438    2.7098503 

# 0.1208711    0.5704448    0.4491545    0.6201209    0.6681976    1.2310201    2.1067936 3.5070977    3.5066438    2.7098503 
