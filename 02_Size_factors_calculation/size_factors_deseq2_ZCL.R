# Load DESeq2 package
library(DESeq2)

#setwd("c:/Bioinformatics/00_Daniocell_data/size_factors_data/")

readCountData <- function(filename, sampTab){
  count_data <- read.csv(filename, sep = "\t", header = TRUE)
  rownames(count_data) <- count_data$gene
  count_data <- count_data[, 2:ncol(count_data)]
  colnames(count_data) <- rownames(sampTab)
  
  return(count_data)
}


# # load and process the sample table - custom for each table
sampleTable <- read.csv("sample_table_ZCL.txt", sep = "\t")
rownames(sampleTable) <- sampleTable$sampleName

# # load the raw count data
count_data <- readCountData("ZCL_count_sums.csv", sampleTable)

saved_rn <- rownames(count_data)
count_data <- as.data.frame(sapply(count_data, ceiling))
rownames(count_data) <- saved_rn

head(count_data)

# # 1. Read in the count data into a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sampleTable,
                              design = ~condition)

dds <- dds[rowSums(counts(dds) == 0) < 2, ]
dds <- dds[rowSums(counts(dds)) > 20, ]


dds <- estimateSizeFactors(dds)

size_factors <- sizeFactors(dds)
print(size_factors)

# 21Day   22Month     24hpf    3Month     72hpf 
# 0.7470865 1.8293085 0.2447139 3.1254952 0.9347619 