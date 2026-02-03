# Load DESeq2 package
library(DESeq2)

#setwd("c:/Users/User/Documents/GitHub/scTEVAz/02_Size_factors_calculation")

readCountData <- function(filename, sampTab){
  count_data <- read.csv(filename, sep = "\t", header = TRUE)
  rownames(count_data) <- count_data$gene
  count_data <- count_data[, 2:ncol(count_data)]
  colnames(count_data) <- rownames(sampTab)
  
  return(count_data)
}


# # load and process the sample table - custom for each table
sampleTable <- read.csv("sample_table_Daniocell.txt", sep = "\t")
rownames(sampleTable) <- sampleTable$sampleName

# # load the raw count data
count_data <- readCountData("daniocell_count_sums.csv", sampleTable)

saved_rn <- rownames(count_data)
count_data <- as.data.frame(sapply(count_data, ceiling))
rownames(count_data) <- saved_rn

head(count_data)

# # 1. Read in the count data into a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sampleTable,
                              design = ~condition)
nrow(dds)
# 36250

dds <- dds[rowSums(counts(dds) == 0) < 6, ]
dds <- dds[rowSums(counts(dds)) > 28, ]

nrow(dds)
#[1] 30346

dds <- estimateSizeFactors(dds)

size_factors <- sizeFactors(dds)
print(size_factors)

# wt3-4     wt5-6     wt7-9   wt10-12   wt14-21   wt24-34   wt36-46   wt48-58   wt60-70   wt72-82   wt84-94  wt96-106 wt108-118     wt120 
# 0.1140180 0.1757729 0.3324738 0.3544112 1.3264766 1.7987805 1.9225959 2.8822294 1.8920302 2.1923608 1.9349128 1.9104694 2.5820813 1.0633326


