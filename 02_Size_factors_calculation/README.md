# single-cell Tissue Expression Values Aggregator in zebrafish - scTEVAz

## Size factors calculation for count sums
scRNA-seq datasets are typically very sparse meaning that each cell contains expression values for relatively small numbers of genes. Therefore, to obtain total expressions values for each gene, it is necessary to first aggregate these values for each gene across all cells of each scRNA-seq library. After such a matrix is obtained, size factors can be easily calculated using standard functions used for RNA-seq. In summary, the overall algorithm consists of the following steps:

1. Aggregate the values across all cells and genes.
2. Calculate size factors using DESeq2.

### 1. Aggregate the values across all cells and genes.
This step generates a data frame containing sums of counts for each dataset. We provide the code for summarizing each dataset.

**Scripts for the initial count sums aggregation:**

1. [Daniocell_data_aggregate.R](https://github.com/SergeyPry/scTEVAz/tree/main/02_Size_factors_calculation/Daniocell_data_aggregate.R).
2. [ZCL_data_aggregate.R](https://github.com/SergeyPry/scTEVAz/tree/main/02_Size_factors_calculation/ZCL_data_aggregate.R).
3. [Zebrahub_data_aggregate.R](https://github.com/SergeyPry/scTEVAz/tree/main/02_Size_factors_calculation/Zebrahub_data_aggregate.R).

This first set of scripts is very time-consuming to run and requires that you load the data as described in the[Data loading section](https://github.com/SergeyPry/scTEVAz/tree/main/01_Loading_raw_data_and_sources/).

### 2. Calculate size factors using DESeq2.
This step takes the output of the previous step and calculates the size factors for each dataset. The required files are included in the same folder as the scripts.

**Scripts for size factors calculation	:**

1. [size_factors_deseq2_Daniocell.R](https://github.com/SergeyPry/scTEVAz/tree/main/02_Size_factors_calculation/size_factors_deseq2_Daniocell.R).
2. [size_factors_deseq2_ZCL.R](https://github.com/SergeyPry/scTEVAz/tree/main/02_Size_factors_calculation/size_factors_deseq2_ZCL.R).
3. [size_factors_deseq2_zebrahub.R](https://github.com/SergeyPry/scTEVAz/tree/main/02_Size_factors_calculation/size_factors_deseq2_zebrahub.R).



