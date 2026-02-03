# single-cell Tissue Expression Values Aggregator in zebrafish - scTEVAz

## Data summarizing
This section contains the scripts summarizing the data for each gene. These scripts are most central to data preparation for this project. The reason why it is important to summarize for each gene is that the complete scRNA-seq datasets are very large and may not be suitable for loading into apps directly. Therefore, it is much more realistic and faster to generate individual files for all genes available in each dataset and then load them whenever they are requested. There are two types of data summaries: count sums and cell counts of gene-expressing cells. The following are the basic algorithms for generating these summary files:

1. Load the whole dataset.
2. Subset the large dataset by gene name and make sure gene name is processed for file name generation.
3. Store the cell metadata together with expression values for individual cells for the selected gene.
4. Exponentiate the expression values to obtain count-like values.
5. Filter out the cells without expression of the selected gene and where tissue or lineage metadata is absent.
6. Divide the resulting values by size factors for each stage. Size factors for count sums and cell counts are somewhat different.
7. Compute summaries at the level of count sums and cell counts.
8. Store the resulting files for future retrieval.

### Scripts for data summarizing
The following scripts implement the algorithms described above and generate the data for the website:

1. [Daniocell_tissue_summarizer.R](https://github.com/SergeyPry/scTEVAz/tree/main/03_Data_summarizing/Daniocell_tissue_summarizer.R).
2. [ZCL_celltypes_data_summarizer.R](https://github.com/SergeyPry/scTEVAz/tree/main/03_Data_summarizing/ZCL_celltypes_data_summarizer.R).
3. [ZCL_lineage_data_summarizer.R](https://github.com/SergeyPry/scTEVAz/tree/main/03_Data_summarizing/ZCL_lineage_data_summarizer.R).
4. [Zebrahub_lineage_data_summarizer.R](https://github.com/SergeyPry/scTEVAz/tree/main/03_Data_summarizing/Zebrahub_lineage_data_summarizer.R).
 
This set of scripts is very time-consuming to run and requires that you load the data as described in the[Data loading section](https://github.com/SergeyPry/scTEVAz/tree/main/01_Loading_raw_data_and_sources/).

