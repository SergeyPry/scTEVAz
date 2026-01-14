# single-cell Tissue Expression Values Aggregator in zebrafish - scTEVAz

Multiple single-cell RNA sequencing (scRNAseq) atlases have been published for zebrafish and many other animal model species. Each such project has their own software tools, websites and interfaces to help make sense of the data. However, many of the visualizations for individual gene expression profiles are not very informative and concise enough to inform the researchers using these tools. Different atlases also understandably did not make their visualizations mutually consistent. Plots for quantitative expression comparisons across different classified cell types are also seriously lacking in many of these atlas projects.

Given these current limitations, we used three zebrafish scRNAseq datasets in a project to allow visualizing normalized expression values for all detectable genes in each dataset or gene-positive cell numbers across the tissues or cell types identified in that dataset. This allows cross-referencing information from different stages of development and adulthood across multiple datasets to obtain highly informative and reliable insights into tissue/lineage or expression patterns of genes, for which little experimental expression pattern evidence exists. However, these expression pattern plots should not be understood as the ultimate truth due to potential limitations of scRNAseq but rather as the guiding light for further experimental work.

## Description

An in-depth paragraph about your project and overview of use.

## Getting Started
To be able to reproduce the results described in this repository and project, a sufficiently modern and powerful PC or a laptop is strongly recommended with 64 GB of RAM and other appropriate specifications.


### Dependencies
All of the required packages and their installation are listed in the README files of the relevant sections of the project.


## Repository structure
This repository for the scTEVAz project describes all the relevant steps performed to obtain the publicly available zebrafish scRNAseq, calculate size factors for individual stages in each dataset, summarize the large datasets into more easily accessible data files for individual genes containing either normalized read count sums or normalized counts of gene-positive cells. These individual data files were further filtered to remove empty data frames or those containing only a single number.

- [01_Loading_raw_data_and_sources](https://github.com/SergeyPry/scTEVAz/tree/main/01_Loading_raw_data_and_sources) - Description and links to the original data sources that were used in this project as well as the code to read the data files for subsequent data extraction and processing.

- [02_Size_factors_calculation](https://github.com/SergeyPry/scTEVAz/tree/main/02_Size_factors_calculation) - Data processing code with explanatory comments to generate summary files for each stage of each dataset. These summary files will then be used to calculate size factors for normalizing the aggregated scRNAseq values.

- [03_Data_summarizing](https://github.com/SergeyPry/scTEVAz/tree/main/03_Data_summarizing) - Detailed description and code for generating normalized and aggregated expression values or cell counts data frames for all three datasets.

- [04_Data_cleanup](https://github.com/SergeyPry/scTEVAz/tree/main/04_Data_cleanup) - Detailed description and code for filtering data frames to keep only those that contain more than a single number. After this filtering, the information on available genes for each dataset was updated.

- [05_Standalone_plotting_scripts](https://github.com/SergeyPry/scTEVAz/tree/main/05_Standalone_plotting_scripts) - The scripts that work with the original large data and allow custom coding to generate plots. 



## Author
[Sergey Prykhozhij](https://github.com/SergeyPry)

