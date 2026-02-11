# single-cell Tissue Expression Values Aggregator in zebrafish - scTEVAz

## Standalone plotting scripts
The scripts provided in this folder visualize the data in Daniocell, Zebrafish Cell Landscape and Zebrahub in an analogous way to the main app, but are much easier to customize to specific needs. They are built on some of the previous parts of the repository, but do not use the small files for individual genes that are used for the app itself. All of the scripts have the following structure:

**The Algorithm:**
1. Dataset loading.
2. Isolation of gene-specific data including both expression data and cell metadata.
3. Converting the expression data to count-like data.
4. Removing cells without known tissue, cell type or lineage assignment.

The next steps are performed separately for count-sum and cell-counts data 
5. Normalize the data by the size factors of the stages.
6. Aggregate the data at the levels of tissue, cell type or lineage.
7. Sort the data by the sum of expression values.
8. Plot the data in the square-root and linear scales for count-sum data and linear-only for the cell-counts data.

**Scripts:**
The scripts can be viewed and downloaded using the links below:
* [Daniocell script](https://github.com/SergeyPry/scTEVAz/tree/main/05_Standalone_plotting_scripts/daniocell_visualize_gene.R)
* [ZCL script](https://github.com/SergeyPry/scTEVAz/tree/main/05_Standalone_plotting_scripts/ZCL_visualize_gene.R)
* [Zebrahub script](https://github.com/SergeyPry/scTEVAz/tree/main/05_Standalone_plotting_scripts/zebrahub_visualize_gene.R)





