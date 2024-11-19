# Description

This repository contains the main scripts used in the analysis of the metastatic mouse lung scRNA-seq data from the publication *Respiratory viral infection promotes the awakening and outgrowth of dormant metastatic breast cancer cells in lungs * (PMID: XXXX). Processed data are deposited to GEO as GSE264175.

# Code overview

Code is categorized and described in this README. 

### Cell Ranger

Scripts utilized during Cell Ranger alignment and counting. Aggregation of samples was performed without normalization. 

Shell scripts:

* `cellranger/cellranger_count.sh`
* `cellranger/cellranger_aggregation.sh`

### Seurat initialization, quality control, and filtering

Generation of initial Seurat object, exploration of data quality, and filtering genes and cells. 

R scripts:

* `r_scripts/1-aggr_to_seurat.R`
* `r_scripts/2-filter.R`

### Downstream analysis

The cell data was normalized, clustered, and cell types were identified. Downstream analysis involved cell proportions, differential expression, and pathway analysis. 

R scripts:

* `r_scripts/3-norm_cluster_celltype.R`

### Data availability

Scripts to process data for GEO submission. 

R scripts:

* `r_scripts/4-GEO_processed.R`

# Releases and changes

* 1.0
  * Initial Release
