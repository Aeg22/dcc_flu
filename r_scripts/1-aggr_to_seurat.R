
################################################################################
####### -------                   Goals                          ------- #######

# Read in cellranger aggr outs/
# Make a seurat object
# Add sample names to metadata

################################################################################
####### -------        Notes on this analysis                    ------- #######

# Normalization was not performed in cellranger aggr

################################################################################
####### -------        Read in R packages                        ------- #######

library(dplyr)
library(plyr)
library(Seurat)
library(openxlsx)

################################################################################
####### -------        Set paths and parameters                  ------- #######

setwd("/Volumes/Partition 2/Core/degregori_continued_from_meher_scRNAseq_Oct2024/") # Main project directory

cellranger_outs_path <- "Data/Meher/cellranger_16samples_run/Aggr_noNorm_DeGregori_16samples/" # Location of the cellranger outs directory 
seurat_outputs <- "Data/Seurat_objects/"
dir.create(seurat_outputs, recursive = T)

min.cells <- 10
min.features <- 100

################################################################################
####### -------                Read in data                      ------- #######

# Load the dataset
data.10x <- Read10X(data.dir = paste0(cellranger_outs_path,
                                       "/filtered_feature_bc_matrix/"))

# Initialize the Seurat object with the raw (non-normalized data)
so <- CreateSeuratObject(counts = data.10x, 
                           min.cells = min.cells, 
                           min.features = min.features)

################################################################################
####### -------        Add sample names to metadata              ------- #######

# You can see that the samples are labels -1, -2, etc.
# This reflects the order of the samples in the CSV used in cellranger aggr
# Read in this CSV to add the sample info
aggr <- read.table(file = paste0(cellranger_outs_path,
                                 "/aggregation.csv"),
                   sep = ",", header = T, stringsAsFactors = F)
aggr$Sample_num <- 1:nrow(aggr)

barcode_order <- cbind.data.frame(Full = rownames(so@meta.data),
                                  Sample_num = sapply(strsplit(rownames(so@meta.data),
                                                               split = "-"), function(x){
                                                                 x[[2]]
                                                               }))
table(barcode_order$Sample_num)
barcode_order <- join(barcode_order,
                      aggr)
table(barcode_order$sample_id) 

# Add more metadata
samples <- read.xlsx("helper_data/samples.xlsx")
barcode_order <- join(barcode_order, samples)
barcode_order <- barcode_order[,-4]

# Add metadata to seurat object
for(i in 2:ncol(barcode_order)){
  so <- AddMetaData(so, 
                    barcode_order[,i], 
                    col.name = colnames(barcode_order)[i])
}

so_meta <- so@meta.data

# Filter out samples removed from the analysis
so
table(so$Group)
so <- subset(so, 
             Group %in% c("excluded_sample"),
             invert = T)
so

# Write RDS object
saveRDS(so, paste0(seurat_outputs, "/AggrNoNorm_minFiltered.RDS"))

