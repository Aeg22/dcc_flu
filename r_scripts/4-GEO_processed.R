

################################################################################
####### -------                     Set paths                    ------- #######

setwd("/Volumes/Partition 2/Core/degregori_continued_from_meher_scRNAseq_Oct2024/")
output <- "GEO/"
dir.create(output)

library(Seurat)
library(ggplot2)
library(plyr)

################################################################################
####### -------        Read in data (merge if necessary)         ------- #######

####### -------        No merging of this object
so <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_cluster_clean_labeled.RDS")
DefaultAssay(so) <- "RNA" # Only RNA has the counts
colnames(so@meta.data) # Will remove some columns later
# Add UMAP coordinates
umapCoord <- as.data.frame(Embeddings(object = so[["umap"]]))
so_meta <- so@meta.data # Will add subsequent metadata columns to this
orig_n <- ncol(so_meta) # Any columns in addition to this will be added back to seurat object
so_meta <- cbind.data.frame(cellbarcode = rownames(so_meta), so_meta)
so_meta <- cbind(so_meta, umapCoord)

# Some considerations
#   Not all cells will have complete labels (ex. T cell subtypes may be NA for other cells)
#   If I want to include multiple UMAP coordinates (ex. all data, T cell subset), I need to be clear

# Add new metadata back to seurat object (will still remove some columns later)
for(i in (orig_n+1):ncol(so_meta)){
  print(i)
  so@meta.data[,colnames(so_meta)[i]] <- so_meta[,i]
}
colnames(so@meta.data)


################################################################################
####### -------       Write counts matrix        ------- #######

####### -------        Primary
counts_matrix <- data.frame(as.matrix(GetAssayData(object = so, slot = "counts")))
counts_matrix <- cbind.data.frame(Gene = rownames(counts_matrix),
                                  counts_matrix)
dim(counts_matrix)
counts_matrix_sub <- counts_matrix[1:10, 1:10]
View(counts_matrix_sub)
write.table(counts_matrix, 
            gzfile(paste0(output, "counts_matrix.txt.gz")), 
            sep = '\t', 
            row.names = F, 
            col.names = T, 
            quote = F)

data <- data.table::fread(paste0(output, "counts_matrix.txt.gz"),
                          data.table = F)
rownames(data) <- data$Gene
datab_ReadIn <- data[1:10, 1:10]


################################################################################
####### -------                    metadata                      ------- #######

####### -------        Primary

metadata <- so@meta.data
colnames(metadata)

# Select columns to keep
colnames(metadata)
keep <- c("Sample", "Sample_num", "sample_id", "probe_barcode_ids", 
          "Group", "Broad.Group",
          "exp", "intervention1", "intervention2", "dpi", 
          "nCount_RNA", "nFeature_RNA", "percent_mito", 
          "S.Score", "G2M.Score", "Phase", "CC.Difference", 
          "seurat_clusters", "Cell.Type", "Broad.Cell.Type", 
          "UMAP_1", "UMAP_2")
keep %in% colnames(metadata)
metadata_sub <- metadata[,keep]

# Add Cell
metadata_sub <- cbind.data.frame(Cell = rownames(metadata_sub),
                                 metadata_sub)
write.table(metadata_sub, 
            paste0(output, "cell_metadata.txt"), 
            sep = '\t', 
            row.names = F, 
            col.names = T, 
            quote = F)

################################################################################
####### -------              Write sessionInfo                   ------- #######

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
