

################################################################################
####### -------                   Goals                          ------- #######

# Explore data quality
# Filter the Seurat object

################################################################################
####### -------        Notes on this analysis                    ------- #######
# Final filtering: 
    # # Remove cells with <200 detected genes
    # # Remove cells with >7500 detected genes (redundant with the UMI filter)
    # # Remove cells with >2.5% mitochondrial reads
    # # Remove cells with < 200 UMIs
    # # Remove cells with > 7500 UMIs
    # # Remove cells predicted to be doublets

umi_low_cut <- 200
umi_high_cut <- 7500
mt_high_cut <- 2.5
gene_low_cut <- 200
gene_high_cut <- 7500

# Vectors for plotting

# Other plotting params
gg_title_text <- 12
gg_axis_text <- 12

# Vectors for plotting
umi_cut_vector <- c(umi_low_cut, umi_high_cut)
umi_cut_vector[is.na(umi_cut_vector)] <- 0

gene_cut_vector <- c(gene_low_cut, gene_high_cut)
gene_cut_vector[is.na(gene_cut_vector)] <- 0

# Other plotting params
gg_title_text <- 12
gg_axis_text <- 12


################################################################################
####### -------                Read in data                      ------- #######
library(dplyr)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(scDblFinder)
library(plyr)

source("functions/SeuratQC.R")

setwd("/Volumes/Partition 2/Core/degregori_continued_from_meher_scRNAseq_Oct2024/")
sobj <- readRDS("Data/Seurat_objects/AggrNoNorm_minFiltered.RDS")

output <- "Analysis/preFilter_QC"
dir.create(output)


################################################################################
####### -------              Prefiltering table                  ------- #######

# Note that very minimal filtering occurred within aggr_to_seurat.R

# Dataset level
table_prefilter <- cbind.data.frame(Metric = c("Number of barcodes", 
                                             "Number of genes",
                                             "Median UMI per cell",
                                             "Median genes per cell"),
                                  Value = c(ncol(sobj),
                                            nrow(sobj),
                                            median(sobj$nCount_RNA),
                                            median(sobj$nFeature_RNA)))
write.table(table_prefilter,
            file = paste0(output, "/prefilter_dataset-level_metrics.txt"),
            sep = "\t", row.names = F, quote = F)

# Sample-level
table_sample <- cbind.data.frame(Sample = c("Number of barcodes", 
                                            "Number of genes",
                                            "Median UMI per cell",
                                            "Median genes per cell"))
for(samp in unique(sobj$Sample)){
  print(samp)
  
  sobj_samp <- subset(sobj,
                      Sample == samp)
  sobj_samp <- sobj_samp[rowSums(sobj_samp) > 0,]
  
  table_sample$samp <- c(ncol(sobj_samp),
                         nrow(sobj_samp),
                         median(sobj_samp$nCount_RNA),
                         median(sobj_samp$nFeature_RNA))
  colnames(table_sample)[ncol(table_sample)] <- samp
}

table_sample <- data.frame(t(table_sample))
write.table(table_sample,
            file = paste0(output, "/prefilter_sample-level_metrics.txt"),
            sep = "\t", row.names = T, col.names = F, quote = F)

################################################################################
####### -------         Add metadata including % MT              ------- #######

# Calculate percentage of mitochondrial reads for each cell
# This calculates the percentage of counts (should be using nCount_RNA)

rownames(sobj)[grep("^mt-", rownames(sobj))]

sobj <- sobj %>%
  PercentageFeatureSet(
    pattern  = "^mt-", # Look for genes that start with MT
    col.name = "percent_mito" # The name of the new column in metadata
  ) 


################################################################################
####### -------                  Plotting                        ------- #######

gg_umi <- sobj %>%
  VlnPlot(
    features = c("nCount_RNA"),
    cols = rep("grey", 12),
    ncol     = 1,
    pt.size  = 0,
    group.by = "Sample"#,
  ) + geom_hline(yintercept=umi_cut_vector,
                 linetype = 2) + 
  ylab("# of UMIs") + ggtitle("Number of UMIs") + xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text=element_text(size = gg_title_text),
        axis.text = element_text(size = gg_axis_text),
        legend.position="none")
gg_umi


sobj$log2UMI <- log2(sobj$nCount_RNA+1)

gg_genes <- sobj %>%
  VlnPlot(
    features = c("nFeature_RNA"),
    cols = rep("grey", 12),
    ncol     = 1,
    pt.size  = 0,
    group.by = "Sample"#,
  ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            text=element_text(size = gg_title_text),
            axis.text = element_text(size = gg_axis_text),
            legend.position="none") + 
  ylab("# of genes") + ggtitle("Number of genes") + xlab("") +
  geom_hline(yintercept=gene_cut_vector,
             linetype = 2)
gg_genes

gg_MT <- sobj %>%
  VlnPlot(
    features = c("percent_mito"),
    cols = rep("grey", 12),
    ncol     = 1,
    pt.size  = 0,
    group.by = "Sample"#,
  )+ geom_hline(yintercept=mt_high_cut,
                linetype = 2) + 
  ylab("UMI (%)") + ggtitle("UMI from mitochondria") + xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text=element_text(size = gg_title_text),
        axis.text = element_text(size = gg_axis_text),
        legend.position="none")
gg_MT

gg_umi_vs_genes <- FeatureScatter(object = sobj, 
                                  feature1 = 'nCount_RNA', 
                                  feature2 = 'nFeature_RNA',
                                  group.by = "Sample",
                                  pt.size = 0.1)
gg_umi_vs_genes


gg_MT_vs_UMI <- FeatureScatter(object = sobj, 
                               feature2 = 'percent_mito', 
                               feature1 = 'nCount_RNA',
                               group.by = "Sample",
                               pt.size = 0.1)

gg_MT_vs_UMI 



#########
# https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html

metadata <- sobj@meta.data

metadata$log10GenesPerUMI <- log10(metadata$nFeature_RNA) / log10(metadata$nCount_RNA)
sobj$log10GenesPerUMI <- log10(sobj$nFeature_RNA) / log10(sobj$nCount_RNA)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
library(plotly)
hist_GenesPerUMI <- metadata %>% 
  ggplot(aes(color=Sample, x=log10GenesPerUMI, fill= Sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  ylab("Log10 Cell density") 
ggplotly(hist_GenesPerUMI ) 

# Visualize the number UMIs/transcripts per cell
hist_UMI <- metadata %>% 
  ggplot(aes(color=Sample, x=nCount_RNA, fill= Sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Log10 UMI density") +
  geom_vline(xintercept = umi_cut_vector)
ggplotly(hist_UMI)

# Visualize the distribution of genes detected per cell via histogram
hist_genes <- metadata %>% 
  ggplot(aes(color=Sample, x=nFeature_RNA, fill= Sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = gene_cut_vector)
ggplotly( hist_genes) 

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent_mito)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  facet_wrap(~Sample)

# Visualize the distribution of mitochondrial gene expression detected per cell
hist_MT <- metadata %>% 
  ggplot(aes(color=Sample, x=percent_mito, fill=Sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  geom_vline(xintercept = mt_high_cut)
ggplotly(hist_MT)


#####################################################################
#####    Save pre-filtering plots
ggsave(paste0(output, "/Genes.pdf"),
       gg_genes)
ggsave(paste0(output, "/UMIs.pdf"),
       gg_umi)
ggsave(paste0(output, "/MT.pdf"),
       gg_MT)
ggsave(paste0(output, "/Hist_UMIs.pdf"),
       hist_UMI)
ggsave(paste0(output, "/Hist_MT.pdf"),
       hist_MT)
ggsave(paste0(output, "/Hist_genes.pdf"),
       hist_genes)
ggsave(paste0(output, "/Hist_GenesPerUMI.pdf"),
       hist_GenesPerUMI)

# Save an aggregated plot for ease of presentation
aggr_vio <- grid.arrange(gg_umi,
                         gg_genes,
                         gg_MT,
                         ncol = 3, nrow = 1)
ggsave(paste0(output, "/Aggregated_violin_plots.pdf"),
       aggr_vio,
       width = 16, height = 4)

prefilter_gg_umi <- gg_umi
prefilter_gg_genes <- gg_genes
prefilter_gg_MT <- gg_MT

aggr_hist <- grid.arrange(  hist_UMI,
                            hist_genes,
                            hist_MT,
                            ncol = 1, nrow = 3)
ggsave(paste0(output, "/Aggregated_histograms.pdf"),
       aggr_hist,
       width = 6, height = 8)


################################################################################
####### -------                Doublet finder                    ------- #######

## Include all samples for fixed assay

# Issue creating SingleCellExperiment object using the full dataset so I will
# subset by sample, convert each to SCE, and then merge before scdblfinder
sobj_type.list <- SplitObject(sobj, split.by = "Sample")
sce_list <- list()
for(i in 1:length(sobj_type.list)){ 
  print(i)
  so.sce <- as.SingleCellExperiment(sobj_type.list[[i]])
  sce_list[[length(sce_list) + 1]] <- so.sce
}
# Merge sample-specific sce objects
library("scMerge")
sce <- sce_cbind( sce_list, method = "intersect", 
                  exprs = c("counts"), 
                  colData_names = T)
# Run dblfinder
so.sce <- scDblFinder(sce)
scDblFinder_results <- cbind.data.frame(Barcode = colnames(so.sce),
                                        scDblFinder.class = so.sce$scDblFinder.class,
                                        scDblFinder.score = so.sce$scDblFinder.score)
table(scDblFinder_results$scDblFinder.class)

# Compare to ind-sample run (these are metrics and cluster filtered with cluster labels)
so_prevDbl <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_cluster_temp.RDS")
so_prevDbl_dbl <- cbind.data.frame(Barcode = colnames(so_prevDbl),
                                   prev_cluster = so_prevDbl$seurat_clusters,
                                   prev_class = so_prevDbl$scDblFinder.class,
                                   prev_score = so_prevDbl$scDblFinder.score)
merge <- join(scDblFinder_results, so_prevDbl_dbl,
              type = "full")
table(merge$scDblFinder.class)
table(merge$prev_class)
new_dblts <- merge[merge$scDblFinder.class == "doublet" & merge$prev_class == "singlet",]
new_dblts <- na.omit(new_dblts)
table(new_dblts$prev_cluster)

scDblFinder_results <- scDblFinder_results[match(colnames(sobj),
                                                 scDblFinder_results$Barcode),]
sobj$scDblFinder.class <- scDblFinder_results$scDblFinder.class
sobj$scDblFinder.score <- scDblFinder_results$scDblFinder.score

table(sobj$scDblFinder.class,
      sobj$Sample)


################################################################################
####### -------                Filter cells                      ------- #######

sobj_metrics_filtered <- sobj %>%
  subset(
    nCount_RNA > umi_low_cut &
      nCount_RNA <= umi_high_cut & 
      nFeature_RNA >= gene_low_cut &   
      nFeature_RNA <= gene_high_cut & 
      percent_mito <= mt_high_cut # & 
    # log10GenesPerUMI > 0.80
  )

table(sobj_metrics_filtered$scDblFinder.class,
      sobj_metrics_filtered$Sample)

VlnPlot(sobj_metrics_filtered,
        features = "scDblFinder.score",
        group.by = "scDblFinder.class")


################################################################################
####### -------   Filter out clusters with high #'s of doublets  ------- #######

# For this dataset, I am filtering out clusters with a high percentage of doublets
# prior to filtering out all doublet cells

# Will repeat if multiple rounds are necessary

threshold <- 0.5 # Will remove clusters with > 50% doublets

## Round 1
so_round1 <- sobj_metrics_filtered
DefaultAssay(so_round1)
so_round1 <- NormalizeData(so_round1, normalization.method = "LogNormalize")
so_round1 <- FindVariableFeatures(so_round1, selection.method = "vst")
so_round1 <- ScaleData(object = so_round1)
so_round1 <- RunPCA(so_round1,
             features = VariableFeatures(so_round1),
             verbose = F,
             ndims.print = 0)
ElbowPlot(so_round1, ndims = 40)
n_dim <- 30 
so_round1 <- FindNeighbors(so_round1, 
                    reduction = "pca",
                    dims = 1:n_dim, 
                    k.param = 20) 
resolutionsCalculated <- c(  3, 4)
so_round1 <- FindClusters(so_round1, resolution = resolutionsCalculated) 
so_round1_meta <- so_round1@meta.data
for(i in resolutionsCalculated){ # Want a resolution with a high number of clusters
  print(paste0("res:", i))
  print(length(unique(so_round1_meta[,paste0("RNA_snn_res." , i)])))
}
so_round1$seurat_clusters <- so_round1$RNA_snn_res.4
prop_df_sample <- data.frame(prop.table(table(so_round1@meta.data[,"scDblFinder.class"], 
                                              so_round1@meta.data[,"seurat_clusters"]), margin = 2))
prop_df_sample_doublets <- prop_df_sample[prop_df_sample$Var1 == "doublet",]
gg_stackbar_sample <- ggplot(prop_df_sample, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +ylab("Proportion") + xlab("New cluster") +
  theme_classic()
gg_stackbar_sample
remove <- as.character(prop_df_sample_doublets$Var2[prop_df_sample_doublets$Freq > threshold])

# FILTER round 1
so_round1_post <- subset(so_round1,
                         seurat_clusters %in% remove,
                         invert=T)

## Round 2
so_round2 <- so_round1_post
DefaultAssay(so_round2)
so_round2 <- NormalizeData(so_round2, normalization.method = "LogNormalize")
so_round2 <- FindVariableFeatures(so_round2, selection.method = "vst")
so_round2 <- ScaleData(object = so_round2)
so_round2 <- RunPCA(so_round2,
                    features = VariableFeatures(so_round2),
                    verbose = F,
                    ndims.print = 0)
ElbowPlot(so_round2, ndims = 40)
n_dim <- 30 
so_round2 <- FindNeighbors(so_round2, 
                           reduction = "pca",
                           dims = 1:n_dim, 
                           k.param = 20) 
resolutionsCalculated <- c(2, 3, 4)
so_round2 <- FindClusters(so_round2, resolution = resolutionsCalculated) 
so_round2_meta <- so_round2@meta.data
for(i in resolutionsCalculated){ # Want a resolution with a high number of clusters
  print(paste0("res:", i))
  print(length(unique(so_round2_meta[,paste0("RNA_snn_res." , i)])))
}
so_round2$seurat_clusters <- so_round2$RNA_snn_res.4
prop_df_sample <- data.frame(prop.table(table(so_round2@meta.data[,"scDblFinder.class"], 
                                              so_round2@meta.data[,"seurat_clusters"]), margin = 2))
prop_df_sample_doublets <- prop_df_sample[prop_df_sample$Var1 == "doublet",]
gg_stackbar_sample <- ggplot(prop_df_sample, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +ylab("Proportion") + xlab("New cluster") +
  theme_classic()
gg_stackbar_sample
remove <- as.character(prop_df_sample_doublets$Var2[prop_df_sample_doublets$Freq > threshold])

# FILTER round 2
so_round2_post <- subset(so_round2,
                         seurat_clusters %in% remove,
                         invert=T)

## Round 3
so_round3 <- so_round2_post
DefaultAssay(so_round3)
so_round3 <- NormalizeData(so_round3, normalization.method = "LogNormalize")
so_round3 <- FindVariableFeatures(so_round3, selection.method = "vst")
so_round3 <- ScaleData(object = so_round3)
so_round3 <- RunPCA(so_round3,
                    features = VariableFeatures(so_round3),
                    verbose = F,
                    ndims.print = 0)
ElbowPlot(so_round3, ndims = 40)
n_dim <- 30 
so_round3 <- FindNeighbors(so_round3, 
                           reduction = "pca",
                           dims = 1:n_dim, 
                           k.param = 20) 
resolutionsCalculated <- c(2, 3, 4)
so_round3 <- FindClusters(so_round3, resolution = resolutionsCalculated) 
so_round3_meta <- so_round3@meta.data
for(i in resolutionsCalculated){ # Want a resolution with a high number of clusters
  print(paste0("res:", i))
  print(length(unique(so_round3_meta[,paste0("RNA_snn_res." , i)])))
}
so_round3$seurat_clusters <- so_round3$RNA_snn_res.4
prop_df_sample <- data.frame(prop.table(table(so_round3@meta.data[,"scDblFinder.class"], 
                                              so_round3@meta.data[,"seurat_clusters"]), margin = 2))
prop_df_sample_doublets <- prop_df_sample[prop_df_sample$Var1 == "doublet",]
gg_stackbar_sample <- ggplot(prop_df_sample, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +ylab("Proportion") + xlab("New cluster") +
  theme_classic()
gg_stackbar_sample
remove <- as.character(prop_df_sample_doublets$Var2[prop_df_sample_doublets$Freq > threshold])

# FILTER round 3
so_round3_post <- subset(so_round3,
                         seurat_clusters %in% remove,
                         invert=T)
so_round3
so_round3_post

## Round 3
so_round4 <- so_round3_post
DefaultAssay(so_round4)
so_round4 <- NormalizeData(so_round4, normalization.method = "LogNormalize")
so_round4 <- FindVariableFeatures(so_round4, selection.method = "vst")
so_round4 <- ScaleData(object = so_round4)
so_round4 <- RunPCA(so_round4,
                    features = VariableFeatures(so_round4),
                    verbose = F,
                    ndims.print = 0)
ElbowPlot(so_round4, ndims = 40)
n_dim <- 30 
so_round4 <- FindNeighbors(so_round4, 
                           reduction = "pca",
                           dims = 1:n_dim, 
                           k.param = 20) 
resolutionsCalculated <- c(2, 3, 4)
so_round4 <- FindClusters(so_round4, resolution = resolutionsCalculated) 
so_round4_meta <- so_round4@meta.data
for(i in resolutionsCalculated){ # Want a resolution with a high number of clusters
  print(paste0("res:", i))
  print(length(unique(so_round4_meta[,paste0("RNA_snn_res." , i)])))
}
so_round4$seurat_clusters <- so_round4$RNA_snn_res.4
prop_df_sample <- data.frame(prop.table(table(so_round4@meta.data[,"scDblFinder.class"], 
                                              so_round4@meta.data[,"seurat_clusters"]), margin = 2))
prop_df_sample_doublets <- prop_df_sample[prop_df_sample$Var1 == "doublet",]
gg_stackbar_sample <- ggplot(prop_df_sample, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +ylab("Proportion") + xlab("New cluster") +
  theme_classic()
gg_stackbar_sample
remove <- as.character(prop_df_sample_doublets$Var2[prop_df_sample_doublets$Freq > threshold])

# FILTER round 4
so_round4_post <- subset(so_round4,
                         seurat_clusters %in% remove,
                         invert=T)
so_round4
so_round4_post

## Round 5
so_round5 <- so_round4_post
DefaultAssay(so_round5)
so_round5 <- NormalizeData(so_round5, normalization.method = "LogNormalize")
so_round5 <- FindVariableFeatures(so_round5, selection.method = "vst")
so_round5 <- ScaleData(object = so_round5)
so_round5 <- RunPCA(so_round5,
                    features = VariableFeatures(so_round5),
                    verbose = F,
                    ndims.print = 0)
ElbowPlot(so_round5, ndims = 40)
n_dim <- 30 
so_round5 <- FindNeighbors(so_round5, 
                           reduction = "pca",
                           dims = 1:n_dim, 
                           k.param = 20) 
resolutionsCalculated <- c( 4)
so_round5 <- FindClusters(so_round5, resolution = resolutionsCalculated) 
so_round5_meta <- so_round5@meta.data
for(i in resolutionsCalculated){ # Want a resolution with a high number of clusters
  print(paste0("res:", i))
  print(length(unique(so_round5_meta[,paste0("RNA_snn_res." , i)])))
}
so_round5$seurat_clusters <- so_round5$RNA_snn_res.4
prop_df_sample <- data.frame(prop.table(table(so_round5@meta.data[,"scDblFinder.class"], 
                                              so_round5@meta.data[,"seurat_clusters"]), margin = 2))
prop_df_sample_doublets <- prop_df_sample[prop_df_sample$Var1 == "doublet",]
gg_stackbar_sample <- ggplot(prop_df_sample, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +ylab("Proportion") + xlab("New cluster") +
  theme_classic()
gg_stackbar_sample
remove <- as.character(prop_df_sample_doublets$Var2[prop_df_sample_doublets$Freq > threshold])

# Ending after round 4 because no clusters > 50% doublets
#saveRDS(so_round4_post, 
#        "Data/Seurat_objects/AggrNoNorm_filtered_temp_round4_mergedblfind.RDS")


# Filter using cells that passed the final round:
sobj_metrics_filtered
so_round4_post
sobj_metrics_filtered_and_cluster <- subset(sobj_metrics_filtered,
                                            cells = colnames(so_round4_post))
sobj_metrics_filtered_and_cluster



################################################################################
####### -------        Filter out all remaining doublets         ------- #######

# Filter by singlet as well
sobj <- subset(sobj_metrics_filtered_and_cluster,
               scDblFinder.class == "singlet")

sobj_meta <- sobj@meta.data

sobj_filtered <- sobj

saveRDS(sobj_filtered, 
        "Data/Seurat_objects/AggrNoNorm_filtered.RDS")
sobj_filtered <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered.RDS")


################################################################################
####### -------            Post filter QC plot                   ------- #######

output_postfilter <- "Analysis/postFilter_QC"
dir.create(output_postfilter,
           recursive = T)

SeuratQC(so = sobj_filtered,
         category = "Sample",
         outdir = output_postfilter,
         violin_pt_size = 0.01)

gg_umi_filtered <- sobj_filtered %>%
  VlnPlot(
    features = c("nCount_RNA"),
    cols = rep("grey", 12),
    ncol     = 1,
    pt.size  = 0,
    group.by = "Sample"#,
  ) + # geom_hline(yintercept=umi_cut_vector) + 
  ylab("# of UMIs") + ggtitle("Number of UMIs") + xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text=element_text(size = gg_title_text),
        axis.text = element_text(size = gg_axis_text),
        legend.position="none")
gg_umi_filtered

gg_genes_filtered <- sobj_filtered %>%
  VlnPlot(
    features = c("nFeature_RNA"),
    cols = rep("grey", 12),
    ncol     = 1,
    pt.size  = 0,
    group.by = "Sample"#,
  ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            text=element_text(size = gg_title_text),
            axis.text = element_text(size = gg_axis_text),
            legend.position="none") + 
  ylab("# of genes") + ggtitle("Number of genes") + xlab("") # +
  # geom_hline(yintercept=gene_cut_vector)
gg_genes_filtered

gg_MT_filtered <- sobj_filtered %>%
  VlnPlot(
    features = c("percent_mito"),
    cols = rep("grey", 12),
    ncol     = 1,
    pt.size  = 0,
    group.by = "Sample"#,
  ) + # geom_hline(yintercept=mt_high_cut) + 
  ylab("UMI (%)") + ggtitle("UMI from mitochondria") + xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text=element_text(size = gg_title_text),
        axis.text = element_text(size = gg_axis_text),
        legend.position="none")
gg_MT_filtered

# Save an aggregated plot for ease of presentation
aggr_vio <- grid.arrange(gg_umi_filtered,
                         gg_genes_filtered,
                         gg_MT_filtered,
                         ncol = 3, nrow = 1)
ggsave(paste0(output_postfilter, "/Aggregated_violin_plots_post_filter.pdf"),
       aggr_vio,
       width = 16, height = 4)



# Save an aggregated plot WITH prefiltered
aggr_vio <- grid.arrange(prefilter_gg_umi,
                         prefilter_gg_genes,
                         prefilter_gg_MT,
                         gg_umi_filtered,
                         gg_genes_filtered,
                         gg_MT_filtered,
                         ncol = 3, nrow = 2)
ggsave(paste0(output_postfilter, "/Aggregated_violin_plots_post_and_pre_filter.pdf"),
       aggr_vio,
       width = 16, height = 8)



################################################################################
####### -------            Post filter tables                    ------- #######

# Here, sobj_filtered and sobj are equivalent 

# Dataset level
table_filter <- cbind.data.frame(Metric = c("Number of cells", 
                                               "Number of genes",
                                               "Median UMI per cell",
                                               "Median genes per cell"),
                                    Value = c(ncol(sobj),
                                              nrow(sobj),
                                              median(sobj$nCount_RNA),
                                              median(sobj$nFeature_RNA)))
write.table(table_filter,
            file = paste0(output_postfilter, "/dataset-level_metrics.txt"),
            sep = "\t", row.names = F, quote = F)

# Sample-level
table_sample <- cbind.data.frame(Sample = c("Number of cells", 
                                            "Number of genes",
                                            "Median UMI per cell",
                                            "Median genes per cell"))
for(samp in unique(sobj$Sample)){
  print(samp)
  
  sobj_samp <- subset(sobj,
                      Sample == samp)
  sobj_samp <- sobj_samp[rowSums(sobj_samp) > 0,]
  
  table_sample$samp <- c(ncol(sobj_samp),
                         nrow(sobj_samp),
                         median(sobj_samp$nCount_RNA),
                         median(sobj_samp$nFeature_RNA))
  colnames(table_sample)[ncol(table_sample)] <- samp
}

table_sample <- data.frame(t(table_sample))
write.table(table_sample,
            file = paste0(output_postfilter, "/sample-level_metrics.txt"),
            sep = "\t", row.names = T, col.names = F, quote = F)



gg_umi_vs_genes <- FeatureScatter(object = sobj_filtered, 
                                  feature1 = 'nCount_RNA', 
                                  feature2 = 'nFeature_RNA',
                                  group.by = "Sample",
                                  pt.size = 0.1)
gg_umi_vs_genes


gg_MT_vs_UMI <- FeatureScatter(object = sobj_filtered, 
                               feature2 = 'percent_mito', 
                               feature1 = 'nCount_RNA',
                               group.by = "Sample",
                               pt.size = 0.1)

gg_MT_vs_UMI 


