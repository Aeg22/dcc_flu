

################################################################################
####### -------                    Goals                         ------- #######

# Additional filtering of multipositive clusters
# Downstream analysis

################################################################################
####### -------                 Load R packages                  ------- #######

library(Seurat)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(msigdbr)
library(plyr)
library(dplyr)
library(clusterProfiler)
library(stringr)
library(openxlsx)
library(reshape2)
library(ggpubr)
library(enrichplot)
library(pheatmap)

source("functions/prop_hm.R")
source("functions/SeuratQC.R")
source("functions/gseaplot2_ag.R")
source("functions/fam_and_ora_functions.R")

################################################################################
####### -------                Update variables                  ------- #######

setwd("/Volumes/Partition 2/Core/degregori_continued_from_meher_scRNAseq_Oct2024/")

input_rds = "Data/Seurat_objects/AggrNoNorm_filtered.RDS"

output <- paste0("Analysis/")
dir.create(output, recursive = T)
output_full <- paste0(output, "Normalized/")
dir.create(output_full)



################################################################################
####### -------           Read in and subset data                ------- #######

# read in Seurat object and split
so <- readRDS(input_rds)

DefaultAssay(so) <- "RNA"

################################################################################
####### -------                 Normalize                        ------- #######

so <- NormalizeData(so, normalization.method = "LogNormalize")


################################################################################
####### -------               Add metrics                        ------- #######

s.genes <- str_to_title(cc.genes$s.genes)
g2m.genes <- str_to_title(cc.genes$g2m.genes)
so <- CellCycleScoring(so, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
so$CC.Difference <- so$S.Score - so$G2M.Score


################################################################################
####### -------               UMAP                               ------- #######

so <- FindVariableFeatures(so, selection.method = "vst")
so <- ScaleData(object = so, 
                vars.to.regress = c("nCount_RNA",
                                    "percent_mito",
                                    "CC.Difference"))
so <- RunPCA(so,
             features = VariableFeatures(so),
             verbose = F,
             ndims.print = 0)
ElbowPlot(so, ndims = 40)
n_dim <- 30 
so <- RunUMAP(so, 
              reduction = "pca", 
              dims = 1:n_dim)

gg_sample <- DimPlot(so, 
                     group.by = "Sample",
                     reduction = "umap")
gg_sample
ggsave(paste0(output_full, "UMAP_sample.pdf"),
       gg_sample)

gg_sample_split <- DimPlot(so, 
                           group.by = "Sample",
                           reduction = "umap",
                           split.by = "Sample")
gg_sample_split
ggsave(paste0(output_full, "UMAP_sample_split.pdf"),
       gg_sample_split, height = 3, width = 25)


################################################################################
####### -------               UMAP with QC metrics               ------- #######

output_umap_qc <- paste0(output_full, "/qc_umap/")
dir.create(output_umap_qc, recursive = T)

gg_MT <- FeaturePlot(so, 
                     features = "percent_mito",
                     reduction = "umap")
gg_MT 
ggsave(paste0(output_umap_qc, "UMAP_MT.pdf"),
       gg_MT)

gg_CCsscore <- FeaturePlot(so, 
                           features = "S.Score",
                           reduction = "umap")
gg_CCsscore 
ggsave(paste0(output_umap_qc, "UMAP_CCsscore.pdf"),
       gg_CCsscore)

gg_CCg2m <- FeaturePlot(so, 
                        features = "G2M.Score",
                        reduction = "umap")
gg_CCg2m
ggsave(paste0(output_umap_qc, "UMAP_CCg2Mscore.pdf"),
       gg_CCg2m)

gg_phase <- DimPlot(so, 
                    group.by = "Phase",
                    reduction = "umap")
gg_phase 
ggsave(paste0(output_umap_qc, "UMAP_CCPhase.pdf"),
       gg_phase)

gg_umi <- FeaturePlot(so, 
                      features = "nCount_RNA",
                      reduction = "umap")
gg_umi 
ggsave(paste0(output_umap_qc, "UMAP_UMIs.pdf"),
       gg_umi)


################################################################################
####### -------              Clustering                         ------- #######

# My intention here is broad cell type clustering

assay <- "RNA"

DefaultAssay(so) <- assay

so <- FindNeighbors(so, 
                    reduction = "pca",
                    dims = 1:n_dim, 
                    k.param = 20) 

resolutionsCalculated <- c(0.5, 1.0, 1.5, 2, 2.5, 3)
so <- FindClusters(so, resolution = resolutionsCalculated) 

for (i in resolutionsCalculated){
  
  # set resolution
  Idents(object = so) <- paste(assay, "_snn_res." , i, sep = "") 
  
  # general umap
  print(DimPlot(so,
                reduction = "umap",
                label = TRUE,
                label.size = 6) + ggtitle(paste(assay, "_snn_res." , i, sep = "")) + 
          NoLegend() + theme(plot.title = element_text(hjust = 0.5)))
}

# Choose a resolution
resolutionChosen <- 1.5
Idents(object = so) <- paste(assay, "_snn_res.", resolutionChosen, sep = "") 
so$seurat_clusters <- so@meta.data[,paste(assay, "_snn_res.", resolutionChosen, sep = "")] 
so$cluster_full <- so$seurat_clusters

umap_cluster <- DimPlot(so, reduction = "umap", label = T)
umap_cluster
ggsave(paste0(output_full, "UMAP_Cluster.PDF"), 
       plot = umap_cluster)

df_cluster_distribution <- data.frame(unclass(table(so$seurat_clusters)))
df_cluster_distribution <- cbind.data.frame(Cluster = rownames(df_cluster_distribution),
                                            Cells = df_cluster_distribution[,1])
write.table(df_cluster_distribution,
            file = paste0(output_full, "Cluster_distribution.txt"),
            sep = "\t", row.names = F, quote = F)
# Compare distribution of new clusters and samples
prop_df_sample <- data.frame(prop.table(table(so@meta.data[,"Sample"], 
                                              so@meta.data[,"seurat_clusters"]), margin = 2))
gg_stackbar_sample <- ggplot(prop_df_sample, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +ylab("Proportion") + xlab("New cluster") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
gg_stackbar_sample
ggsave(paste0(output_full, "/", "Sample_stacked_barplot.pdf"),
       gg_stackbar_sample)
write.table(prop_df_sample,
            file = paste0(output_full, "/", "Sample_stacked_barplot_table.txt"),
            sep = "\t", row.names = F, col.names = T, quote = F)

as.character(prop_df_sample$Var2[prop_df_sample$Freq > 0.5])

################################################################################
####### -------      Saving clustered RDS object temporarily     ------- #######

#saveRDS(so,
#        "Data/Seurat_objects/AggrNoNorm_filtered_cluster_temp.RDS")
#so <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_cluster_temp.RDS")



################################################################################
####### -------          Dotplot of expected markers             ------- #######

DotPlot(so,
        features = str_to_title(c(
          
          # Inflammatory monocyte (13)
          "CD14", "S100a8",
          
          # Complement macrophage (14)
          "ITGAM", "Adgre1", "CD68", "C1qa", "C1qb",
          
          # Mesothelial (23)
          "LRRN4", "UPK3B",
          
          # Aveolar type I (21)
          "Ager", "Gprc5a",
          
          # Aveolar type II (24, 27, 29)
          "Sftpc", "Sftpa1",
          
          # Dendritic (20)
          "Xcr1", "Clec9a",
          
          # Platelet (26)
          "Mmrn1", "ITGB3",
          
          # Fibroblast (22)
          "DCN", "COL1A1",
          
          # Smooth muscle (12, 17)
          "DES", "ITGA8", 
          
          # Endothelial (5, 18)
          "TMEM100", "Adgrf5",
          
          # CD4+CD8+ proliferating T (25)
          # none special, see above
          
          # CD4+CD8+ T (8)
          # none special, see above
          
          # General T
          "CD3E",
          
          # Common T
          "CD4", "CD8A",
          
          # CD4+IL7R+ T (0, 2, 4)
          "IL7R",
          
          # B (6,7,9)
          "Cd79a", "Ighd", 
          
          # ITGAM+ macrophage (15)
          # none special, see above
          
          # ITGAM- macrophage (19)
          # none special, see above
          
          # Basophil (28)
          "Mcpt8", "Fcer1a",
          
          # NK (1)
          "GZMB", "KLRC1",
          
          # CD4+ Effector T (10)
          "ICOS", "Ifng", 
          
          # CD8+ Proliferating Effector T (16)
          "Mki67", "Top2a",
          
          # CD8+ Effector T (3, 11)
          "CCL5", "ID2"
          
        )),
        group.by = "seurat_clusters") + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))



DotPlot(so,
        features = str_to_title(c(
          
          # Myocytes
          "Acta2",
          
          ## Additional fibroblast markers
          "Dcn","Col1a2", "Col5a1", "Loxl1", "Lum", "Fbln2", "Cd34","Pdgfra","S100a4",
          
          # Adventitial fibroblast
          "Col1a1","FBLN1","SCARA5","OGN","CCDC80","PI16","MFAP5",
          
          # Alveolar fibroblasts
          "NPNT","BMP5","MYH10","ANGPT1","FGFR4","MACF1","LIMCH1",
          
          # Myofibroblasts
          "DES","FABP3","TPM1","Cox7a1","Tnnt2",
          
          # Myeloid derived suppressor cell
          "Itgam", "S100a8","S100a9","Cd33",
          
          # Inflammatory monocyte
          #"Cd14",
          
          # Macrophage
          "Adgre1", "Cd68", "Fcer1g","C1qa","C1qb","Stat1","Irf1",
          
          #Dendritic
          "H2-Ab1","Clec9a",
          
          #Plasmacytoid dendritic
          "Siglech","Irf7", "CCL2","CCL17",
          
          #ILCs
          "Il1rl1","Cd200r1","Itgb3",
          
          #Mast
          "Ccl4","Ccl3","Hdc","Cyp11a1",
          
          #NK cells
          #"Gzmb",
          "Klrd1","Klrk1","Ncr1","Klrg1",
          
          #Mesothelial
          "Lrrn4", "Upk3b",
          
          #AT2
          "Sftpc","Sftpa1",
          
          #AT1
          "Ager","Gprc5a",
          
          # B (6,7,9)
          "Cd79a", "Ighd",
          
          # Megakaryocyte
          "Mmrn1","Ppbp", "Thbs1","Pf4","Itga2b",
          
          # Endothelial
          "Calcrl","Cdh2","Cdh5","Cldn5","Pecam1","Kdr","Eng",
          
          # Epithelial
          "Epcam","Cdh1",
          
          # Proliferating
          "Mki67","Top2a",
          
          #gCap
          "Epas1", "Kit","Gpihbp1","Ptprb",
          
          #aCap
          "Car4","Emp2","Ednrb",
          
          # known mural cell markers 
          "Mcam", "Tagln","Notch3", "Pdgfrb","CD13",
          
          # Additional T cell genes
          "Sell", # "Ccr7",
          
          # Additional NK
          "Nkg7", "Gzma",
          
          # cluster 16
          "Id2", "Cd5",
          
          # Fig 1D
          # https://pmc.ncbi.nlm.nih.gov/articles/PMC10187412/pdf/nihpp-rs2887159v1.pdf
          "Vwf", "Gja5", "Eln", "Ackr1", "N2rf2", # "Mmrn1",
          "Aplnr", # "Kit", 
          "Clic4", "Flot1", "Scn7A", "Vipr1",
          "Lingo2", "Fabp4", # "Ednrb",
          "Apln",
          
          # Cluster 28
          "Cxcr6", "Tmem64", "Zbtb16",
          
          # More cluster 28
          # https://www.sciencedirect.com/science/article/pii/S0896841121000615
          # Seem most closely activated double negative (cluster 1 but not the IL17)
          "Icos", "Il17a", "Il4", "Lta", "Serpinb6a", "Zbtb32", "Rora",
          # "Cxcr6",  Id2, Rora, Zbtb16 and Zbtb32
          # Same pub but cluster 0
          #"Xcl1", 
          "Ikzf2", "C1qbp", "Crtam", "Ncl",
          # Our double negatives could be a mix of these clusters because they are Ikzf2+
          
          # Cluster 5+16
          "Ifng", "Tnf",
          
          # cluster 24
          "Postn", "Bgn", "Myh11",
          "Ng2", "Cspg4", "Rgs5", "Cnn1",
          "Hhip", "Rgs10", "Kcnj8",
          
          
          # potential mural marker
          "Anpep",
          
          # Cluster 32
          "Ccl21a", "Flt4", "Gng11",
          
          # Cluster 37
          "Mgp",
          
          # Cluster 17,
          "Ccr1", "Ccr5", "Fcgr3", # This is mouse Cd16
          "Cd14", "XAF1",
          
          
          # T cell 
          "Cd3e","Cd4","Cd8a", "Il7r",   
          
          # T cell panel
          # https://www.nature.com/articles/s41467-021-23324-4
          #"Foxp3", "Ccr7", "Slamf6", "Ifngr1", "Pdcd1", "Ctla4",
          #"Xcl1", "Gzmb", "Cx3cr1", "Tox"
          # Re-ordered and more
          "Tcf7", "Ccr7", "Pdcd1", "Tnfrsf9", "Tnfrsf4",
          "Gzmb", "Gzmk",
          "Ctla4", "Lag3", "Tigit", "Havcr2", "Tox",
          "Ifngr1", "Fasl",
          "Cxcr5", "Slamf6",
          "Foxp3",
          "Xcl1"
          
        )),
        group.by = "seurat_clusters", cluster.idents = T) + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))

ggsave(filename = paste0(output_full, "dotplot_cluster_markers.pdf"),
       width = 22, height = 10)

so_one_rep <- subset(so_temp_lymph,
                     Sample == "Her2+IAV+anti-CD4_dpi15_rep2")
table(so_one_rep$seurat_clusters)

so_temp_lymph <- subset(so,
                        seurat_clusters %in% c(0,
                                               2,
                                               3,
                                               4,
                                               5,
                                               6,
                                               8,
                                               11,
                                               15,
                                               16,
                                               18,
                                               19,
                                               21,
                                               26,
                                               28,
                                               29,
                                               35))
DotPlot(so_temp_lymph,
        features = str_to_title(c(
          
          #ILCs
          "Il1rl1","Cd200r1","Itgb3",
          
          
          # Proliferating
          "Mki67","Top2a",
          
          # Additional T cell genes
          "Sell", # "Ccr7",
          
          # Additional NK
          "Nkg7", "Gzma",
          
          # Cluster 28
          "Cxcr6", "Tmem64", "Zbtb16",
          
          # More cluster 28
          # https://www.sciencedirect.com/science/article/pii/S0896841121000615
          # Seem most closely activated double negative (cluster 1 but not the IL17)
          "Icos", "Il17a", "Il4", "Lta", "Serpinb6a", "Zbtb32", "Rora",
          # "Cxcr6",  Id2, Rora, Zbtb16 and Zbtb32
          # Same pub but cluster 0
          #"Xcl1", 
          "Ikzf2", "C1qbp", "Crtam", "Ncl",
          # Our double negatives could be a mix of these clusters because they are Ikzf2+
          
          # Cluster 5+16
          "Ifng", "Tnf",
          
          # T cell 
          "Cd3e","Cd4","Cd8a", "Il7r",   
          
          "Lef1",
          
          # T cell panel
          # https://www.nature.com/articles/s41467-021-23324-4
          #"Foxp3", "Ccr7", "Slamf6", "Ifngr1", "Pdcd1", "Ctla4",
          #"Xcl1", "Gzmb", "Cx3cr1", "Tox"
          # Re-ordered and more
          "Tcf7", "Ccr7", "Pdcd1", "Tnfrsf9", "Tnfrsf4",
          "Gzmb", "Gzmk",
          "Ctla4", "Lag3", "Tigit", "Havcr2", "Tox",
          "Ifngr1", "Fasl",
          "Cxcr5", "Slamf6",
          "Foxp3",
          "Xcl1",
          
          # CD4 states
          # https://elifesciences.org/articles/76339
          "Lyc62", "Izumo1r", "Ifit1",
          
          "Eomes",
          "Id2", "Klrd1",
          
          "Il18r1", "Cd5", "Il2rb"
          
        )),
        group.by = "seurat_clusters", cluster.idents = T) + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))

ggsave(filename = paste0(output_full, "dotplot_cluster_markers_temp_lymph.pdf"),
       width = 12, height = 6)



so_temp_myeloid <- subset(so,
                          seurat_clusters %in% c(14,
                                                 17,
                                                 22,
                                                 23,
                                                 25))
DotPlot(so_temp_myeloid,
        features = str_to_title(c(
          
          # Myeloid derived suppressor cell
          #"Itgam", 
          "S100a8","S100a9","Cd33",
          
          # Inflammatory monocyte
          #"Cd14",
          
          # Macrophage
          #"Adgre1", "Cd68", 
          "Fcer1g","C1qa","C1qb","Stat1","Irf1",
          
          #Dendritic
          "H2-Ab1","Clec9a",
          
          #Plasmacytoid dendritic
          "Siglech","Irf7", "CCL2","CCL17",
          
          #Mast
          "Ccl4","Ccl3","Hdc","Cyp11a1",
          
          # Current
          "Adgre1", "Cd68", "Itgam", 
          "Cd14", "Fcgr3", "Ccr1",
          
          # M1/M2
          # https://pmc.ncbi.nlm.nih.gov/articles/PMC6542613/
          "Ptprc", "Lyz2",
          # M1:
          "Nos2", "Cd86", "Tnf",  "CXCL9", "CXCL10", "CXCL11", "CXCL16", "IL15", "GBP1",
          "GBP5", "OAS3", "CD38", "Il12b", "Oasl1", "Ptges", "Saa3", "H2-Aa", "Irf8",
          "Psme1", 
          # M2: 
          "Chil3", "Arg1", "Retnla", "Mgl2", "PPARG", "VWF", "PTGS1", "Mrc1",
          "Dut", "Tmem37", "Cd79a", "Ebf1", "Rela",
          
          
          # M1/M2
          # https://www.nature.com/articles/s41467-023-40156-6
          # M0/M1:
          "SELENOP", "EREG", "CD300E", "VCAN", "FCN1", "INHBA",
          # M2:
          "CD209", "CD163L1", "FUCA1", "FOLR2", "A2M",
          # MRC1,
          
          
          # https://www.nature.com/articles/s41467-018-07387-4
          # Basically, loss of Itgam (Cd11b) promotes protumor processes
          #  immune suppression and angiogenesis, such as Arg1, Tgfb, Il10, Il6, 
          # and Pdgfb, and significantly lower expression of genes associated with 
          # immune stimulation, such as Ifng, Nos2, and Tnfa
          
          "Tgfb1", "Il10", "Il6", "Pdgfb",
          "Ifng"
          
          
        )),
        group.by = "seurat_clusters", cluster.idents = T) + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))

ggsave(filename = paste0(output_full, "dotplot_cluster_markers_temp_myeloid.pdf"),
       width = 12, height = 6)

DotPlot(so_temp_myeloid,
        features = list(
          
          # General myeloid
          Mye = c("Ptprc", "Lyz2"),
          
          # Monocyte
          Mono = c("Cd14", "Fcgr3", "Ccr1"),
          
          # Macrophage
          Mac = c("Cd68", "Adgre1"),
          
          # Select M1 markers
          M1 = str_to_title(c("Nos2", "Cd86", "CXCL9", "CXCL10", "CXCL16", "GBP5", "OAS3", 
                              "CD38", "Oasl1",
                              "CD300E", "VCAN")),
          
          # Select M2 markers
          M2 = str_to_title(c("Chil3", "PPARG",  "Mrc1", "FUCA1")),
          
          # Higher in 22,25 compared to 14
          Higher_22_25 = c("Ifitm6", "Ace", "Nr4a1", "Klf2") , 
          # All have anti-tumor and/or anti-inflammation effects
          
          # Higher in 22,25 
          # KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION
          Higher_22_24_migration = c("Itgal",  "Itga4" , "Rassf5" ,"Ncf2",   "Rap1a",
                                     "Itgb2" , "Msn"  ,  "Plcg2" , "Rac2"  , "Mapk14",
                                     "Rap1b",  "Gnai2",  "Rhoa",   "Cybb" ),
          
          
          # Higher in 14 compared to 22,25
          Higher_14 = c("C1qa","C1qb", "Ccr5"),
          # Seems more inflammatory
          
          # Higher in 14
          
          Higher_14_interferon = c("H2-Aa" , "Cfb" , "Cd74" ,"Fgl2" ,  "Lgals3bp",
                                   "Isg15", "Il18bp" ,  "Fcgr1", "Stat1",
                                   "Bst2" , "Rtp4" , "Irf7" , "Cd274", "Nampt",
                                   "Rnf213" ,  "Il4ra" )
          
        ),
        group.by = "seurat_clusters", cluster.idents = T) + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))

ggsave(filename = paste0(output_full, "dotplot_cluster_markers_temp_myeloid_trim.pdf"),
       width = 15, height = 6)


################################################################################
####### -------  Cluster-level QC (whole dataset and by sample)  ------- #######

output_qc_cluster <- paste0(output_full, "/qc_by_cluster/")
dir.create(output_qc_cluster)

# Whole dataset
SeuratQC(so = so,
         category = "seurat_clusters",
         outdir = output_qc_cluster,
         violin_pt_size = 0)

# By sample
for(sample in unique(so$Sample)){
  print(sample)
  temp_dir <- paste0(output_qc_cluster, sample, "/")
  dir.create(temp_dir)
  so_temp <- subset(so, Sample == sample)
  SeuratQC(so = so_temp,
           category = "seurat_clusters",
           outdir = temp_dir,
           violin_pt_size = 0)
}

################################################################################
####### -------       Help with cell type identity               ------- #######

####### -------       ORA-based method

Idents(so) <- so$seurat_clusters

so.markers_orig <- FindAllMarkers(so, 
                                  only.pos = TRUE, 
                                  min.pct = 0.5, 
                                  logfc.threshold = 0.5)
so.markers <- so.markers_orig[so.markers_orig$p_val_adj < 0.05,] # Only keep sig genes
table(so.markers$cluster) 

output_ORA <- paste0(output_full, "/ORA_clusters/")
dir.create(output_ORA)

# Writing the sig tables
write.table(so.markers,
            file=paste0(output_ORA,
                        "Clusters_FindAllMarkers_sigOnly.txt"),
            sep = "\t",
            row.names = F)
so.markers <- read.table(file=paste0(output_ORA,
                                     "Clusters_FindAllMarkers_sigOnly.txt"),
                         sep = "\t", header = T, stringsAsFactors = F)

# Subsetting to top genes
so.markers_ORA <- so.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC) # 

all = msigdbr(species = "Mus musculus") 

specific_output <- "ORA_C8"

gene_sets <- all[all$gs_cat %in% "C8",] # C8 only
length(unique(gene_sets$gs_name))
gene_sets <- gene_sets %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

# Cycle through each cluster and perform ORA
clusters <- unique(so.markers_ORA$cluster)
background <- rownames(so) # All gene names in dataset
output_ORA_collection <- paste0(output_ORA, specific_output, "/")
dir.create(output_ORA_collection,
           recursive = T)
for(cluster in clusters){
  print(cluster)
  
  genes <- so.markers_ORA$gene[so.markers_ORA$cluster %in% cluster]
  
  ora_enricher <- enricher(genes, TERM2GENE=gene_sets, pAdjustMethod = "fdr", 
                           minGSSize = 5, maxGSSize = 1000, 
                           qvalueCutoff = 1, pvalueCutoff =1,
                           universe = background) 
  ora_result <- data.frame(ora_enricher)
  
  write.table(ora_result,
              file = paste0(output_ORA_collection,
                            cluster, ".txt"),
              row.names = F, sep = "\t", quote = F)
}


# Testing out PangloaDB (plus with technical gene sets)
panglao_gene_set_plus <- read.table(file = paste0("helper_data/PanglaoDB_plus_mouse_GeneSets.csv"),
                                    sep = ",", header = T, stringsAsFactors = F)
panglao_gene_set_plus$Collection_Shortname <- "PanglaoDB_plus"
panglao_gene_set_plus$Gene_Set_Shortname <- panglao_gene_set_plus$gs_name

# Get gene sets in clusterprofiler format
gene_sets_to_analyze = panglao_gene_set_plus[,c("Gene_Set_Shortname", "gene_symbol")]

specific_output <- "PanglaoDB_plus"
output_ORA_collection <- paste0(output_ORA, specific_output, "/")
dir.create(output_ORA_collection,
           recursive = T)
for(cluster in clusters){
  print(cluster)
  
  genes <- so.markers_ORA$gene[so.markers_ORA$cluster %in% cluster]
  
  ora_enricher <- enricher(genes, TERM2GENE=gene_sets_to_analyze, pAdjustMethod = "fdr", 
                           minGSSize = 5, maxGSSize = 1000, 
                           qvalueCutoff = 1, pvalueCutoff =1,
                           universe = background) 
  ora_result <- data.frame(ora_enricher)
  
  write.table(ora_result,
              file = paste0(output_ORA_collection,
                            cluster, ".txt"),
              row.names = F, sep = "\t", quote = F)
}


################################################################################
####### -------    Name cells (prior to categorizing T cells)    ------- #######

annotation <- read.table(file = "helper_data/maincluster_cell_types.txt",
                         sep = "\t", header = T, stringsAsFactors = F)
cluster_type <- as.list(annotation$Cell.Type)
names(cluster_type) <- annotation$Cluster
cluster_broad <- as.list(annotation$Broad.Cell.Type)
names(cluster_broad) <- annotation$Cluster

# Discrete cell type
so$Cell.Type <- "cell"
for(cluster in unique(names(cluster_type))){
  so$Cell.Type[so$seurat_clusters == cluster] <- cluster_type[[cluster]]
}
# Broad cell type
so$Broad.Cell.Type <- "cell"
for(cluster in unique(names(cluster_broad))){
  so$Broad.Cell.Type[so$seurat_clusters == cluster] <- cluster_broad[[cluster]]
}


################################################################################
## -- All clusters including those to be removed - Cell type UMAP and marker dot plot -- ##

umap_celltype_all <- DimPlot(so, 
                             reduction = "umap", 
                             label = T,
                             group.by = "Cell.Type",
                             repel = T,
                             label.size = 3)
umap_celltype_all
ggsave(paste0(output_full, "UMAP_celltype_all_clusters.PDF"), 
       plot = umap_celltype_all, width = 14, height = 8)

# Labeling this dot plot with cluster.Cell.Type
so$cluster_Cell.Type <- as.character(paste0(so$seurat_clusters, ".",
                                            so$Cell.Type))
so$cluster_Cell.Type <- factor(so$cluster_Cell.Type,
                               levels = c("30.AT1",
                                          "36.REMOVE_AT2 + collagen",
                                          "1.B",
                                          "7.B",
                                          "10.B",
                                          "9.Capillary endothelial",
                                          "12.Capillary endothelial",
                                          "20.Car4+ endothelial",
                                          "32.Lymphatic endothelial",
                                          "31.Venous endothelial",
                                          "27.Adventitial fibroblast",
                                          "13.Alveolar fibroblast",
                                          "37.Inflammatory adventitial fibroblast",
                                          "14.M1 macrophage",
                                          "23.M2 macrophage",
                                          "22.Migratory macrophage",
                                          "25.Migratory macrophage",
                                          "38.REMOVE_Macrophage-Myofibroblast",
                                          "33.Mesothelial",
                                          "17.Monocyte",
                                          "6.Natural Killer",
                                          "11.Natural Killer",
                                          "24.Pericyte",
                                          "34.Smooth muscle",
                                          "16.T",
                                          "5.T",
                                          "0.T",
                                          "2.T",
                                          "19.T",
                                          "21.T",
                                          "3.T",
                                          "15.T",
                                          "18.T",
                                          "8.T",
                                          "28.T",
                                          "4.REMOVE_CD4+CD8+",
                                          "26.REMOVE_CD4+CD8+",
                                          "29.REMOVE_CD4+CD8+",
                                          "35.REMOVE_CD4+CD8+"))
DotPlot(so,
        features = str_to_title(c(
          
          #AT1
          "Ager","Gprc5a",
          
          #AT2
          "Sftpc","Sftpa1",
          
          # B (6,7,9)
          "Cd79a", "Ighd",
          
          # General endothelial
          "Cdh5","Cldn5",
          
          # Capillary endothelial
          "Gpihbp1", "Ptprb", "Kit",
          
          # Car4+ endothelial
          "Car4", "Cd34",
          
          # Lymphatic endothelial
          "Mmrn1", "Ccl21a", "Flt4",
          
          # Venous endothelial
          "Vwf", "Ackr1",
          
          # General fibroblast (specifics are a combo of the following markers)
          "Col1a1", 
          
          # Adventitial fibroblast
          "Dcn", "Pdgfra", "Pi16", "Fbln1",
          
          # Alveolar fibroblast
          "Npnt", "Fgfr4",
          
          # Inflammatory adventitial fibroblast
          "C3", "C1qa",
          
          # Myeloid
          "Ptprc",
          
          # Monocyte
          "Cd14", "Fcgr3", "Ccr1",
          
          # General macrophage
          "Adgre1", "Cd68",
          
          # M1 macrophage
          "Nos2", "Cd86", "CXCL9",
          
          # M2 macrophage
          "Pparg",  "Mrc1", "Fuca1",
          
          # Migratory macrophage
          "Itgal",  "Itga4" ,"Ncf2", "Rhoa",
          
          # Mesothelial
          "Lrrn4", "Upk3b",
          
          # Natural Killer
          #"Nkg7", 
          "Gzma", "Gzmb", "Klrd1", "Klrk1", 
          
          # Pericyte
          "Myh11", "Pdgfrb", "Cspg4",
          
          # Smooth muscle
          "Acta2", "Hhip",
          
          # General T
          "Cd3e", "Cd4", "Cd8a",
          
          # CD4+ effector
          "Nkg7", "Id2", "Cxcr6",
          
          # CD4+ effector intermediate
          
          # CD4+ memory or naive
          "Il7r", "Ccr7", "Sell",
          
          # Proliferating T
          "Mki67", "Top2a",
          
          # CD8+ Effector T
          
          # CD8+ memory or naive
          
          # Double negative
          "Zbtb16", "Rora", "Ikzf2",
          
          # Inhibitory/exhausted T markers
          "Pdcd1", "Ctla4", "Lag3", "Tigit", "Havcr2", "Tox"
          
          
        )),
        group.by = "cluster_Cell.Type", cluster.idents = F) + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))

ggsave(filename = paste0(output_full, "dotplot_all_clusters_celltype_markers.pdf"),
       width = 20, height = 10)


################################################################################
####### -------             Subclustering T cells                ------- #######

so_t <- subset(so,
               Broad.Cell.Type == "T")

so_t <- FindVariableFeatures(so_t, selection.method = "vst")
so_t <- ScaleData(object = so_t, 
                  vars.to.regress = c("nCount_RNA",
                                      "percent_mito",
                                      "CC.Difference"))
so_t <- RunPCA(so_t,
               features = VariableFeatures(so_t),
               verbose = F,
               ndims.print = 0)
ElbowPlot(so_t, ndims = 40)
n_dim_t <- 20 
so_t <- RunUMAP(so_t, 
                reduction = "pca", 
                dims = 1:n_dim_t)

## Cluster
assay <- "RNA"
DefaultAssay(so_t) <- assay
so_t <- FindNeighbors(so_t, 
                      reduction = "pca",
                      dims = 1:n_dim_t, 
                      k.param = 20) 
resolutionsCalculated <- c(0.5, 1.5, 2, 3)
so_t <- FindClusters(so_t, resolution = resolutionsCalculated) 
for (i in resolutionsCalculated){
  # set resolution
  Idents(object = so_t) <- paste(assay, "_snn_res." , i, sep = "") 
  # general umap
  print(DimPlot(so_t,
                reduction = "umap",
                label = TRUE,
                label.size = 6) + ggtitle(paste(assay, "_snn_res." , i, sep = "")) + 
          NoLegend() + theme(plot.title = element_text(hjust = 0.5)))
}
# Choose a resolution
resolutionChosen <- 1.5
Idents(object = so_t) <- paste(assay, "_snn_res.", resolutionChosen, sep = "") 
so_t$seurat_clusters <- so_t@meta.data[,paste(assay, "_snn_res.", resolutionChosen, sep = "")] 

t_dist <- data.frame(unclass(table(so_t$seurat_clusters, so_t$Sample)))
table(so_t$Sample[so_t$seurat_clusters == "16"])

FeaturePlot(so_t, 
            features = "Foxp3",
            reduction = "umap")

FeaturePlot(so_t, 
            features = "Cd8a",
            reduction = "umap")

FeaturePlot(so_t, 
            features = "Cd4",
            reduction = "umap")

DotPlot(so_t,
        features = str_to_title(c(
          
          # Natural Killer
          #"Nkg7", 
          "Gzma", "Gzmb", "Klrd1", "Klrk1", 
          
          # General T
          "Cd3e", "Cd4", "Cd8a",
          
          # CD4+ effector
          "Nkg7", "Id2", "Cxcr6",
          
          # CD4+ effector intermediate
          
          # CD4+ memory or naive
          "Il7r", "Ccr7", "Sell",
          
          # Proliferating T
          "Mki67", "Top2a",
          
          # CD8+ Effector T
          
          # CD8+ memory or naive
          
          # Double negative
          "Zbtb16", "Rora", "Ikzf2",
          
          # Inhibitory/exhausted T markers
          "Pdcd1", "Lag3", "Tigit", "Havcr2", "Tox",
          
          
          # Treg
          str_to_title(c("Foxp3", "Il2ra", # "Cd25", 
                         "Tnfrsf18", "CD127", "CTLA4", "CD62L")),
          
          # Nuocytes are typically  IL1RL1+ Icos+ and negative for CD4, CD8, CD3
          # cluster 21 nuocyte markers
          "Il1rl1", "Icos", "Il17rb", "Gata3",
          "Cd81", "Arg1", "Il13", "Il5",
          
          "Il17a", "Il17f", "Il17re", "Il17rc", "Il17ra", "Il17c", "Il17rd",  "Il17d",  "Il17b" 
          
          
        )),
        group.by = "seurat_clusters", cluster.idents = T) + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
ggsave(filename = paste0(output_full, "t_subcluster_cluster_dotplot.pdf"))

t_tab <- data.frame(unclass(table(as.character(so_t$Cell.Type),
                                  so_t$seurat_clusters)))

DimPlot(so_t, 
        group.by = "seurat_clusters",
        reduction = "umap", label = T)
ggsave(filename = paste0(output_full, "t_subcluster_cluster_umap.pdf"))


dge_21 <- FindMarkers(so_t,
                      ident.1 = "21",
                      only.pos = T)
dge_21 <- cbind.data.frame(Gene = rownames(dge_21), dge_21)
write.table(dge_21, paste0(output_full, "T_21_vs_rest_of_T.txt"),
            sep = "\t", row.names = F, quote = F)

t_barcodes_21 <- colnames(so_t)[so_t$seurat_clusters == "21"]
so_temp <- so
so_temp$t21 <- "else"
so_temp$t21[colnames(so_temp) %in% t_barcodes_21] <- "t21"
dge_21 <- FindMarkers(so_temp,
                      ident.1 = "t21",
                      only.pos = T,
                      group.by = "t21")
dge_21 <- cbind.data.frame(Gene = rownames(dge_21), dge_21)
write.table(dge_21, paste0(output_full, "T_21_vs_rest_of_dataset.txt"),
            sep = "\t", row.names = F, quote = F)
dge_21 <- read.table(paste0(output_full, "T_21_vs_rest_of_dataset.txt"),
                     sep = "\t", header = T, stringsAsFactors = F)

dge_21 <- dge_21[order(dge_21$avg_log2FC, decreasing = T),]
dge_21 <- dge_21[dge_21$p_val_adj < 0.05,]
t21_genes <- dge_21$Gene[1:200]
all = msigdbr(species = "Mus musculus") 
background <- rownames(so) # All gene names in dataset

# C8 only
gene_sets <- all[all$gs_cat %in% "C8",] 
gene_sets <- gene_sets %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
ora_enricher <- enricher(t21_genes, TERM2GENE=gene_sets, pAdjustMethod = "fdr", 
                         minGSSize = 5, maxGSSize = 1000, 
                         qvalueCutoff = 1, pvalueCutoff =1,
                         universe = background) 
ora_result <- data.frame(ora_enricher)
write.table(ora_result,
            file = paste0(output_full, "ORA_C8_T21_vs_rest_of_dataset.txt"),
            row.names = F, sep = "\t", quote = F)

# PangloaDB (plus with technical gene sets)
panglao_gene_set_plus <- read.table(file = paste0("helper_data/PanglaoDB_plus_mouse_GeneSets.csv"),
                                    sep = ",", header = T, stringsAsFactors = F)
panglao_gene_set_plus$Collection_Shortname <- "PanglaoDB_plus"
panglao_gene_set_plus$Gene_Set_Shortname <- panglao_gene_set_plus$gs_name
gene_sets_to_analyze = panglao_gene_set_plus[,c("Gene_Set_Shortname", "gene_symbol")]
ora_enricher <- enricher(t21_genes, TERM2GENE=gene_sets_to_analyze, pAdjustMethod = "fdr", 
                         minGSSize = 5, maxGSSize = 1000, 
                         qvalueCutoff = 1, pvalueCutoff =1,
                         universe = background) 
ora_result <- data.frame(ora_enricher)
write.table(ora_result,
            file = paste0(output_full, "ORA_PangloaDb_T21_vs_rest_of_dataset.txt"),
            row.names = F, sep = "\t", quote = F)

# Add T labels
t_annotation <- read.table(file = "helper_data/t_cell_types.txt",
                         sep = "\t", header = T, stringsAsFactors = F)
Tcluster_type <- as.list(t_annotation$Cell.Type)
names(Tcluster_type) <- t_annotation$Cluster
Tcluster_broad <- as.list(t_annotation$Broad.Cell.Type)
names(Tcluster_broad) <- t_annotation$Cluster

# Discrete cell type
so_t$Cell.Type <- "cell"
for(cluster in unique(names(Tcluster_type))){
  so_t$Cell.Type[so_t$seurat_clusters == cluster] <- Tcluster_type[[cluster]]
}
# Broad cell type
so_t$Broad.Cell.Type <- "cell"
for(cluster in unique(names(Tcluster_broad))){
  so_t$Broad.Cell.Type[so_t$seurat_clusters == cluster] <- Tcluster_broad[[cluster]]
}

table(so_t$Cell.Type)
table(so_t$Broad.Cell.Type)

so_t$cluster_Cell.Type <- paste0(so_t$Cell.Type, ".", so_t$seurat_clusters)
table(so_t$cluster_Cell.Type)
so_t$cluster_Cell.Type <- factor(so_t$cluster_Cell.Type)

# proportions of these t cells (and nuocytes)
dir.create(paste0(output_full, "proportions_temp_t_by_sample/"))
prop_hm(so = so_t,
        group1 = "Sample",
        group2 = "cluster_Cell.Type",
        order_by = NULL,
        outdir = paste0(output_full, "proportions_temp_t_by_sample/"),
        breaks = seq(0, 1, by = 0.01), # Sets max color to 1
        height = 8,
        width = 10,
        decimals = 4
)

# Save T cell subtype
#saveRDS(so_t, "Data/Seurat_objects/AggrNoNorm_filtered_T_labeled.RDS")
#so_t <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_T_labeled.RDS")



################################################################################
####### -------      Merge T labels with the full dataset        ------- #######

# If I need T clusters or UMAP projections, I can get them from the T RDS object above

# Have two meta data file
# If one barcode is in the so_t object, use that broad and cell type label
# also, replace seurat_object with example t1, t2, t3 to denote it is the cluster from t subset

subclustered_cells <- colnames(so_t)

so_meta_pre <- so@meta.data

# Precompute indices for main and subcluster datasets
main_cells <- colnames(so)
subcluster_cells <- colnames(so_t)

# Use match() to find the indices in a vectorized way
n_main_indices <- match(subclustered_cells, main_cells)
n_subcluster_indices <- match(subclustered_cells, subcluster_cells)

# Check if any cells from subclustered_cells are not found in so
missing_cells <- which(is.na(n_main_indices))
if (length(missing_cells) > 0) {
  cat("Warning: The following cells from subclustered_cells are not found in 'so':\n")
  print(subclustered_cells[missing_cells])
}

# Now update only those cells that have matches in 'so'
so$old_seurat_clusters <- so$seurat_clusters
so$Cell.Type <- as.character(so$Cell.Type)
so$seurat_clusters <- as.character(so$seurat_clusters)
so$Broad.Cell.Type <- as.character(so$Broad.Cell.Type)
valid_cells <- which(!is.na(n_main_indices))  # Get indices of valid matches
if (length(valid_cells) > 0) {
  # Update the 'so' object for valid cells
  so$seurat_clusters[n_main_indices[valid_cells]] <- paste0("t", so_t$seurat_clusters[n_subcluster_indices[valid_cells]])
  so$Cell.Type[n_main_indices[valid_cells]] <- so_t$Cell.Type[n_subcluster_indices[valid_cells]]
  so$Broad.Cell.Type[n_main_indices[valid_cells]] <- so_t$Broad.Cell.Type[n_subcluster_indices[valid_cells]]
} else {
  cat("No valid cells found in 'so' to update.\n")
}

so_meta_post <- so@meta.data

table(so$seurat_clusters)
table(so$Broad.Cell.Type, so$seurat_clusters)

clusters_df <- data.frame(unclass(table(so$old_seurat_clusters, so$seurat_clusters)))

table(so$Cell.Type)

table(so$Broad.Cell.Type)


################################################################################
####### -------          Removing troublesome clusters           ------- #######

# Include removal of 1 T cluster as well

so
so <- subset(so,
             Broad.Cell.Type == "REMOVE",
             invert = T)
so



################################################################################
####### -------        Cell type UMAP and marker dot plot        ------- #######

# https://github.com/satijalab/seurat/issues/2264
getDiscretePalette <- function(ident.used = GetClusteringRuns()[1],
                               obj = combined.obj
                               , palette.used = c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[1]
                               , show.colors = F) {
  n.clusters <-  nrow(unique(obj[[ident.used]]))
  colz <- DiscretePalette(n = n.clusters, palette = palette.used)
  if (show.colors) Color_Check(colz)
  return(colz)
}

umap_celltype <- DimPlot(so, 
                         cols = getDiscretePalette(obj = so, ident.used = "Cell.Type"),
                         reduction = "umap", 
                         label = T,
                         group.by = "Cell.Type",
                         repel = T,
                         label.size = 3)
umap_celltype
ggsave(paste0(output_full, "UMAP_celltype.PDF"), 
       plot = umap_celltype,
       width = 13.5, height = 8)

# Sample plot by broad group
so$Broad.Group <- so$Group
umap_celltype_facet <- DimPlot(so, 
                               cols = getDiscretePalette(obj = so, ident.used = "Cell.Type"),
                               reduction = "umap", 
                               label = T,
                               group.by = "Cell.Type",
                               repel = T,
                               label.size = 1.5,
                               split.by = "Broad.Group",
                               ncol = 4)
ggsave(paste0(output_full, "UMAP_celltype_facet_sample.PDF"), 
       plot = umap_celltype_facet,
       width = 16, height = 8)

# Select samples by day and replicate
unique(so$Sample)
so_sub <- subset(so,
                 Broad.Group %in% c("Her2+IAV+anti-CD4_dpi9",
                                    "Her2+IAV+anti-CD4_dpi15",
                                    "Her2+IAV_dpi9",
                                    "Her2+IAV_dpi15"))
so_sub$replicate <- "Replicate 1"
so_sub$replicate[grep("rep2", so_sub$Sample)] <- "Replicate 2"
umap_replicate_facet <- DimPlot(so_sub, 
                                reduction = "umap", 
                                label = F,
                                group.by = "replicate",
                                split.by = "Broad.Group",
                                shuffle=T)
ggsave(paste0(output_full, "UMAP_celltype_facet_replicate.PDF"), 
       plot = umap_replicate_facet,
       width = 12, height = 4)

# Labeling this dot plot with cluster.Cell.Type
so$cluster_Cell.Type <- as.character(paste0(so$seurat_clusters, ".",
                                            so$Cell.Type))
unique(so$cluster_Cell.Type)
unique(so$cluster_Cell.Type)[grep("^t", unique(so$cluster_Cell.Type))]
so$cluster_Cell.Type <- factor(so$cluster_Cell.Type,
                               levels = c("30.AT1",
                                          "1.B",
                                          "7.B",
                                          "10.B",
                                          "9.Capillary endothelial",
                                          "12.Capillary endothelial",
                                          "20.Car4+ endothelial",
                                          "32.Lymphatic endothelial",
                                          "31.Venous endothelial",
                                          "27.Adventitial fibroblast",
                                          "13.Alveolar fibroblast",
                                          "37.Inflammatory adventitial fibroblast",
                                          "24.Pericyte",
                                          "34.Smooth muscle",
                                          "17.Monocyte",
                                          "14.M1 macrophage",
                                          "23.M2 macrophage",
                                          "22.Migratory macrophage",
                                          "25.Migratory macrophage",
                                          "33.Mesothelial",
                                          "6.Natural Killer",
                                          "11.Natural Killer",
                                          "t3.CD4+ effector T",
                                          "t5.CD4+ effector T",
                                          "t0.CD4+ memory T",
                                          "t1.CD4+ memory T",
                                          "t7.CD4+ memory T",
                                          "t14.CD4+ memory T",
                                          "t15.Proliferating T",
                                          "t20.Proliferating T",
                                          "t4.CD8+ effector T",
                                          "t6.CD8+ effector T",
                                          "t8.CD8+ effector T",
                                          "t9.CD8+ effector T",
                                          "t12.CD8+ effector T",
                                          "t2.CD8+ memory T",
                                          "t10.CD8+ memory T",
                                          "t13.CD8+ memory T",
                                          "t11.Double negative T",
                                          "t17.Double negative T", 
                                          "t19.Double negative T",
                                          "t18.Regulatory T",
                                          "t21.Nuocyte"                                  
                               ))
DotPlot(so,
        features = str_to_title(c(
          
          #AT1
          "Ager","Gprc5a",
          
          #AT2
          "Sftpc","Sftpa1",
          
          # B (6,7,9)
          "Cd79a", "Ighd",
          
          # General endothelial
          "Cdh5","Cldn5",
          
          # Capillary endothelial
          "Gpihbp1", "Ptprb", "Kit",
          
          # Car4+ endothelial
          "Car4", "Cd34",
          
          # Lymphatic endothelial
          "Mmrn1", "Ccl21a", "Flt4",
          
          # Venous endothelial
          "Vwf", "Ackr1",
          
          # General fibroblast (specifics are a combo of the following markers)
          "Col1a1", 
          
          # Adventitial fibroblast
          "Dcn", "Pdgfra", "Pi16", "Fbln1",
          
          # Alveolar fibroblast
          "Npnt", "Fgfr4",
          
          # Inflammatory adventitial fibroblast
          "C3", "C1qa",
          
          # Pericyte
          "Myh11", "Pdgfrb", "Cspg4",
          
          # Smooth muscle
          "Acta2", "Hhip",
          
          # Myeloid
          "Ptprc",
          
          # Monocyte
          "Cd14", "Fcgr3", "Ccr1",
          
          # General macrophage
          "Adgre1", "Cd68",
          
          # M1 macrophage
          "Nos2", "Cd86", "CXCL9",
          
          # M2 macrophage
          "Pparg",  "Mrc1", "Fuca1",
          
          # Migratory macrophage
          "Itgal",  "Itga4" ,"Ncf2", "Rhoa",
          
          # Mesothelial
          "Lrrn4", "Upk3b",
          
          # Natural Killer
          #"Nkg7", 
          "Gzma", "Gzmb", "Klrd1", "Klrk1", 
          
          # General T
          "Cd3e", "Cd4", "Cd8a",
          
          # CD4+ effector
          "Nkg7", "Id2", "Cxcr6",
          
          # CD4+ effector intermediate
          
          # CD4+ memory or naive
          "Il7r", "Ccr7", "Sell",
          
          # Proliferating T
          "Mki67", "Top2a",
          
          # CD8+ Effector T
          
          # CD8+ memory or naive
          
          # Double negative
          "Zbtb16", "Rora", "Ikzf2",
          
          # T reg
          "Foxp3", "Ctla4",
          
          # Nuocyte
          "Icos", "Il17rb", "Gata3",
          
          # Inhibitory/exhausted T markers
          "Pdcd1", "Lag3", "Tigit", "Havcr2", "Tox"
          
          
        )),
        group.by = "cluster_Cell.Type", cluster.idents = F) + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))

ggsave(filename = paste0(output_full, "dotplot_cluster_celltype_markers.pdf"),
       width = 18, height = 8)


# Cell type-level
unique(so$Cell.Type)
so$Cell.Type <- factor(so$Cell.Type,
                       levels = c("AT1",
                                  "B",
                                  "Capillary endothelial",
                                  "Car4+ endothelial",
                                  "Lymphatic endothelial",
                                  "Venous endothelial",
                                  "Adventitial fibroblast",
                                  "Alveolar fibroblast",
                                  "Inflammatory adventitial fibroblast",
                                  "Pericyte",
                                  "Smooth muscle",
                                  "Monocyte",
                                  "M1 macrophage",
                                  "M2 macrophage",
                                  "Migratory macrophage",
                                  "Mesothelial",
                                  "Natural Killer",
                                  "CD4+ effector T",
                                  "CD4+ memory T",
                                  "Proliferating T",
                                  "CD8+ effector T",
                                  "CD8+ memory T",
                                  "Double negative T",
                                  "Regulatory T",
                                  "Nuocyte"
                       ))
table(so$Cell.Type)
DotPlot(so,
        features = str_to_title(c(
          
          #AT1
          "Ager","Gprc5a",
          
          #AT2
          "Sftpc","Sftpa1",
          
          # B (6,7,9)
          "Cd79a", "Ighd",
          
          # General endothelial
          "Cdh5","Cldn5",
          
          # Capillary endothelial
          "Gpihbp1", "Ptprb", "Kit",
          
          # Car4+ endothelial
          "Car4", "Cd34",
          
          # Lymphatic endothelial
          "Mmrn1", "Ccl21a", "Flt4",
          
          # Venous endothelial
          "Vwf", "Ackr1",
          
          # General fibroblast (specifics are a combo of the following markers)
          "Col1a1", 
          
          # Adventitial fibroblast
          "Dcn", "Pdgfra", "Pi16", "Fbln1",
          
          # Alveolar fibroblast
          "Npnt", "Fgfr4",
          
          # Inflammatory adventitial fibroblast
          "C3", "C1qa",
          
          # Pericyte
          "Myh11", "Pdgfrb", "Cspg4",
          
          # Smooth muscle
          "Acta2", "Hhip",
          
          # Myeloid
          "Ptprc",
          
          # Monocyte
          "Cd14", "Fcgr3", "Ccr1",
          
          # General macrophage
          "Adgre1", "Cd68",
          
          # M1 macrophage
          "Nos2", "Cd86", "CXCL9",
          
          # M2 macrophage
          "Pparg",  "Mrc1", "Fuca1",
          
          # Migratory macrophage
          "Itgal",  "Itga4" ,"Ncf2", "Rhoa",
          
          # Mesothelial
          "Lrrn4", "Upk3b",
          
          # Natural Killer
          #"Nkg7", 
          "Gzma", "Gzmb", "Klrd1", "Klrk1", 
          
          # General T
          "Cd3e", "Cd4", "Cd8a",
          
          # CD4+ effector
          "Nkg7", "Id2", "Cxcr6",
          
          # CD4+ effector intermediate
          
          # CD4+ memory or naive
          "Il7r", "Ccr7", "Sell",
          
          # Proliferating T
          "Mki67", "Top2a",
          
          # CD8+ Effector T
          
          # CD8+ memory or naive
          
          # Double negative
          "Zbtb16", "Rora", "Ikzf2",
          
          # T reg
          "Foxp3", "Ctla4",
          
          # Nuocyte
          "Icos", "Il17rb", "Gata3",
          
          # Inhibitory/exhausted T markers
          "Pdcd1", "Lag3", "Tigit", "Havcr2", "Tox"
          
          
        )),
        group.by = "Cell.Type", cluster.idents = F) + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))

ggsave(filename = paste0(output_full, "dotplot_celltype_markers.pdf"),
       width = 18, height = 8)

# Same plot as above but by sample
output_by_sample_dotplot <- paste0(output_full, "/dotplot_by_sample/")
dir.create(output_by_sample_dotplot)
for(sample in unique(so$Sample)){
  so_sample <- subset(so, Sample == sample)
  ggd <- DotPlot(so_sample,
                 features = str_to_title(c(
                   
                   #AT1
                   "Ager","Gprc5a",
                   
                   #AT2
                   "Sftpc","Sftpa1",
                   
                   # B (6,7,9)
                   "Cd79a", "Ighd",
                   
                   # General endothelial
                   "Cdh5","Cldn5",
                   
                   # Capillary endothelial
                   "Gpihbp1", "Ptprb", "Kit",
                   
                   # Car4+ endothelial
                   "Car4", "Cd34",
                   
                   # Lymphatic endothelial
                   "Mmrn1", "Ccl21a", "Flt4",
                   
                   # Venous endothelial
                   "Vwf", "Ackr1",
                   
                   # General fibroblast (specifics are a combo of the following markers)
                   "Col1a1", 
                   
                   # Adventitial fibroblast
                   "Dcn", "Pdgfra", "Pi16", "Fbln1",
                   
                   # Alveolar fibroblast
                   "Npnt", "Fgfr4",
                   
                   # Inflammatory adventitial fibroblast
                   "C3", "C1qa",
                   
                   # Pericyte
                   "Myh11", "Pdgfrb", "Cspg4",
                   
                   # Smooth muscle
                   "Acta2", "Hhip",
                   
                   # Myeloid
                   "Ptprc",
                   
                   # Monocyte
                   "Cd14", "Fcgr3", "Ccr1",
                   
                   # General macrophage
                   "Adgre1", "Cd68",
                   
                   # M1 macrophage
                   "Nos2", "Cd86", "CXCL9",
                   
                   # M2 macrophage
                   "Pparg",  "Mrc1", "Fuca1",
                   
                   # Migratory macrophage
                   "Itgal",  "Itga4" ,"Ncf2", "Rhoa",
                   
                   # Mesothelial
                   "Lrrn4", "Upk3b",
                   
                   # Natural Killer
                   #"Nkg7", 
                   "Gzma", "Gzmb", "Klrd1", "Klrk1", 
                   
                   # General T
                   "Cd3e", "Cd4", "Cd8a",
                   
                   # CD4+ effector
                   "Nkg7", "Id2", "Cxcr6",
                   
                   # CD4+ effector intermediate
                   
                   # CD4+ memory or naive
                   "Il7r", "Ccr7", "Sell",
                   
                   # Proliferating T
                   "Mki67", "Top2a",
                   
                   # CD8+ Effector T
                   
                   # CD8+ memory or naive
                   
                   # Double negative
                   "Zbtb16", "Rora", "Ikzf2",
                   
                   # T reg
                   "Foxp3", "Ctla4",
                   
                   # Nuocyte
                   "Icos", "Il17rb", "Gata3",
                   
                   # Inhibitory/exhausted T markers
                   "Pdcd1", "Lag3", "Tigit", "Havcr2", "Tox"
                   
                   
                 )),
                 group.by = "Cell.Type", cluster.idents = F) + 
    theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
  
  ggsave(filename = paste0(output_by_sample_dotplot, sample, "_dotplot_cluster_celltype_markers.pdf"),
         width = 22, height = 10,
         plot = ggd)
}

################################################################################
####### -------             Save labeled RDS object               ------- #######

#saveRDS(so,
#         "Data/Seurat_objects/AggrNoNorm_filtered_cluster_clean_labeled.RDS")
# so <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_cluster_clean_labeled.RDS")



################################################################################
####### -------            Post filter tables                    ------- #######

# Recalculating these tables because cells were removed

output_postfilter <- "Analysis/postFilter_QC"

# Dataset level
table_filter <- cbind.data.frame(Metric = c("Number of cells", 
                                            "Number of genes",
                                            "Median UMI per cell",
                                            "Median genes per cell"),
                                 Value = c(ncol(so),
                                           nrow(so),
                                           median(so$nCount_RNA),
                                           median(so$nFeature_RNA)))
write.table(table_filter,
            file = paste0(output_postfilter, "/dataset-level_metrics.txt"),
            sep = "\t", row.names = F, quote = F)

# Sample-level
table_sample <- cbind.data.frame(Sample = c("Number of cells", 
                                            "Number of genes",
                                            "Median UMI per cell",
                                            "Median genes per cell"))
for(samp in unique(so$Sample)){
  print(samp)
  
  so_samp <- subset(so,
                    Sample == samp)
  so_samp <- so_samp[rowSums(so_samp) > 0,]
  
  table_sample$samp <- c(ncol(so_samp),
                         nrow(so_samp),
                         median(so_samp$nCount_RNA),
                         median(so_samp$nFeature_RNA))
  colnames(table_sample)[ncol(table_sample)] <- samp
}

table_sample <- data.frame(t(table_sample))
write.table(table_sample,
            file = paste0(output_postfilter, "/sample-level_metrics.txt"),
            sep = "\t", row.names = T, col.names = F, quote = F)


################################################################################
####### -------     Cell type proportions by sample              ------- #######

# Compare distribution of new clusters and samples
prop_sample_celltype <- data.frame(prop.table(table(so@meta.data[,"Cell.Type"], 
                                                    so@meta.data[,"Sample"]), margin = 2))
gg_stackbar_sample <- ggplot(prop_sample_celltype, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +ylab("Cell type proportion") + 
  xlab("") + labs(fill="Cell type") +
  theme_classic() + theme(axis.text.x=element_text(angle=90,
                                                   hjust=1, vjust = 0.5, size = 12))
gg_stackbar_sample
ggsave(paste0(output_full, "/", "Celltype_stacked_barplot.pdf"),
       gg_stackbar_sample,
       width = 7, height = 6)

# Same but styled by proportion heatmap
dir.create(paste0(output_full, "proportions_celltype_by_sample/"))
prop_hm(so = so,
        group1 = "Sample",
        group2 = "Cell.Type",
        order_by = T,
        outdir = paste0(output_full, "proportions_celltype_by_sample/"),
        breaks = seq(0, 1, by = 0.01), # Sets max color to 1
        height = 8,
        width = 10,
        decimals = 4
)

dir.create(paste0(output_full, "proportions_broadcelltype_by_sample/"))
prop_hm(so = so,
        group1 = "Sample",
        group2 = "Broad.Cell.Type",
        order_by = T,
        outdir = paste0(output_full, "proportions_broadcelltype_by_sample/"),
        breaks = seq(0, 1, by = 0.01), # Sets max color to 1
        height = 8,
        width = 8,
        decimals = 4
)


#### Long version of cell type propotions 
prop_df <- read.table(file = paste0(output_full, "proportions_celltype_by_sample/proportion_table.txt"),
                      sep = "\t", header = T, stringsAsFactors = F,
                      check.names = F)
prop_df <- melt(prop_df)
colnames(prop_df) <- c("Cell.Type", "Sample", "Proportion")
prop_df$Group <- gsub("_rep1|_rep2", "", prop_df$Sample)


################################################################################
####### -------  Myeloid cell type proportions by select samples  ------ #######

dir.create(paste0(output_full, "proportions_myeloid_celltypes_by_select_samples/"))
so_mye <- subset(so,
                 Broad.Cell.Type %in% c("Macrophage", "Monocyte"))
unique(so_mye$Group)
so_mye <- subset(so_mye,
                 Group %in% c("Her2+IAV+anti-CD4_dpi15", "Her2+IAV+anti-CD4_dpi9",
                              "Her2+IAV+ant-CD4_dpi9"),
                 invert=T)
unique(so_mye$Group)
unique(so_mye$Sample)
so_mye$Cell.Type <- as.character(so_mye$Cell.Type)
table(so_mye$Cell.Type)
so_mye$Sample <- factor(so_mye$Sample,
                        levels = c(
                          "Her2+PBS_rep1",
                          "Her2+PBS_rep2",
                          "WT+IAV_dpi9_rep1",
                          "Her2+IAV_dpi9_rep1",
                          "Her2+IAV_dpi9_rep2",
                          "WT+IAV_dpi15_rep1",
                          "Her2+IAV_dpi15_rep1",
                          "Her2+IAV_dpi15_rep2"
                        ))
prop_hm(so = so_mye,
        group1 = "Sample",
        group2 = "Cell.Type",
        order_by = T,
        outdir = paste0(output_full, "proportions_myeloid_celltypes_by_select_samples/"),
        breaks = seq(0, 1, by = 0.01), # Sets max color to 1
        height = 4,
        width = 7,
        decimals = 3
)


################################################################################
####### -------     T cell type proportions by select samples    ------- #######

dir.create(paste0(output_full, "proportions_T_celltypes_by_select_samples/"))
so_t <- subset(so,
               Broad.Cell.Type %in% c("T"))
unique(so_t$Group)
so_t <- subset(so_t,
               Group %in% c("WT+IAV_dpi9", "WT+IAV_dpi15",
                            "Her2+IAV_dpi9", "Her2+IAV_dpi15")
)
unique(so_t$Group)
unique(so_t$Sample)
so_t$Cell.Type <- as.character(so_t$Cell.Type)
table(so_t$Cell.Type)
so_t$Sample <- factor(so_t$Sample,
                      levels = c(
                        "WT+IAV_dpi9_rep1",
                        "Her2+IAV_dpi9_rep1",
                        "Her2+IAV_dpi9_rep2",
                        "WT+IAV_dpi15_rep1",
                        "Her2+IAV_dpi15_rep1",
                        "Her2+IAV_dpi15_rep2"
                      ))
prop_hm(so = so_t,
        group1 = "Sample",
        group2 = "Cell.Type",
        order_by = T,
        outdir = paste0(output_full, "proportions_T_celltypes_by_select_samples/"),
        breaks = seq(0, 1, by = 0.01), # Sets max color to 1
        height = 4,
        width = 7,
        decimals = 3
)


################################################################################
####### -------      Select comparisons in select cell types     ------- #######

output_comparisons <- paste0(output_full, "select_comparisons/")
dir.create(output_comparisons)

# Define groups to compare and in what cell types
# Here, I will be subsetting to 2 groups and 1 set of BROAD cell types
# comparisons will be performed in more specific cell types within that broad subset
unique(so$Group)
table(so$Broad.Cell.Type)
group1s <- c("Her2+IAV_dpi9", "Her2+IAV_dpi15", 
             "Her2+IAV_dpi9", "Her2+IAV_dpi15",
             "Her2+IAV+anti-CD4_dpi9", "Her2+IAV+anti-CD4_dpi15")
group2s <- c("WT+IAV_dpi9", "WT+IAV_dpi15", 
             "WT+IAV_dpi9", "WT+IAV_dpi15",
             "Her2+IAV_dpi9", "Her2+IAV_dpi15")
cell_subsets <- c("myeloid", "myeloid", 
                  "t", "t", 
                  "t", "t") # These need to be defined in the if statements below

# Define gene sets for GSEA
gseapathways1 <- c("HALLMARK_IL2_STAT5_SIGNALING",
                   "HALLMARK_MTORC1_SIGNALING",
                   "HALLMARK_MYC_TARGETS_V1",
                   "KEGG_OXIDATIVE_PHOSPHORYLATION")
gseapathways2 <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE",
                   "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                   "KEGG_OXIDATIVE_PHOSPHORYLATION")
gene_sets <- all[all$gs_cat %in% "H" | all$gs_subcat %in% c("CP:KEGG"),] %>% 
  dplyr::select(gs_name, gene_symbol) %>% as.data.frame()


# Loop through each comparison
for(i in 1:length(group1s)){
  
  print(i)
  group1 <- group1s[i]
  group2 <- group2s[i]
  cell_subset <- cell_subsets[i]
  
  comparison_name <- paste0(group1, "_vs_", group2)
  
  ####### -------    
  # Subset object to comparison and select cell types
  ####### -------    
  so_comp <- subset(so,
                    Group %in% c(group1, group2))
  if(cell_subset == "myeloid"){
    so_comp <- subset(so_comp,
                      Broad.Cell.Type %in% c("Monocyte", "Macrophage"))
  }
  if(cell_subset == "t"){
    so_comp <- subset(so_comp,
                      Broad.Cell.Type %in% "T")
  }
  
  ####### -------    
  # Looping through each Cell.Type within the Broad.Cell.Type
  ####### -------    
  so_comp$Cell.Type <- as.character(so_comp$Cell.Type)
  for(celltype in unique(so_comp$Cell.Type)){
    
    
    ####### -------   Comparison output
    print(celltype)
    output_comparisons_specific <- paste0(output_comparisons,
                                          comparison_name, "_in_",
                                          gsub(" ", "_", celltype), "/")
    dir.create(output_comparisons_specific)
    
    so_comp_celltype <- subset(so_comp,
                               Cell.Type == celltype)
    
    
    ####### -------    DEG 
    comp_celltype_distribution <- data.frame(unclass(table(so_comp_celltype$Group)))
    comp_celltype_distribution <- cbind.data.frame(Group = rownames(comp_celltype_distribution),
                                                   Cells = comp_celltype_distribution[,1])
    dge <- FindMarkers(so_comp_celltype, 
                       ident.1 = group1, 
                       ident.2 = group2, 
                       only.pos = F, 
                       min.pct = 0.1, 
                       logfc.threshold = 0.1, # Lowered 
                       group.by = "Group"
    )
    dge <- cbind.data.frame(Gene = rownames(dge),
                            dge)
    dge_sig <- dge[dge$p_val_adj < 0.05,]
    write.table(dge,
                file = paste0(output_comparisons_specific, "dge_full_", 
                              comparison_name, "_in_",
                              gsub(" ", "_", celltype), ".txt"),
                sep = "\t", row.names = F, quote = F)
    write.table(dge_sig,
                file = paste0(output_comparisons_specific, "dge_sig_", 
                              comparison_name, "_in_",
                              gsub(" ", "_", celltype), ".txt"),
                sep = "\t", row.names = F, quote = F)
    write.table(comp_celltype_distribution,
                file = paste0(output_comparisons_specific, "cell_numbers__", 
                              comparison_name, "_in_",
                              gsub(" ", "_", celltype), ".txt"),
                sep = "\t", row.names = F, quote = F)
    
    
    ####### -------    Also calculate and save logFC after removing FC and pct filters
    fc_df <- FoldChange(so_comp_celltype, 
                        ident.1 = group1, 
                        ident.2 = group2, 
                        group.by = "Group",
                        mean.fxn = function(x) {
                          return(log(x = rowMeans(x = expm1(x = x)) + 1, base =2))}
                        # See https://github.com/satijalab/seurat/issues/6701
    )
    fc_df <- cbind.data.frame(Gene = rownames(fc_df),
                              fc_df)
    # Removing gene sets where there were 0 counts (would be ties in the middle)
    fc_df <- fc_df[!(fc_df$avg_log2FC == 0 & fc_df$pct.1 == 0 & fc_df$pct.2 == 0),]
    
    
    ####### -------    Run GSEA and plot two sets of gene sets
    ## Ranked vector
    ranked <- as.numeric(fc_df$avg_log2FC)
    names(ranked) <- fc_df$Gene
    ranked <- ranked[order(ranked,
                           decreasing = T)]
    ## Run GSEA
    edo2 <- GSEA(geneList = ranked,
                 TERM2GENE = gene_sets,
                 pvalueCutoff = 1,eps = 0)
    edo2_df <- data.frame(edo2)
    edo2_df <- edo2_df[,-2]
    write.table(edo2_df,
                paste0(output_comparisons_specific, "GSEA_results_",
                       comparison_name, "_in_",
                       gsub(" ", "_", celltype), ".txt"),
                sep = "\t", row.names = F, quote = F)
    ## GSEA multi plot 1
    gg_multi1 <- gseaplot2_ag(edo2, 
                              geneSetID = which(edo2@result$ID %in% gseapathways1), 
                              base_size = 20, subplots = c(1:3),
                              pvalue_table = TRUE, 
                              title = paste0(comparison_name, " in ", celltype),
                              fc_high = group1, fc_low = group2,
                              fc_label_size = 6,
                              title_size = 14
    )
    ggsave(file = paste0(output_comparisons_specific, "GSEA_plot1_",
                         comparison_name, "_in_",
                         gsub(" ", "_", celltype), ".pdf"),
           print(gg_multi1), width=12, height=6)
    ## GSEA multi plot 2
    gg_multi2 <- gseaplot2_ag(edo2, 
                              geneSetID = which(edo2@result$ID %in% gseapathways2), 
                              base_size = 20, subplots = c(1:3),
                              pvalue_table = TRUE, 
                              title = paste0(comparison_name, " in ", celltype),
                              fc_high = group1, fc_low = group2,
                              fc_label_size = 6,
                              title_size = 14
    )
    ggsave(file = paste0(output_comparisons_specific, "GSEA_plot2_",
                         comparison_name, "_in_",
                         gsub(" ", "_", celltype), ".pdf"),
           print(gg_multi2), width=12, height=6)
    
    
    ####### -------    Heatmap of top 10 genes in each direction
    deg_very_sig <- dge_sig[which(dge_sig$p_val_adj < 0.05),] # Meher had 10^-20
    if(nrow(deg_very_sig) > 0){
      deg_very_sig$cluster <- group1
      deg_very_sig$cluster[deg_very_sig$avg_log2FC < 0] <- group2
      dge_sig_top <- deg_very_sig %>%
        group_by(cluster) %>%
        slice_max(n = 10, order_by = abs(avg_log2FC) )
      so_comp_celltype_scale <- ScaleData(object = so_comp_celltype, 
                                          features = unique(dge_sig_top$Gene))
      Idents(so_comp_celltype_scale) <- so_comp_celltype_scale$Group
      p2 <- DoHeatmap(so_comp_celltype_scale, 
                      features = unique(dge_sig_top$Gene), group.by = "Group", slot = "scale.data",
                      angle = 0,hjust = 0.5,size = 4, raster = F) +
        theme(plot.margin = margin(0.8,0.8,0.8,0.8, "cm")) +
        theme(title=element_text(size=10,hjust = 0.5),legend.title = element_text(size=8),legend.text = element_text(size=7), legend.key.size = unit(0.5, 'cm')) + guides(color="none") +
        ggtitle(paste0(comparison_name, " in ", celltype)) +
        scale_fill_gradient2(low = rev(c('skyblue','skyblue3','skyblue4')), mid = "white", high = rev(c('sienna3','sienna2','sienna1')), midpoint = 0, guide = "colourbar", aesthetics = "fill")
      p2
      ggsave(paste0(output_comparisons_specific, "Heatmap_top10_",
                    comparison_name, "_in_",
                    gsub(" ", "_", celltype), ".pdf"),
             plot = print(p2))
    }
    
    
  }
}


################################################################################
####### -------              DEG tables for Afshin               ------- #######

output_comparisons <- paste0(output_full, "select_comparisons_Afshin/")
dir.create(output_comparisons)

# DEGs for the following comparisons:
# Her2 flu vs Her2 PBS (9 and 15 dpi)
# Her2 flu vs Her2 WT (9 and 15 dpi)
# Her2 flu anti-CD4 vs Her2 flu (9 and 15 dpi)

unique(so$Sample)
so$Sample_drop_rep <- so$Sample
so$Sample_drop_rep <- gsub("_rep1|_rep2", "", so$Sample_drop_rep)
unique(so$Sample_drop_rep)
group1s <- c("Her2+IAV_dpi9", "Her2+IAV_dpi15", 
             "Her2+IAV_dpi9", "Her2+IAV_dpi15",
             "Her2+IAV+anti-CD4_dpi9", "Her2+IAV+anti-CD4_dpi15")
group2s <- c("Her2+PBS", "Her2+PBS", 
             "WT+IAV_dpi9", "WT+IAV_dpi15",
             "Her2+IAV_dpi9", "Her2+IAV_dpi15")

####### -------
# Loop through each comparison
####### -------
for(i in 1:length(group1s)){
  
  print(i)
  group1 <- group1s[i]
  group2 <- group2s[i]
  
  comparison_name <- paste0(group1, "_vs_", group2)
  
  comparison_output <- paste0(output_comparisons,
                              comparison_name, "/")
  dir.create(comparison_output)
  
  ####### -------    
  # Subset object to comparison and select cell types
  ####### -------    
  so_comp <- subset(so,
                    Sample_drop_rep %in% c(group1, group2))
  
  ####### -------    
  # Looping through each Cell.Type
  ####### -------    
  so_comp$Cell.Type <- as.character(so_comp$Cell.Type)
  for(celltype in unique(so_comp$Cell.Type)){
    
    so_comp_celltype <- subset(so_comp,
                               Cell.Type == celltype)
    
    # Skip comparison if there isn't at least 1 cell in both groups
    cell_numbers <- data.frame(unclass(table(so_comp_celltype$Sample_drop_rep)))
    if(sum(cell_numbers[,1] >= 3) == 2){
      dge <- FindMarkers(so_comp_celltype, 
                         ident.1 = group1, 
                         ident.2 = group2, 
                         only.pos = F, 
                         min.pct = 0, # Lowered 
                         logfc.threshold = 0, # Lowered 
                         min.cells.feature = 0,
                         min.cells.group = 3,
                         group.by = "Sample_drop_rep"
      )
      dge <- cbind.data.frame(Gene = rownames(dge),
                              dge)
      write.table(dge,
                  file = paste0(comparison_output, gsub(" ", "_", celltype), ".txt"),
                  sep = "\t", row.names = F, quote = F)
    }
    
    
  }
}


####### -------
# Reformat the above files for a collaborator
####### -------
merged_table <- data.frame()
for(i in 1:length(group1s)){
  
  print(i)
  group1 <- group1s[i]
  group2 <- group2s[i]
  
  comparison_name <- paste0(group1, "_vs_", group2)
  
  comparison_output <- paste0(output_comparisons,
                              comparison_name, "/")
  
  ####### -------    
  # Looping through each Cell.Type (file)
  ####### -------    
  celltypes <- list.files(comparison_output)
  celltypes <- gsub(".txt", "", celltypes)
  for(celltype in celltypes){
    
    dge <- read.table(paste0(comparison_output, gsub(" ", "_", celltype), ".txt"),
                      sep = "\t", header = T, stringsAsFactors = F)
    dge <- dge[,c("Gene", "p_val", "avg_log2FC", "p_val_adj")]
    colnames(dge) <- c("Mouse.Gene", "pval", "log2FC", "adj.p")
    colnames(dge)[2:4] <- paste0(celltype, ":", comparison_name, 
                                 ":", colnames(dge)[2:4])
    
    if(nrow(merged_table) == 0){
      merged_table <- dge
    }else{
      merged_table <- join(merged_table, dge)
    }
    
  }
}

# Create metadata table
merged_table_metadata <- cbind.data.frame(column_name = 
                                            colnames(merged_table)[2:ncol(merged_table)])
merged_table_metadata$Cell.Type <- gsub("_", " ", sapply(strsplit(merged_table_metadata$column_name,
                                                                  split = ":"),
                                                         function(x){x[[1]]}))
merged_table_metadata$Comparison <- gsub("_vs_", " vs ", sapply(strsplit(merged_table_metadata$column_name,
                                                                         split = ":"),
                                                                function(x){x[[2]]}))
merged_table_metadata$Metric <- sapply(strsplit(merged_table_metadata$column_name,
                                                split = ":"),
                                       function(x){x[[3]]})
merged_table_metadata$Day <- "dpiX"
merged_table_metadata$Day[grep("dpi9", merged_table_metadata$column_name)] <- "dpi9"
merged_table_metadata$Day[grep("dpi15", merged_table_metadata$column_name)] <- "dpi15"
merged_table_metadata$Day2 <- merged_table_metadata$Day
merged_table_metadata$Day2 <- gsub("dpi9", "9dpi", merged_table_metadata$Day2)
merged_table_metadata$Day2 <- gsub("dpi15", "15dpi", merged_table_metadata$Day2)

# Order by metric 
merged_table_metadata$Metric <- factor(merged_table_metadata$Metric,
                                       levels = c("log2FC", "pval", "adj.p"))
merged_table_metadata <- merged_table_metadata[order(merged_table_metadata$Metric),]
merged_table <- merged_table[,match(c("Mouse.Gene", merged_table_metadata$column_name),
                                    colnames(merged_table))]

write.xlsx(merged_table,
           file = paste0(output_comparisons, "/merged_comparisons.xlsx"),
           overwrite = T, row.names =F, col.names=T)
write.xlsx(merged_table_metadata,
           file = paste0(output_comparisons, "/merged_comparisons_column_metadata.xlsx"),
           overwrite = T, row.names =F, col.names=T)

# Another version where I try to match the format as closely as possible
merged_table_match_format <- rbind(c("", as.character(merged_table_metadata$Metric)),
                                   c("", paste0(merged_table_metadata$Cell.Type, ":",
                                                merged_table_metadata$Comparison, ":",
                                                merged_table_metadata$Day2)),
                                   c("", paste0(as.character(merged_table_metadata$Metric), ":",
                                                merged_table_metadata$Cell.Type, ":",
                                                merged_table_metadata$Comparison, ":",
                                                merged_table_metadata$Day2)),
                                   merged_table)
write.xlsx(merged_table_match_format,
           file = paste0(output_comparisons, "/merged_comparisons_formatted.xlsx"),
           overwrite = T, row.names =F, col.names=T)


################################################################################
####### -------    IFNa and IFNy ORA dpi15, many cell types      ------- #######

output_ifn <- paste0(output_full, "ORA_IFNa_IFNy_dpi15/")
dir.create(output_ifn)

# Subset to 3 groups
unique(so$Broad.Group)
groups_of_interest <- c("Her2+PBS", "WT+IAV_dpi15", "Her2+IAV_dpi15")
so_sub <- subset(so,
                 Broad.Group %in% groups_of_interest)
unique(so_sub$Broad.Group)
so_sub$Broad.Group <- factor(so_sub$Broad.Group,
                             levels = groups_of_interest)

# Prep gene sets
genesets_of_interest <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                          "HALLMARK_INTERFERON_GAMMA_RESPONSE")
# Could be interesting: HALLMARK_TNFA_SIGNALING_VIA_NFKB, GOBP_CELL_ACTIVATION
# GOBP_IMMUNE_EFFECTOR_PROCESS, GOBP_RESPONSE_TO_CYTOKINE
# gene_sets <- all[all$gs_name %in% genesets_of_interest,] 
gene_sets <- all[all$gs_cat %in% "H" | all$gs_subcat %in% c("CP:KEGG", "GO:BP"),]
length(unique(gene_sets$gs_name))
gene_sets <- gene_sets %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

# ORA for 2 gene sets in each cell type
ora_merged <- data.frame() # Collecting all ORA results for merged heatmaps
for(celltype in unique(so_sub$Cell.Type)){
#for(celltype in unique(so_sub$Cell.Type)[14:24]){
  
  # Subset by celltype
  print(celltype)
  so_sub_celltype <- subset(so_sub,
                            Cell.Type == celltype)
  output_ifn_celltype <- paste0(output_ifn, celltype, "/")
  dir.create(output_ifn_celltype)
  
  # Find all markers across Broad.Group
  fam <- enriched_genes_by_group(so = so_sub_celltype,
                                 group = "Broad.Group",
                                 outdir = output_ifn_celltype,
                                 n_genes = 200)
  
  # ORA
  ora_table <- ora( gene_sets = gene_sets,
                    gene_list = fam,
                    background = rownames(so),
                    outdir = output_ifn_celltype,
                    targeted_gene_sets_heatmap = genesets_of_interest,
                    cluster_cols = F,
                    top_n_gene_sets_heatmap = 5,
                    cluster_rows = F,
                    fontsize = 4
  )
  
  # Subset table to genesets_of_interest
  ora_table_sub <- ora_table[ora_table$ID %in% genesets_of_interest,]
  if(nrow(ora_table_sub) == 0){
    ora_table_sub <- cbind.data.frame(celltype = celltype,
                                      group = unique(so_sub_celltype$Broad.Group),
                                      geneset = genesets_of_interest[1],
                                      neg_log10_fdr = 0)
  }else{
    ora_table_sub <- cbind.data.frame(celltype = celltype,
                                    group = ora_table_sub$Group,
                                    geneset = ora_table_sub$ID,
                                    neg_log10_fdr = ora_table_sub$neg_log10_fdr)
  }
  
  # Fill in no gene overlap with 0
  missing_combinations <- apply(expand.grid(groups_of_interest, genesets_of_interest), 1, paste, collapse=".")
  missing_combinations <- missing_combinations[!(missing_combinations %in% 
                                                   paste0(ora_table_sub$group, ".",
                                                          ora_table_sub$geneset))]
  if(length(missing_combinations) > 0){
    ora_table_sub <- rbind(ora_table_sub,
                           cbind.data.frame(celltype = celltype,
                                            group = sapply(strsplit(missing_combinations, split = "\\."),
                                                           function(x){x[[1]]}),
                                            geneset = sapply(strsplit(missing_combinations, split = "\\."),
                                                             function(x){x[[2]]}),
                                            neg_log10_fdr = 0))
  }
  
  # Merge ORA
  if(nrow(ora_merged) == 0){
    ora_merged <- ora_table_sub
  }else{
    ora_merged <- rbind(ora_merged, ora_table_sub)
  }
  
}

# Function to save heatmap
save_pheatmap_pdf <- function(x, filename, width=width, height=height) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  grid.text("-log10(FDR)" , x= 0.93, y = 1.1,
            gp=gpar(fontsize=6, col="black"),
            just = "left") # AG added 4/10/24. Not sure how generalizeable it is. 
  dev.off()
}

# Heatmap across cell type for each gene set
for(gs in genesets_of_interest){
  print(gs)
  
  ora_merged_gs <- ora_merged[ora_merged$geneset == gs,]
  
  ora_merged_gs_wide <- dcast(ora_merged_gs, celltype ~ group,
                              value.var = "neg_log10_fdr")
  rownames(ora_merged_gs_wide) <- ora_merged_gs_wide$celltype
  ora_merged_gs_wide <- ora_merged_gs_wide[,groups_of_interest]
  
  targeted_heatmap_table_sig <- ora_merged_gs_wide[] > 1.30103 # Equiv to FDR 0.05
  targeted_heatmap_table_sig <- apply(targeted_heatmap_table_sig, 2, 
                                      function(x){gsub("TRUE", "*", x)})
  targeted_heatmap_table_sig <- apply(targeted_heatmap_table_sig, 2, 
                                      function(x){gsub("FALSE", "", x)})
  
  # Plot as heatmap
  ph <- pheatmap(ora_merged_gs_wide,
                 fontsize = 6,
                 color = c("white", colorRampPalette(c("#F3F2FE", "#0E0AFA"))(100)),
                 display_numbers = targeted_heatmap_table_sig,
                 fontsize_number = 20,
                 number_color = "black",
                 cluster_cols = F,
                 cluster_rows = F,
                 main = gsub("_", " ", gs))
  save_pheatmap_pdf(ph,
                    paste0(output_ifn, "ORA_heatmap_", gs, "_across_cell_type.pdf"))
  write.table(cbind.data.frame(celltype = rownames(ora_merged_gs_wide),
                               ora_merged_gs_wide), 
              file = paste0(output_ifn, "table_ORA_heatmap_", gs, "_across_cell_type.txt"),
              sep = "\t", row.names = F, quote = F)
}



################################################################################
####### -------    IFNa and IFNy ORA dpi9, many cell types      ------- #######

output_ifn <- paste0(output_full, "ORA_IFNa_IFNy_dpi9/")
dir.create(output_ifn)

# Subset to 3 groups
unique(so$Broad.Group)
groups_of_interest <- c("Her2+PBS", "WT+IAV_dpi9", "Her2+IAV_dpi9")
so_sub <- subset(so,
                 Broad.Group %in% groups_of_interest)
unique(so_sub$Broad.Group)
so_sub$Broad.Group <- factor(so_sub$Broad.Group,
                             levels = groups_of_interest)

# Prep gene sets
genesets_of_interest <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                          "HALLMARK_INTERFERON_GAMMA_RESPONSE")
# Could be interesting: HALLMARK_TNFA_SIGNALING_VIA_NFKB, GOBP_CELL_ACTIVATION
# GOBP_IMMUNE_EFFECTOR_PROCESS, GOBP_RESPONSE_TO_CYTOKINE
# gene_sets <- all[all$gs_name %in% genesets_of_interest,] 
gene_sets <- all[all$gs_cat %in% "H" | all$gs_subcat %in% c("CP:KEGG", "GO:BP"),]
length(unique(gene_sets$gs_name))
gene_sets <- gene_sets %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

# ORA for 2 gene sets in each cell type
ora_merged <- data.frame() # Collecting all ORA results for merged heatmaps
for(celltype in unique(so_sub$Cell.Type)){
  
  # Subset by celltype
  print(celltype)
  so_sub_celltype <- subset(so_sub,
                            Cell.Type == celltype)
  output_ifn_celltype <- paste0(output_ifn, celltype, "/")
  dir.create(output_ifn_celltype)
  
  # Find all markers across Broad.Group
  fam <- enriched_genes_by_group(so = so_sub_celltype,
                                 group = "Broad.Group",
                                 outdir = output_ifn_celltype,
                                 n_genes = 200)
  
  # ORA
  ora_table <- ora( gene_sets = gene_sets,
                    gene_list = fam,
                    background = rownames(so),
                    outdir = output_ifn_celltype,
                    targeted_gene_sets_heatmap = genesets_of_interest,
                    cluster_cols = F,
                    top_n_gene_sets_heatmap = 5,
                    cluster_rows = F,
                    fontsize = 4
  )
  
  # Subset table to genesets_of_interest
  ora_table_sub <- ora_table[ora_table$ID %in% genesets_of_interest,]
  if(nrow(ora_table_sub) == 0){
    ora_table_sub <- cbind.data.frame(celltype = celltype,
                                      group = unique(so_sub_celltype$Broad.Group),
                                      geneset = genesets_of_interest[1],
                                      neg_log10_fdr = 0)
  }else{
    ora_table_sub <- cbind.data.frame(celltype = celltype,
                                      group = ora_table_sub$Group,
                                      geneset = ora_table_sub$ID,
                                      neg_log10_fdr = ora_table_sub$neg_log10_fdr)
  }
  
  # Fill in no gene overlap with 0
  missing_combinations <- apply(expand.grid(groups_of_interest, genesets_of_interest), 1, paste, collapse=".")
  missing_combinations <- missing_combinations[!(missing_combinations %in% 
                                                   paste0(ora_table_sub$group, ".",
                                                          ora_table_sub$geneset))]
  if(length(missing_combinations) > 0){
    ora_table_sub <- rbind(ora_table_sub,
                           cbind.data.frame(celltype = celltype,
                                            group = sapply(strsplit(missing_combinations, split = "\\."),
                                                           function(x){x[[1]]}),
                                            geneset = sapply(strsplit(missing_combinations, split = "\\."),
                                                             function(x){x[[2]]}),
                                            neg_log10_fdr = 0))
  }
  
  # Merge ORA
  if(nrow(ora_merged) == 0){
    ora_merged <- ora_table_sub
  }else{
    ora_merged <- rbind(ora_merged, ora_table_sub)
  }
  
}

## d9 did not produce any significance for this cell type
#ora_merged <- rbind(ora_merged,
#                    cbind.data.frame(celltype = "Lymphatic endothelial",
#                                     group = c(groups_of_interest, groups_of_interest),
#                                     geneset = c(rep(genesets_of_interest[1], length(groups_of_interest)),
#                                                 rep(genesets_of_interest[2], length(groups_of_interest))),
#                                     neg_log10_fdr = 0))

# Function to save heatmap
save_pheatmap_pdf <- function(x, filename, width=width, height=height) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  grid.text("-log10(FDR)" , x= 0.93, y = 1.1,
            gp=gpar(fontsize=6, col="black"),
            just = "left") # AG added 4/10/24. Not sure how generalizeable it is. 
  dev.off()
}

# Heatmap across cell type for each gene set
for(gs in genesets_of_interest){
  print(gs)
  
  ora_merged_gs <- ora_merged[ora_merged$geneset == gs,]
  
  ora_merged_gs_wide <- dcast(ora_merged_gs, celltype ~ group,
                              value.var = "neg_log10_fdr")
  rownames(ora_merged_gs_wide) <- ora_merged_gs_wide$celltype
  ora_merged_gs_wide <- ora_merged_gs_wide[,groups_of_interest]
  
  targeted_heatmap_table_sig <- ora_merged_gs_wide[] > 1.30103 # Equiv to FDR 0.05
  targeted_heatmap_table_sig <- apply(targeted_heatmap_table_sig, 2, 
                                      function(x){gsub("TRUE", "*", x)})
  targeted_heatmap_table_sig <- apply(targeted_heatmap_table_sig, 2, 
                                      function(x){gsub("FALSE", "", x)})
  
  # Plot as heatmap
  ph <- pheatmap(ora_merged_gs_wide,
                 fontsize = 6,
                 color = c("white", colorRampPalette(c("#F3F2FE", "#0E0AFA"))(100)),
                 display_numbers = targeted_heatmap_table_sig,
                 fontsize_number = 20,
                 number_color = "black",
                 cluster_cols = F,
                 cluster_rows = F,
                 main = gsub("_", " ", gs))
  save_pheatmap_pdf(ph,
                    paste0(output_ifn, "ORA_heatmap_", gs, "_across_cell_type.pdf"))
  write.table(cbind.data.frame(celltype = rownames(ora_merged_gs_wide),
                               ora_merged_gs_wide), 
              file = paste0(output_ifn, "table_ORA_heatmap_", gs, "_across_cell_type.txt"),
              sep = "\t", row.names = F, quote = F)
}


################################################################################
####### ------  See which cells express higher levels of Il6/Il6ra ------ #######

output_il6 <- paste0(output_full, "/IL6/")
dir.create(output_il6)

# IL6-related genes
rownames(so)[grep("Il6", rownames(so))]
genes <- c("Il6ra", "Il6", "Il6st")

# FindAllMarkers across cell type in full dataset
Idents(so) <- so$Cell.Type
fam <- FindAllMarkers(so,
               features = genes,
               only.pos = T,
               return.thresh = 1)
fam <- fam[fam$p_val_adj < 0.05,] # Only 1 was nonsig
write.table(fam, file = paste0(output_il6, "Il6_findallmarkers.txt"),
            sep = "\t", row.names = F, quote = F)

# dotplots
DotPlot(so,
        features = genes,
        group.by = "Cell.Type") + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
ggsave(filename = paste0(output_il6, "Il6_dotplot_all_cell_types.pdf"))


################################################################################
####### -------       Distribution of cell types by sample       ------- #######

# Compare distribution of new clusters and samples
celltype_cells <- data.frame(t(unclass(table(so@meta.data[,"Sample"], 
                                              so@meta.data[,"Cell.Type"]))))
celltype_cells <- cbind.data.frame(Cell.Type = rownames(celltype_cells), celltype_cells)
write.table(celltype_cells,
            file = paste0(output_full, "/", "celltype_by_sample.txt"),
            sep = "\t", row.names = F, col.names = T, quote = F)


################################################################################
####### -------            Saving 'final' RDS object             ------- #######

saveRDS(so,
         "Data/Seurat_objects/AggrNoNorm_filtered_cluster_clean_labeled.RDS")
so <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_cluster_clean_labeled.RDS")









