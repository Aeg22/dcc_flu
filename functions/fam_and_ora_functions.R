### Two functions:

# 1. enriched_genes_by_group() generates a list of enriched genes per group
#    Running this function is optional if you already have a list of gene vectors

# 2. ora() runs over representation analysis on the above list

### Function 1. enriched_genes_by_group() generates a list of enriched genes per group

enriched_genes_by_group <- function(...,
                                    so = NULL,
                                    group = NULL,
                                    outdir = NULL,
                                    min.pct = 0.1,
                                    logfc.threshold = 0.1,
                                    return.thresh = 0.05,
                                    n_genes = 200
                                    ){
  
  # This function will use Seurat::FindAllMarkers() to determine enriched genes
  # per group. All positively enriched genes with an adjusted p-value < `return.thresh`
  # will be considered. If the number of significant genes within a group is > `n_genes`, 
  # the top `n_genes` are chosen by highest fold-change. 
  #
  # so = NULL              # Seurat object. Prefilter if necessary
  # group = NULL           # Group column in `so` to calculate enrichment on
  # outdir = NULL          # (optional) Directory to save files sig_table.txt, full_table.txt, and dge_top.rds
  # min.pct=0.1            # Param in FindAllMarkers(): `min.pct`
  # logfc.threshold=0.1    # Param in FindAllMarkers(): `logfc.threshold`
  # return.thresh=0.05     # FDR cutoff to determine significant genes
  # n_genes=200            # Maximum number of genes per group to include in list output
  # ...additional arguments for FindAllMarkers()
  #
  # Output
  # 
  # List of enriched genes per group
  # Optional files to `outdir`
  #     dge_sig_table.txt       # txt of significant genes per group
  #     dge_full_table.txt      # txt of all tested genes per group
  #     dge_top.rds             # rds of gene enrichment list (main function output)
  
  library(Seurat)
  library(dplyr)
  
  if (is.null(so) | is.null(group)){
    stop("Both so and group are required")
  }
  
  # DGE and adj pval filter
  Idents(so) <- so@meta.data[,group]
  full_table <- FindAllMarkers(so, 
                               only.pos = TRUE, 
                               min.pct = min.pct, 
                               logfc.threshold = logfc.threshold,
                               verbose = T,
                               return.thresh = return.thresh,
                               ...)
  sig_table <- full_table[full_table$p_val_adj < return.thresh,] # Only keep sig genes
  
  print("Significant genes per group: ")
  print(table(sig_table$cluster))
  
  # Subset to top genes
  sig_table_top <- sig_table %>%
    group_by(cluster) %>%
    slice_max(n = n_genes, order_by = avg_log2FC) 
  
  # List
  sig_top_list <- split(sig_table_top$gene, sig_table_top$cluster)
  
  # Write files
  if(!is.null(outdir)){
    write.table(sig_table, file = paste0(outdir, "/dge_sig_table.txt"), 
                sep = "\t", row.names = F, quote = F)
    write.table(full_table, file = paste0(outdir, "/dge_full_table.txt"), 
                sep = "\t", row.names = F, quote = F)
    saveRDS(sig_top_list, paste0(outdir, "/dge_top.rds"))
    }
  
  # Return list
  return(sig_top_list)
  
}


### Function 2. ora() runs over representation analysis on the above list

ora <- function(...,
                gene_list = NULL,
                gene_sets = NULL,
                background = NULL,
                outdir = NULL,
                targeted_gene_sets_heatmap = NULL,
                top_n_gene_sets_heatmap = 2, # 0 to skip
                include_ties_for_top_n_gene_sets = T,
                adj_pval_table = 1,
                adj_pval_heatmap = 0.05,
                # Changed several pheatmap() defaults
                fontsize = 6,
                color = c("white", colorRampPalette(c("#F3F2FE", "#0E0AFA"))(100)),
                # Size of heatmap
                width = 4,
                height = 2,
                file_prefix = ""
                ){
  
  # This function will take in a list of vectors corresponding to genes to perform
  # over representation analysis on. Optionally, heatmaps displaying the -log10(fdr)
  # of select or top gene sets can also be generated. 
  #
  # Parameters
  # gene_list = NULL                      # Named list of vectors of genes
  # gene_sets = NULL                      # Data.frame of gene sets (columns: gs_name, gene_symbol)
  # background = NULL                     # Background genes for ORA testing
  # outdir = NULL                         # Directory to output file(s)
  # targeted_gene_sets_heatmap = NULL     # (optional) vector of gene sets with gene_sets for targeted heatmap 
  # top_n_gene_sets_heatmap = 2           # (optional, use 0 to skip) Number of top gene sets to display in heatmap per group
  # include_ties_for_top_n_gene_sets = T  # Include gene sets with ties for p-value (TRUE) or randomly select from the top (FALSE)
  # adj_pval_table = 1                    # FDR cutoff of merged ORA table. Default is to return all tests. 
  # adj_pval_heatmap = 0.05               # FDR cutoff to label significant gene sets in the heatmap
  # fontsize = 6                          # Heatmap fontsize
  # color = c("white", colorRampPalette(c("#F3F2FE", "#0E0AFA"))(100))
  #                                       # Heatmap color scheme (two colors makes more sense than three to me)
  # width = 4                             # Heatmap pdf width
  # height = 2                            # Heatmap pdf height
  # file_prefix                           # Prefix to add to all file outputs
  # ...additional arguments for pheatmap()
  #
  # Output
  # 
  # ORA_table.txt                             # txt of significant gene sets per group
  # ORA_heatmap_targeted_gene_sets.pdf        # Optional, heatmap of targeted gene sets
  # ORA_heatmap_targeted_gene_sets_values.txt # Optional, -log10(fdr) values of targeted heatmap
  # ORA_heatmap_top_gene_sets.pdf             # Optional, heatmap of top gene sets
  # ORA_heatmap_top_gene_sets_values.txt      # Optional, -log10(fdr) values of top heatmap
  
  # Default coloring used to be colorRampPalette(c("white", "orange", "darkorange", "red", "firebrick3", "darkred"))(100)
  
  library(clusterProfiler)
  library(pheatmap)
  library(grid)
  
  if (is.null(gene_list) | is.null(gene_sets) | is.null(background)| is.null(outdir)){
    stop("The following parameters are required: gene_list, gene_sets, background, outdir")
  }
  
  # Function to save heatmap
  save_pheatmap_pdf <- function(x, filename, width=width, height=height) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    grid.text("-log10(FDR)" , x= 0.93, y = 0.9,
              gp=gpar(fontsize=fontsize, col="black"),
              just = "left") # AG added 4/10/24. Not sure how generalizeable it is. 
    dev.off()
  }
  
  adj_pval_log10 <- -log10(adj_pval_heatmap)
  
  merged_ora_table <- data.frame()
  for(i in 1:length(gene_list)){
    
    genes <- gene_list[[i]]
    
    print(paste0("Calculating ORA for ", names(gene_list)[i]))
    
    ora_enricher <- enricher(genes, TERM2GENE=gene_sets, pAdjustMethod = "fdr", 
                             minGSSize = 1, maxGSSize = 1000000, 
                             qvalueCutoff = 1, pvalueCutoff = 1,
                             universe = background) 
    if(is.null(ora_enricher)){ # Skip to next i if there are no enricher results
      next
    }
    ora_result <- data.frame(ora_enricher)
    ora_result <- ora_result[,-2]
    ora_result <- cbind.data.frame(Group = names(gene_list)[i], ora_result)
    if(nrow(merged_ora_table) == 0){
      merged_ora_table <- ora_result
    }else{
      merged_ora_table <- rbind(merged_ora_table, ora_result)
    }
  }
  
  # Calculate -log10(FDR)
  merged_ora_table$neg_log10_fdr <- -log10(merged_ora_table$p.adjust)
  
  # Subset to significance 
  merged_ora_table_sig <- merged_ora_table[merged_ora_table$p.adjust <= adj_pval_table,]
  write.table(merged_ora_table_sig, file = paste0(outdir, file_prefix, "ORA_table.txt"),
              sep = "\t", row.names = F, quote = F)

  #### If a targeted gene set heatmap is desired
  if(!is.null(targeted_gene_sets_heatmap)){
    
    targeted_gene_sets_heatmap_in_data <- targeted_gene_sets_heatmap[targeted_gene_sets_heatmap %in% gene_sets$gs_name]
    if(!(identical(targeted_gene_sets_heatmap_in_data, targeted_gene_sets_heatmap))){
      warning(paste0("The following gene sets are expected to be in gene_sets but are not: ",
                     paste0(targeted_gene_sets_heatmap[!(targeted_gene_sets_heatmap %in% gene_sets$gs_name)], collapse = " ")))
    }
    
    targeted_heatmap_table <- cbind.data.frame(GeneSet = targeted_gene_sets_heatmap_in_data,
                                               holder = NA)
    merged_ora_table_sub <- merged_ora_table[,c("Group", "neg_log10_fdr", "ID")]
    for(group in names(gene_list)){
      merged_ora_table_sub_group <- merged_ora_table_sub[merged_ora_table_sub$Group == group,]
      targeted_heatmap_table$group <- merged_ora_table_sub_group$neg_log10_fdr[match(targeted_heatmap_table$GeneSet,
                                                     merged_ora_table_sub_group$ID)]
      colnames(targeted_heatmap_table)[ncol(targeted_heatmap_table)] <- group
    }
    targeted_heatmap_table_save <- targeted_heatmap_table[,-2]
    rownames(targeted_heatmap_table) <- targeted_heatmap_table$GeneSet
    targeted_heatmap_table <- targeted_heatmap_table[,-c(1,2)]
    targeted_heatmap_table[is.na(targeted_heatmap_table)] <- 0 # If it wasn't tested, there were no overlapping genes
    
    # Label significance 
    targeted_heatmap_table_sig <- targeted_heatmap_table[] > adj_pval_log10
    targeted_heatmap_table_sig <- apply(targeted_heatmap_table_sig, 2, 
                                        function(x){gsub("TRUE", "*", x)})
    targeted_heatmap_table_sig <- apply(targeted_heatmap_table_sig, 2, 
                                        function(x){gsub("FALSE", "", x)})
    
    # Plot as heatmap
    # Skip if all values are the same (happens in cases where the gene sets have 
    #   no overlap with any of the groups, meaning they are not in the ORA table at all,
    #   but note that they are in the tested gene sets)
    if(length(unique(as.vector(unlist(targeted_heatmap_table)))) > 1){
      ph <- pheatmap(targeted_heatmap_table,
                     fontsize = fontsize,
                     color = color,
                     display_numbers = targeted_heatmap_table_sig,
                     fontsize_number = 20,
                     number_color = "black",
                     ...)
      save_pheatmap_pdf(ph,
                        paste0(outdir, file_prefix, "ORA_heatmap_targeted_gene_sets.pdf")) 
    }
    write.table(targeted_heatmap_table_save, file = paste0(outdir, file_prefix, "ORA_heatmap_targeted_gene_sets_values.txt"),
                sep = "\t", row.names = F, quote = F)
    
  }
  
  
  
  #### If heatmap based on top gene sets is desired
  if(top_n_gene_sets_heatmap > 0){
    
    merged_ora_table_top <- merged_ora_table %>%
      group_by(Group) %>%
      slice_max(n = top_n_gene_sets_heatmap, order_by = neg_log10_fdr,
                with_ties = include_ties_for_top_n_gene_sets) 
    
    targeted_gene_sets_heatmap_in_data <- unique(merged_ora_table_top$ID)
    
    targeted_heatmap_table <- cbind.data.frame(GeneSet = targeted_gene_sets_heatmap_in_data,
                                               holder = NA)
    merged_ora_table_sub <- merged_ora_table[,c("Group", "neg_log10_fdr", "ID")]
    for(group in names(gene_list)){
      merged_ora_table_sub_group <- merged_ora_table_sub[merged_ora_table_sub$Group == group,]
      targeted_heatmap_table$group <- merged_ora_table_sub_group$neg_log10_fdr[match(targeted_heatmap_table$GeneSet,
                                                                                     merged_ora_table_sub_group$ID)]
      colnames(targeted_heatmap_table)[ncol(targeted_heatmap_table)] <- group
    }
    targeted_heatmap_table_save <- targeted_heatmap_table[,-2]
    rownames(targeted_heatmap_table) <- targeted_heatmap_table$GeneSet
    targeted_heatmap_table <- targeted_heatmap_table[,-c(1,2)]
    targeted_heatmap_table[is.na(targeted_heatmap_table)] <- 0 # If it wasn't tested, there were no overlapping genes
    
    # Label significance 
    targeted_heatmap_table_sig <- targeted_heatmap_table[] > adj_pval_log10
    targeted_heatmap_table_sig <- apply(targeted_heatmap_table_sig, 2, 
                                           function(x){gsub("TRUE", "*", x)})
    targeted_heatmap_table_sig <- apply(targeted_heatmap_table_sig, 2, 
                                        function(x){gsub("FALSE", "", x)})
    
    # Plot as heatmap
    ph <- pheatmap(targeted_heatmap_table,
                   fontsize = fontsize,
                   color = color,
                   display_numbers = targeted_heatmap_table_sig,
                   fontsize_number = 20,
                   number_color = "black",
                   ...)
    save_pheatmap_pdf(ph,
                      paste0(outdir, file_prefix, "ORA_heatmap_top_gene_sets.pdf"))

    write.table(targeted_heatmap_table_save, file = paste0(outdir, file_prefix, "ORA_heatmap_top_gene_sets_values.txt"),
                sep = "\t", row.names = F, quote = F)
    
  }
  
  # Return table of significant gene sets
  return(merged_ora_table_sig)
  
}






#### Example usage

# Real life AG example here: /Volumes/Partition 2/Core/Cittelly_immune_scRNAseq_Feb2024/Scripts/norm_cluster_celltype.R

### Get data for step 1
#setwd("/Volumes/Partition 2/Reproducible and test scripts/scRNAseq")
#pbmc <- readRDS("scRNAseq_examples/PBMC/Data/Seurat/PBMC_Seurat_filtered_clustered.RDS")

### 1. Generate a list of enriched genes per group
#gene_list <- enriched_genes_by_group(so = pbmc,
#                        group = "seurat_clusters",
#                        outdir = "fam_and_ora/test_output/")

### Get data for step 2
## Need a gene set input as a data frame with two columns: gs_name, gene_symbol
## Below is an example using msigdb but any set of properly formatted genes will work
#library(msigdbr)
#all = msigdbr(species = "human") 
#gene_sets <- all[all$gs_cat %in% "H" | all$gs_subcat %in% "GO:BP",] 
#gene_sets <- gene_sets %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
#background <- rownames(pbmc)
#targeted_gene_sets_heatmap <- c("GOBP_COMPLEMENT_ACTIVATION", 
#                                "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY_INVOLVED_IN_HEART_DEVELOPMENT",
#                                "GOBP_CELL_ACTIVATION")

### 2. ORA on the above list
#ora_table <- ora(gene_sets = gene_sets,
#    gene_list = gene_list,
#    background = background,
#    outdir = "fam_and_ora/test_output/",
#    targeted_gene_sets_heatmap = targeted_gene_sets_heatmap,
#    cluster_cols = TRUE # Shows I can override and add new pheatmap() parameters
#    )


