

#####################################
# Trajectory analysis with monocle3
# 
# srun --pty -c20 -t5-0:00:00 --mem=200000 /bin/bash
# singularity exec --bind /data/projects/p735_Endometriosis /mnt/apps/centos7/R_v4.3.1_scRNAseq.sif R
#####################################

library(Seurat)          
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(ggplot2)


workPath <- "/data/projects/p735_Endometriosis"            ### ADJUST
analysisFolder <- paste0(workPath, "/downstream_analyses/trajectoryAnalysis/monocle3") # This will be correct, if files were produced by scRNAseq_qc-filtering-normalization.rmd

out <- "E044A050_seed621_filteredWO_ciliatedMesothelialMUC5B"
rawData <- paste0(workPath, "/raw_data/Heidi_timetrajectories/E044A050_seed621_filteredWO_ciliatedMesothelialMUC5B.rds")
metadataFile <- ""

subsetClusters <- c()  
#-------------

dir.create(paste0(analysisFolder, "/", out))

data <- readRDS(rawData)

if(metadataFile != ""){
  metadata <- readRDS(metadataFile)
  data@meta.data$Symphony_Refined <- metadata[rownames(data@meta.data), "Symphony_Refined"]
  data@meta.data$CyclePhase <- metadata[rownames(data@meta.data), "CyclePhase"]
  data@meta.data$CycleDay <- metadata[rownames(data@meta.data), "CycleDay"]
}

DefaultAssay(data) <- "SCT"
Idents(data) <- data@meta.data$Symphony_Refined

pdf(paste0(analysisFolder, "/", out, "/UMAP_Symphony_Refined.pdf"))
  DimPlot(data, reduction = "umap", label=TRUE, group.by="Symphony_Refined") + labs(title="Symphony_Refined")
dev.off()

pdf(paste0(analysisFolder, "/", out, "/UMAP_CycleDay.pdf"))
  DimPlot(data, reduction = "umap", label=TRUE, group.by="CycleDay") + labs(title="CycleDay")
dev.off()

pdf(paste0(analysisFolder, "/", out, "/UMAP_CyclePhase.pdf"))
  DimPlot(data, reduction = "umap", label=TRUE, group.by="CyclePhase") + labs(title="CyclePhase")
dev.off()

if(is.null(subsetClusters)){
  subsetCellnames <- WhichCells(data, idents = subsetClusters)
  data <- subset(data, cells=subsetCellnames)
}

data_monocl <- as.cell_data_set(data)

# retrive clustering information from Seurat object
recreate.partitions <- c(rep(1, length(data_monocl@colData@rownames)))
names(recreate.partitions) <- data_monocl@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
data_monocl@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
list.cluster <- data@active.ident
data_monocl@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
#assign UMAP coordinates
data_monocl@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- data@reductions$umap@cell.embeddings



# trajectory analysis ##############################
out2 <- paste0(out, "/trajectory_day5-6-7")
dir.create(paste0(analysisFolder, "/", out2))
data_monocl <- learn_graph(data_monocl)

pdf(paste0(analysisFolder, "/", out2, "/trajectory_clusters.pdf"), width = 10, height = 10)
plot_cells(data_monocl, color_cells_by="cluster",
           label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE,
           group_label_size=4, cell_size=1,trajectory_graph_segment_size=2,graph_label_size=3)
dev.off()

pdf(paste0(analysisFolder, "/", out2, "/trajectory_CycleDay.pdf"), width = 10, height = 10)
plot_cells(data_monocl, color_cells_by="CycleDay",
           label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE,
           group_label_size=4, cell_size=1,trajectory_graph_segment_size=2,graph_label_size=3)
dev.off()

pdf(paste0(analysisFolder, "/",out2, "/trajectory_CyclePhase.pdf"), width = 10, height = 10)
plot_cells(data_monocl, color_cells_by="CyclePhase",
           label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE,
           group_label_size=4, cell_size=1,trajectory_graph_segment_size=2,graph_label_size=3)
dev.off()

saveRDS(data_monocl, paste0(analysisFolder, "/", out2, "/trajectory.rds"))

#Set timepoint 0:
data_monocl <- order_cells(data_monocl)
out2 <- paste0(out, "/trajectory_day5-6-7_greenCircle")
dir.create(paste0(analysisFolder, "/", out2))
saveRDS(data_monocl, paste0(analysisFolder, "/", out2, "/trajectory_timepoint0.rds"))

pdf(paste0(analysisFolder, "/", out2, "/trajectory_pseudotime.pdf"), width = 10, height = 10)
  plot_cells(data_monocl, color_cells_by = "pseudotime",
                     label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE,
                     graph_label_size=3, cell_size = 1)
dev.off()



# Differential expression along trajectory --------
# Perform the test
data_monocl <- estimate_size_factors(data_monocl)
pr_graph_test_res <- graph_test(data_monocl, neighbor_graph="principal_graph", cores = 17)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

rowData(data_monocl)$gene_name <- rownames(data_monocl)
rowData(data_monocl)$gene_short_name <- rowData(data_monocl)$gene_name

interesting_genes <- pr_graph_test_res %>% 
  filter(q_value < 0.05) %>% 
  arrange(desc(morans_I)) 

write.table(interesting_genes, paste0(analysisFolder, "/", out2, "/trajectory_pseudotime_diffExpression.txt"), sep="\t", quote=FALSE, row.names = TRUE)

pdf(paste0(analysisFolder, "/", out2, "/trajectory_diffExpression_diffExpression.pdf"), width = 10, height = 10)
plot_cells(data_monocl, genes=rownames(interesting_genes)[1:9],
           show_trajectory_graph=TRUE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, cell_size=1)
dev.off()

saveRDS(pr_graph_test_res, paste0(analysisFolder, "/", out2, "/trajectory_diffExpression_graphTestRes.rds"))
saveRDS(data_monocl, paste0(analysisFolder, "/", out2, "/trajectory_diffExpression.rds"))



##Grouping these genes into modules can reveal fate specific genes or those that are activate immediate prior to or following the branch point:
#collect the trajectory-variable genes into modules

data_monocl <- preprocess_cds(data_monocl)

gene_module_df <- find_gene_modules(data_monocl[rownames(interesting_genes),], resolution=c(10^seq(-6,-1)))
write.table(gene_module_df, paste0(analysisFolder, "/", out2, "/modules.txt"), sep="\t", quote=FALSE, row.names = TRUE)
saveRDS(gene_module_df, paste0(analysisFolder, "/", out2, "/modules.rds"))

start <- 1
while(start <= length(table(gene_module_df$module))){
  if(start+8 <= length(table(gene_module_df$module))){
    pdf(paste0(analysisFolder, "/", out2, "/modules", start, "-", start+8, ".pdf"), width = 10, height = 10)
  } else {
    if(start == length(table(gene_module_df$module))){
      start <- start - 1
    }
    pdf(paste0(analysisFolder, "/", out2, "/modules", start, "-", length(table(gene_module_df$module)), ".pdf"), width = 10, height = 10)
  }
  print(plot_cells(data_monocl,
                   genes=gene_module_df %>% filter(module %in% c(start:(start+8))),
                   label_cell_groups=FALSE,
                   show_trajectory_graph=TRUE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE))
  dev.off()
  start <- start + 9
}

# associate ciclyDay with module
cluster_CycleDay <- data.frame(cell = row.names(colData(data_monocl)),
                             CycleDay_group = colData(data_monocl)$CycleDay)
agg_mat <- aggregate_gene_expression(data_monocl, gene_module_df, cluster_CycleDay)
row.names(agg_mat) <- paste0("Module ", row.names(agg_mat))

pdf(paste0(analysisFolder, "/", out2, "/modules-CycleDay_heatmap.pdf"), width = 10, height = 10)
pheatmap::pheatmap(agg_mat,
                   scale="column",
                   cluster_cols=FALSE,
                   clustering_method="ward.D2")
dev.off()

# associate ciclyPhase with module
cluster_CyclePhase <- data.frame(cell = row.names(colData(data_monocl)),
                               CycleDay_group = colData(data_monocl)$CyclePhase)
agg_mat <- aggregate_gene_expression(data_monocl, gene_module_df, cluster_CyclePhase)
row.names(agg_mat) <- paste0("Module ", row.names(agg_mat))

pdf(paste0(analysisFolder, "/", out2, "/modules-CyclePhase_heatmap.pdf"), width = 10, height = 10)
pheatmap::pheatmap(agg_mat,
                   scale="column",
                   cluster_cols=FALSE,
                   clustering_method="ward.D2")
dev.off()


for(module in unique(gene_module_df$module)){
  genes <- gene_module_df[gene_module_df$module == module, ]
  
  pdf(paste0(analysisFolder, "/", out2, "/module", module, "_violinPlot_top6.pdf"), width = 10, height = 10)
  print(plot_genes_violin(data_monocl[rowData(data_monocl)$gene_short_name %in% genes$id[1:6],], group_cells_by="CycleDay", ncol=2) +
          theme(axis.text.x=element_text(angle=45, hjust=1)))
  dev.off()
  
  pdf(paste0(analysisFolder, "/", out2, "/module", module, "_ExpressionPlot_top4_CycleDay.pdf"), width = 10, height = 10)
  print(plot_genes_in_pseudotime(data_monocl[genes$id[1:4],],
                           color_cells_by="CycleDay",
                           min_expr=0.5))
  dev.off()
  
  pdf(paste0(analysisFolder, "/", out2, "/module", module, "_ExpressionPlot_top4_CyclePhase.pdf"), width = 10, height = 10)
  print(plot_genes_in_pseudotime(data_monocl[genes$id[1:4],],
                                 color_cells_by="CyclePhase",
                                 min_expr=0.5))
  dev.off()
  
  pdf(paste0(analysisFolder, "/", out2, "/module", module, "_heatmap_scaledGenes.pdf"), width = 10, height = 10)
   try({print(DoHeatmap(object = data, features = genes$id, group.by="CycleDay"))})
  dev.off()
  #GSEA and enrich GO not possible due to gene annotation mix
}

# differential expression along trajectory separate endo/control ##########
data_monocl.ctl <- data_monocl[, rownames(data_monocl@colData[data_monocl@colData$treatment == "CTL",])]
data_monocl.endo <- data_monocl[, rownames(data_monocl@colData[data_monocl@colData$treatment == "Endo",])]

pdf(paste0(analysisFolder, "/", out2, "/trajectory_pseudotime_ctl.pdf"), width = 10, height = 10)
plot_cells(data_monocl.ctl, color_cells_by = "pseudotime",
           label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE,
           graph_label_size=3, cell_size = 1)
dev.off()

pdf(paste0(analysisFolder, "/", out2, "/trajectory_pseudotime_endo.pdf"), width = 10, height = 10)
plot_cells(data_monocl.endo, color_cells_by = "pseudotime",
           label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE,
           graph_label_size=3, cell_size = 1)
dev.off()


pr_graph_test_res.ctl <- graph_test(data_monocl.ctl, neighbor_graph="principal_graph", cores = 17)
pr_deg_ids.ctl <- row.names(subset(pr_graph_test_res.ctl, q_value < 0.05))

rowData(data_monocl.ctl)$gene_name <- rownames(data_monocl.ctl)
rowData(data_monocl.ctl)$gene_short_name <- rowData(data_monocl.ctl)$gene_name

interesting_genes.ctl <- pr_graph_test_res.ctl %>% 
  filter(q_value < 0.05) %>% 
  arrange(desc(morans_I)) 

write.table(interesting_genes.ctl, paste0(analysisFolder, "/", out2, "/trajectory_pseudotime_diffExpression_ctl.txt"), sep="\t", quote=FALSE, row.names = TRUE)

pr_graph_test_res.endo <- graph_test(data_monocl.endo, neighbor_graph="principal_graph", cores = 17)
pr_deg_ids.endo <- row.names(subset(pr_graph_test_res.endo, q_value < 0.05))

rowData(data_monocl.endo)$gene_name <- rownames(data_monocl.endo)
rowData(data_monocl.endo)$gene_short_name <- rowData(data_monocl.endo)$gene_name

interesting_genes.endo <- pr_graph_test_res.endo %>% 
  filter(q_value < 0.05) %>% 
  arrange(desc(morans_I)) 

write.table(interesting_genes.endo, paste0(analysisFolder, "/", out2, "/trajectory_pseudotime_diffExpression_endo.txt"), sep="\t", quote=FALSE, row.names = TRUE)

#in endo, but not ctl
interesting_genes.endo.notCtl <- interesting_genes.endo[!rownames(interesting_genes.endo) %in% rownames(interesting_genes.ctl),]
write.table(interesting_genes.endo.notCtl, paste0(analysisFolder, "/", out2, "/trajectory_pseudotime_diffExpression_endo.notCtl.txt"), sep="\t", quote=FALSE, row.names = TRUE)

pdf(paste0(analysisFolder, "/", out2, "/trajectory_diffExpression_diffExpression.endo.notCtl_endo.pdf"), width = 10, height = 10)
plot_cells(data_monocl.endo, genes=rownames(interesting_genes.endo.notCtl)[1:9],
           show_trajectory_graph=TRUE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, cell_size=1)
dev.off()

pdf(paste0(analysisFolder, "/", out2, "/trajectory_diffExpression_diffExpression.endo.notCtl_ctl.pdf"), width = 10, height = 10)
plot_cells(data_monocl.ctl, genes=rownames(interesting_genes.endo.notCtl)[1:9],
           show_trajectory_graph=TRUE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, cell_size=1)
dev.off()

#in ctl, but not endo
interesting_genes.ctl.notEndo <- interesting_genes.ctl[!rownames(interesting_genes.ctl) %in% rownames(interesting_genes.endo),]
write.table(interesting_genes.ctl.notEndo, paste0(analysisFolder, "/", out2, "/trajectory_pseudotime_diffExpression_ctl.notEndo.txt"), sep="\t", quote=FALSE, row.names = TRUE)

pdf(paste0(analysisFolder, "/", out2, "/trajectory_diffExpression_diffExpression.ctl.notEndo_endo.pdf"), width = 10, height = 10)
plot_cells(data_monocl.endo, genes=rownames(interesting_genes.ctl.notEndo)[1:9],
           show_trajectory_graph=TRUE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, cell_size=1)
dev.off()

pdf(paste0(analysisFolder, "/", out2, "/trajectory_diffExpression_diffExpression.ctl.notEndo_ctl.pdf"), width = 10, height = 10)
plot_cells(data_monocl.ctl, genes=rownames(interesting_genes.ctl.notEndo)[1:9],
           show_trajectory_graph=TRUE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, cell_size=1)
dev.off()



# use Monocle to find genes that are differentially expressed according xxx ###########################
regressionVariable <- "CycleDay"
regressionVariable <- "CyclePhase"

dir.create(paste0(analysisFolder, "/", out, "/diffExpression_", regressionVariable))

gene_fits <- fit_models(data_monocl, model_formula_str = paste0("~", regressionVariable))
fit_coefs <- coefficient_table(gene_fits)
terms <- fit_coefs %>% filter(term == regressionVariable)

#get genes that have a significant time component
interesting_genes <- terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate) %>% 
  arrange(q_value) 

write.table(interesting_genes, paste0(analysisFolder, "/", out, "/diffExpression_", regressionVariable, "/diffExpression_", regressionVariable, ".txt"), sep="\t", quote=FALSE, row.names = TRUE)

#saveRDS(gene_fits, paste0(analysisFolder, "/", out, "/diffExpression_", regressionVariable, "/diffExpression_", regressionVariable, "_geneFits.rds"))


#heatmap of top 50 scaled genes
selectedGenes <- interesting_genes$gene_short_name[interesting_genes$gene_short_name %in% rownames(GetAssayData(object = data, slot = "scale.data"))]


pdf(paste0(analysisFolder, "/", out, "/diffExpression_", regressionVariable, "/diffExpression_", regressionVariable, "_heatmap_top50.pdf"))
  DoHeatmap(object = data, features = selectedGenes[1:50], group.by=regressionVariable )
dev.off()

##Grouping these genes into modules can reveal fate specific genes 
#collect the CycleDay-variable genes into modules
gene_module_df <- find_gene_modules(data_monocl[interesting_genes$gene_short_name,], resolution=c(10^seq(-6,-1)))
write.table(gene_module_df, paste0(analysisFolder, "/", out, "/diffExpression_", regressionVariable, "/modules.txt"), sep="\t", quote=FALSE, row.names = TRUE)
saveRDS(gene_module_df, paste0(analysisFolder, "/", out, "/diffExpression_", regressionVariable, "/modules.rds"))

start <- 1
while(start <= length(table(gene_module_df$module))){
  if(start+8 <= length(table(gene_module_df$module))){
    pdf(paste0(analysisFolder, "/", out, "/diffExpression_", regressionVariable, "/modules", start, "-", start+8, ".pdf"), width = 10, height = 10)
  } else {
    pdf(paste0(analysisFolder, "/", out, "/diffExpression_", regressionVariable, "/modules", start, "-", length(table(gene_module_df$module)), ".pdf"), width = 10, height = 10)
  }
  if(start == length(table(gene_module_df$module))){
    start <- start - 1
  }
  print(plot_cells(data_monocl,
             genes=gene_module_df %>% filter(module %in% c(start:(start+8))),
             label_cell_groups=FALSE,
             show_trajectory_graph=FALSE))
  dev.off()
  start <- start + 9
}


# associate xxx with module
cluster <- data.frame(cell = row.names(colData(data_monocl)),
                               CycleDay_group = colData(data_monocl)[, regressionVariable])
agg_mat <- aggregate_gene_expression(data_monocl, gene_module_df, cluster)
row.names(agg_mat) <- paste0("Module ", row.names(agg_mat))

pdf(paste0(analysisFolder, "/", out, "/diffExpression_", regressionVariable, "/modules-", regressionVariable, "_heatmap.pdf"), width = 10, height = 10)
pheatmap::pheatmap(agg_mat,
                   scale="column",
                   cluster_cols=FALSE,
                   clustering_method="ward.D2")

dev.off()

for(module in unique(gene_module_df$module)){
  genes <- gene_module_df[gene_module_df$module == module, ]

  pdf(paste0(analysisFolder, "/", out, "/diffExpression_", regressionVariable, "/module", module, "_violinPlot_top6.pdf"), width = 10, height = 10)
  print(plot_genes_violin(data_monocl[rowData(data_monocl)$gene_short_name %in% genes$id[1:6],], group_cells_by=regressionVariable, ncol=2) +
    theme(axis.text.x=element_text(angle=45, hjust=1)))
  dev.off()
  
  #heatmap of top 50 genes
  pdf(paste0(analysisFolder, "/", out, "/diffExpression_", regressionVariable, "/module", module, "_heatmap_scaledGenes.pdf"))
  try({print(DoHeatmap(object = data, features = genes$id, group.by=regressionVariable ))})
  dev.off()
}


get_citations(data_monocl)

