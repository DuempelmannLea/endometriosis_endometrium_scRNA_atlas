fig.size <- function (height, width) {
  options(repr.plot.height = height, repr.plot.width = width)
}

plotBasic = function(umap_labels,                # metadata, with UMAP labels in UMAP1 and UMAP2 slots
                     title = 'Query',         # Plot title
                     color.by = 'cell_type',  # metadata column name for coloring
                     facet.by = NULL,         # (optional) metadata column name for faceting
                     color.mapping = NULL,    # custom color mapping
                     legend.position = 'right') {  # Show cell type legend
  
  p = umap_labels %>%
    dplyr::sample_frac(1L) %>% # permute rows randomly
    ggplot(aes(x = UMAP1, y = UMAP2)) + 
    geom_point_rast(aes(col = get(color.by)), size = 0.3, stroke = 0.2, shape = 16)
  if (!is.null(color.mapping)) { p = p + scale_color_manual(values = color.mapping) }
  
  # Default formatting
  p = p + theme_bw() +
    labs(title = title, color = color.by) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position=legend.position) +
    theme(legend.text = element_text(size=8), legend.title=element_text(size=12)) + 
    guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = 'none')
  
  if(!is.null(facet.by)) {
    p = p + facet_wrap(~get(facet.by)) +
      theme(strip.text.x = element_text(size = 12)) }    
  return(p)
}

symphony <- function(CELLTYPEref, 
                     CELLTYPEquery) {
  
  ########################################################
  ##### 00. Create dir_out
  ######################################################## 
  
  #create directory
  dir_out <- paste0("../_Data/02_Cell_annotation/",CELLTYPEref,"ref_",CELLTYPEquery,"query","/")
  dir.create(dir_out, recursive = TRUE)
  setwd(dir_out)
  
  ########################################################
  ##### 01. Load and prepare ref and query datasets
  ########################################################
  
  #Load ref data set
  cells_ref <- readRDS(paste0("../_Data/02_Cell_annotation/Tan_",CELLTYPEref,".rds"))
  cells_ref@meta.data$dataset <- "ref"
  cells_ref <- NormalizeData(cells_ref, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
  ref_exp_full <- cells_ref@assays$RNA@data
  ref_metadata <- cells_ref@meta.data
  ref_metadata <- ref_metadata[, c("sample","stage","subtypes","sample_type_rename","phase")]
  ref_metadata[,"phase"] <- NA
  ref_metadata <- setnames(ref_metadata, old = c("subtypes", "phase"), new = c("cell_type", "batch"))
  
  #Load query data set
  cells_query <- readRDS(paste0(dir_out,"/ENDO_",CELLTYPEquery,".rds"))
  DefaultAssay(cells_query) <- "RNA"
  cells_query <- DietSeurat(cells_query, assay="RNA", dimreducs = "umap")
  cells_query@meta.data$dataset <- "query"
  cells_query <- NormalizeData(cells_query, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
  query_exp <- cells_query@assays$RNA@data
  query_metadata <- cells_query@meta.data[,c("sample","stage","cell_type","sample_type_rename","batch")]
  query_metadata$sample <- as.factor(query_metadata$sample)
  
  ########################################################
  ##### 02. Build Symphony Reference
  ########################################################
  ##Option 1: Build from Harmony object (preferred method)
  # Sparse matrix with the normalized genes x cells matrix
  ref_exp_full[1:5, 1:2]
  
  #Select variable genes ofand subset reference expression by variable genes
  var_genes <- row.names(cells_ref@assays$RNA@meta.features)[cells_ref@assays$RNA@meta.features$highly_variable]
  ref_exp = ref_exp_full[var_genes, ]
  dim(ref_exp)
  
  #Calculate and save the mean and standard deviations for each gene
  vargenes_means_sds = tibble(symbol = var_genes, mean = Matrix::rowMeans(ref_exp))
  vargenes_means_sds$stddev = singlecellmethods::rowSDs(ref_exp, vargenes_means_sds$mean)
  head(vargenes_means_sds)
  
  #Scale data using calculated gene means and standard deviations
  ref_exp_scaled = singlecellmethods::scaleDataWithStats(ref_exp, vargenes_means_sds$mean, vargenes_means_sds$stddev, 1)
  
  #Run SVD, save gene loadings (s$u)
  set.seed(1)
  s = irlba(ref_exp_scaled, nv = 20) #can result in error due to NAs
  Z_pca_ref = diag(s$d) %*% t(s$v) # [pcs by cells]
  loadings = s$u
  
  #Run Harmony integration
  set.seed(1)
  if (CELLTYPEref == "Tan_epithelial") {
  TanIntegrationVars <- "sample"
  THETA_REF <- 1.5
  } else if (CELLTYPEref %in% c("Tan_global", "Tan_lymphocyte", "Tan_myeloid", "Tan_endothelial", "Tan_stromal")) {
  TanIntegrationVars <- c("sample", "stage")
  THETA_REF <- c(2, 1)
  }   #Same parameters as in Tan et al. 2022 paper
  ref_harmObj = harmony::HarmonyMatrix(
    data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
    meta_data = ref_metadata, ## dataframe with cell labels
    npcs = 50,                ## RunHarmony on Seurat object uses all (50), not default HarmonyMatrix 20
    theta = THETA_REF,        ## cluster diversity enforcement
    sigma = 0.1,
    vars_use = TanIntegrationVars,   ## variable to integrate out
    nclust = 100,             ## number of clusters in Harmony model
    max.iter.harmony = 20,
    return_object = TRUE,     ## return the full Harmony model object
    do_pca = FALSE,           ## don't recompute PCs
    plot_convergence = TRUE   ## plot convergence
  )
  ggsave(filename = paste0(dir_out,'Haromy_convergence_sample.png'))
  
  # Compress a Harmony object into a Symphony reference
  reference = symphony::buildReferenceFromHarmonyObj(
    ref_harmObj,            # output object from HarmonyMatrix()
    ref_metadata,           # reference cell metadata
    vargenes_means_sds,     # gene names, means, and std devs for scaling
    loadings,               # genes x PCs matrix
    verbose = TRUE,         # verbose output
    do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
    save_uwot_path = './testing_uwot_model_1')
  
  # Optionally, you can specify which normalization method was
  # used to build the reference as a custom slot inside the Symphony object to 
  # help record this information for future query users
  reference$normalization_method = 'log(CP10k+1)'
  
  #Save Symphony reference (modify with your desired output path)
  #saveRDS(reference, './reference.rds')
  str(reference)
  dim(reference$Z_corr)
  reference$Z_corr[1:5, 1:5]
  
  #Visualize reference UMAP
  umap_labels_ref = cbind(ref_metadata, reference$umap$embedding)
  fig.size(3, 5)
  plotBasic(umap_labels_ref, title = 'Reference', color.by = 'cell_type')
  
  ########################################################
  ##### 03. Map query
  #########################################################
  
  #The query dataset is assumed to have been normalized in the same manner as the reference cells (here, default is log(CP10k+1) normalization).
  # Read in Symphony reference to map to
  #reference = readRDS('./testing_reference1.rds')
  # Map query
  query = mapQuery(query_exp,             # query gene expression (genes x cells)
                   query_metadata,        # query metadata (cells x attributes)
                   reference,             # Symphony reference object
                   vars = c("batch", "sample"),      
                   sigma = 0.1,
                   do_normalize = FALSE,  # perform log(CP10k+1) normalization on query
                   do_umap = TRUE)        # project query cells into reference UMAP
  str(query)
  
  #Predict query cell types using k-NN
  query = knnPredict(query, reference, reference$meta_data$cell_type, k = 5)
  #saveRDS(query, './query.rds')
  
  #Query cell type predictions are now in the cell_type_pred_knn column. The cell_type_pred_knn_prob column reports the proportion of nearest neighbors with the winning vote (can help identify query cells that fall "on the border" between 2 reference cell types).
  head(query$meta_data)
  
  #Visualize query UMAP
  umap_labels_query = cbind(query$meta_data, query$umap)
  fig.size(3, 5)
  plotBasic(umap_labels_query, title = 'Reference', color.by = 'cell_type_pred_knn')
  
  saveRDS(query$meta_data, paste0("prolif_secretory_",CELLTYPEquery,"_Symphony.RDS"))
  ########################################################
  ##### 04. Visualization of mapping
  #########################################################
  # Sync the column names for both data frames
  reference$meta_data$cell_type_pred_knn = NA
  reference$meta_data$cell_type_pred_knn_prob = NA
  reference$meta_data$ref_query = 'reference'
  query$meta_data$ref_query = 'query'
  
  # Add the UMAP coordinates to the metadata
  meta_data_combined = rbind(query$meta_data, reference$meta_data)
  umap_combined = rbind(query$umap, reference$umap$embedding)
  umap_combined_labels = cbind(meta_data_combined, umap_combined)
  
  # Plot UMAP visualization of all cells
  p0 <- plotBasic(umap_combined_labels, title = paste0(CELLTYPEref," reference and ",CELLTYPEquery," query"), color.by = 'ref_query')
  p1 <- plotBasic(umap_labels_ref, title = 'Reference', color.by = 'cell_type')
  p2 <- plotBasic(umap_labels_query, title = 'Query', color.by = 'cell_type_pred_knn')
  p3 <- plotBasic(umap_labels_query, title = 'Query', color.by = 'cell_type_pred_knn_prob')
  p4 <- plotBasic(umap_labels_ref, title = 'Reference', color.by = 'sample_type_rename')
  ggsave(plot_grid(p0,p1,p2,p3,p4, ncol = 2, align = "v"), filename = paste0(dir_out,'summary_p5_2.png'),
         width = 15, height = 17)
  
  plot(density(query$meta_data$cell_type_pred_knn_prob), main = paste0("Density ", CELLTYPEref," reference and ",CELLTYPEquery," query"), )
  ggsave(filename = paste0(dir_out,'celltypepredknnprob_density.png'))
  
  fig.size(3, 7)
  plotBasic(umap_combined_labels, title = paste0(CELLTYPEref," reference and ",CELLTYPEquery," query"), 
            color.by = 'cell_type', facet.by = 'ref_query')
  ggsave(filename = paste0(dir_out,'Reference and query cells_facet.png'))
  plotBasic(umap_combined_labels, title = paste0(CELLTYPEref," reference and ",CELLTYPEquery," query"),
            color.by = 'sample_type_rename', facet.by = 'ref_query')
  ggsave(filename = paste0(dir_out,'Reference_query_sampletype.png'))
  plotBasic(umap_combined_labels, title = paste0(CELLTYPEref," reference and ",CELLTYPEquery," query"), 
            color.by = 'cell_type_pred_knn', facet.by = 'ref_query')
  ggsave(filename = paste0(dir_out,'Reference_query_celltypepredknn.png'))
  plotBasic(umap_combined_labels, title = paste0(CELLTYPEref," reference and ",CELLTYPEquery," query"), 
            color.by = 'cell_type_pred_knn_prob', facet.by = 'ref_query')
  ggsave(filename = paste0(dir_out,'Reference_query_celltypepredknnprob.png'))
  
  rm(list=ls())
  writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
}
