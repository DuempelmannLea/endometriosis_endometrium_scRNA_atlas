buildReferenceFromSeurat <- function(
    obj, assay = 'RNA', verbose = TRUE, save_umap = TRUE, save_uwot_path = NULL
) {
    if(!assay %in% c('RNA', 'SCT')) {
        stop('Only supported assays are RNA or SCT.')
    }
    res <- list()
    ## TODO: check that these objects are all correctly initialized
    res$Z_corr <- t(obj@reductions$harmony@cell.embeddings)
    res$Z_orig <- t(obj@reductions$pca@cell.embeddings)
    message('Saved embeddings')
    
    res$R <- t(obj@reductions$harmony@misc$R[colnames(obj),]) # seurat does not subset misc slot
    message('Saved soft cluster assignments')
    
    if (assay == 'RNA') {
        vargenes_means_sds <- tibble(
            symbol = obj@assays[[assay]]@var.features, 
            mean = Matrix::rowMeans(obj@assays[[assay]]@data[obj@assays[[assay]]@var.features, ])
        )
        vargenes_means_sds$stddev <- rowSDs(
            obj@assays[[assay]]@data[obj@assays[[assay]]@var.features, ], 
            vargenes_means_sds$mean
        )
    } else if (assay == 'SCT') {
        vargenes_means_sds <- tibble(
            symbol = obj@assays[[assay]]@var.features, 
            mean = Matrix::rowMeans(obj@assays[[assay]]@scale.data[obj@assays[[assay]]@var.features, ])
        )
        asdgc = Matrix(obj@assays[[assay]]@scale.data[obj@assays[[assay]]@var.features, ], sparse = TRUE)
        vargenes_means_sds$stddev <- rowSDs(
            asdgc, 
            vargenes_means_sds$mean
        )
    }
    
    res$vargenes_means_sds <- vargenes_means_sds
    message('Saved variable gene information for ', nrow(vargenes_means_sds), ' genes.')
    
    res$loadings <- obj@reductions$pca@feature.loadings
    message('Saved PCA loadings.')
    
    res$meta_data <- obj@meta.data
    message('Saved metadata.')
    
    ## Check UMAP 
    if (save_umap) {
        if (is.null(save_uwot_path)) {
            error('Please provide a valid path to save_uwot_path in order to save uwot model.')
        }
        if (is.null(obj@reductions$umap@misc$model)) {
            error('uwot model not initialiazed in Seurat object. Please do RunUMAP with umap.method=\'uwot\', return.model=TRUE first.')
        }
        res$umap <- obj@reductions$umap@misc$model
        res$save_uwot_path <- save_uwot_path
        if (file.exists(res$save_uwot_path)) {
            file.remove(res$save_uwot_path)    
        }
        uwot::save_uwot(res$umap, save_uwot_path)
    }
    
    ## Build Reference! 
    if (verbose) 
        message("Calculate final L2 normalized reference centroids (Y_cos)")
    res$centroids = t(cosine_normalize_cpp(res$R %*% t(res$Z_corr), 1))
    if (verbose) 
        message("Calculate reference compression terms (Nr and C)")
    res$cache = compute_ref_cache(res$R, res$Z_corr)
    colnames(res$Z_orig) = row.names(res$metadata)
    rownames(res$Z_orig) = paste0("PC_", seq_len(nrow(res$Z_corr)))
    colnames(res$Z_corr) = row.names(res$metadata)
    rownames(res$Z_corr) = paste0("harmony_", seq_len(nrow(res$Z_corr)))
        
    if (verbose) 
        message("Finished nicely.")
    return(res)    
}

environment(buildReferenceFromSeurat) <- environment(symphony::buildReference)

RunHarmony.Seurat <- function(
  object,
  group.by.vars,
  reduction = 'pca',
  dims.use = NULL,
  theta = NULL,
  lambda = NULL,
  sigma = 0.1,
  nclust = NULL,
  tau = 0,
  block.size = 0.05,
  max.iter.harmony = 10,
  max.iter.cluster = 20,
  epsilon.cluster = 1e-5,
  epsilon.harmony = 1e-4,
  plot_convergence = FALSE,
  verbose = TRUE,
  reference_values = NULL,
  reduction.save = "harmony",
  assay.use = 'RNA',
  project.dim = TRUE,
  ...
) {
  if (reduction == "pca") {
    tryCatch(
      embedding <- Seurat::Embeddings(object, reduction = "pca"),
      error = function(e) {
        if (verbose) {
          message("Harmony needs PCA. Trying to run PCA now.")
        }
        tryCatch(
          object <- Seurat::RunPCA(
            object,
            assay = assay.use, verbose = verbose
          ),
          error = function(e) {
            stop("Harmony needs PCA. Tried to run PCA and failed.")
          }
        )
      }
    )
  } else {
    available.dimreduc <- names(methods::slot(object = object, name = "reductions"))
    if (!(reduction %in% available.dimreduc)) {
      stop("Requested dimension reduction is not present in the Seurat object")
    }
    embedding <- Seurat::Embeddings(object, reduction = reduction)
  }
  if (is.null(dims.use)) {
    dims.use <- seq_len(ncol(embedding))
  }
  dims_avail <- seq_len(ncol(embedding))
  if (!all(dims.use %in% dims_avail)) {
    stop("trying to use more dimensions than computed. Rereun dimension reduction
         with more dimensions or run Harmony with fewer dimensions")
  }
  if (length(dims.use) == 1) {
    stop("only specified one dimension in dims.use")
  }
  metavars_df <- Seurat::FetchData(object, group.by.vars)
    
  harmonyObject <- HarmonyMatrix(
    embedding,
    metavars_df,
    group.by.vars,
    FALSE,
    0,
    theta,
    lambda,
    sigma,
    nclust,
    tau,
    block.size,
    max.iter.harmony,
    max.iter.cluster,
    epsilon.cluster,
    epsilon.harmony,
    plot_convergence,
    TRUE,
    verbose,
    reference_values
  )

  harmonyEmbed <- t(as.matrix(harmonyObject$Z_corr))
  rownames(harmonyEmbed) <- row.names(embedding)
  colnames(harmonyEmbed) <- paste0(reduction.save, "_", seq_len(ncol(harmonyEmbed)))

  harmonyClusters <- t(harmonyObject$R)
  rownames(harmonyClusters) <- row.names(embedding)
  colnames(harmonyClusters) <- paste0('R', seq_len(ncol(harmonyClusters)))
  
  suppressWarnings({
    harmonydata <- Seurat::CreateDimReducObject(
      embeddings = harmonyEmbed,
      stdev = as.numeric(apply(harmonyEmbed, 2, stats::sd)),
      assay = assay.use,
      key = reduction.save,
      misc=list(R=harmonyClusters)
    )
  })

  object[[reduction.save]] <- harmonydata
  if (project.dim) {
    object <- Seurat::ProjectDim(
      object,
      reduction = reduction.save,
      overwrite = TRUE,
      verbose = FALSE
    )
  }
  return(object)
}

environment(RunHarmony.Seurat) <- environment(harmony::HarmonyMatrix)

RunUMAP2 <- function (object, reduction.key = "UMAP_", assay = NULL, reduction.model = NULL, 
    return.model = FALSE, umap.method = "uwot", n.neighbors = 30L, 
    n.components = 2L, metric = "cosine", n.epochs = NULL, learning.rate = 1, 
    min.dist = 0.3, spread = 1, set.op.mix.ratio = 1, local.connectivity = 1L, 
    repulsion.strength = 1, negative.sample.rate = 5, a = NULL, 
    b = NULL, uwot.sgd = FALSE, seed.use = 42, metric.kwds = NULL, 
    angular.rp.forest = FALSE, verbose = TRUE, ...) 
{
    CheckDots(...)
    if (!is.null(x = seed.use)) {
        set.seed(seed = seed.use)
    }
    if (umap.method != "umap-learn" && getOption("Seurat.warn.umap.uwot", 
        TRUE)) {
        warning("The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric", 
            "\nTo use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'", 
            "\nThis message will be shown once per session", 
            call. = FALSE, immediate. = TRUE)
        options(Seurat.warn.umap.uwot = FALSE)
    }
    if (umap.method == "uwot-learn") {
        warning("'uwot-learn' is deprecated. Set umap.method = 'uwot' and return.model = TRUE")
        umap.method <- "uwot"
        return.model <- TRUE
    }
    if (return.model) {
        if (verbose) {
            message("UMAP will return its model")
        }
        umap.method = "uwot"
    }
    if (inherits(x = object, what = "Neighbor")) {
        object <- list(idx = Indices(object), dist = Distances(object))
    }
    if (!is.null(x = reduction.model)) {
        if (verbose) {
            message("Running UMAP projection")
        }
        umap.method <- "uwot-predict"
    }
    umap.output <- switch(EXPR = umap.method, `umap-learn` = {
        if (!py_module_available(module = "umap")) {
            stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
        }
        if (!is.null(x = seed.use)) {
            py_set_seed(seed = seed.use)
        }
        if (typeof(x = n.epochs) == "double") {
            n.epochs <- as.integer(x = n.epochs)
        }
        umap_import <- import(module = "umap", delay_load = TRUE)
        umap <- umap_import$UMAP(n_neighbors = as.integer(x = n.neighbors), 
            n_components = as.integer(x = n.components), metric = metric, 
            n_epochs = n.epochs, learning_rate = learning.rate, 
            min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio, 
            local_connectivity = local.connectivity, repulsion_strength = repulsion.strength, 
            negative_sample_rate = negative.sample.rate, a = a, 
            b = b, metric_kwds = metric.kwds, angular_rp_forest = angular.rp.forest, 
            verbose = verbose)
        umap$fit_transform(as.matrix(x = object))
    }, uwot = {
        if (metric == "correlation") {
            warning("UWOT does not implement the correlation metric, using cosine instead", 
                call. = FALSE, immediate. = TRUE)
            metric <- "cosine"
        }
        if (is.list(x = object)) {
            umap(X = NULL, nn_method = object, n_threads = nbrOfWorkers(), 
                n_components = as.integer(x = n.components), 
                metric = metric, n_epochs = n.epochs, learning_rate = learning.rate, 
                min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio, 
                local_connectivity = local.connectivity, repulsion_strength = repulsion.strength, 
                negative_sample_rate = negative.sample.rate, 
                a = a, b = b, fast_sgd = uwot.sgd, verbose = verbose, 
                ret_model = return.model)
        } else {
            umap(X = object, n_threads = nbrOfWorkers(), n_neighbors = as.integer(x = n.neighbors), 
                n_components = as.integer(x = n.components), 
                metric = metric, n_epochs = n.epochs, learning_rate = learning.rate, 
                min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio, 
                local_connectivity = local.connectivity, repulsion_strength = repulsion.strength, 
                negative_sample_rate = negative.sample.rate, 
                a = a, b = b, fast_sgd = uwot.sgd, verbose = verbose, 
                ret_model = return.model)
        }
    }, `uwot-predict` = {
        if (metric == "correlation") {
            warning("UWOT does not implement the correlation metric, using cosine instead", 
                call. = FALSE, immediate. = TRUE)
            metric <- "cosine"
        }
        if (is.null(x = reduction.model) || !inherits(x = reduction.model, 
            what = "DimReduc")) {
            stop("If running projection UMAP, please pass a DimReduc object with the model stored to reduction.model.", 
                call. = FALSE)
        }
        model <- Misc(object = reduction.model, slot = "model")
        if (length(x = model) == 0) {
            stop("The provided reduction.model does not have a model stored. Please try running umot-learn on the object first", 
                call. = FALSE)
        }
        if (is.list(x = object)) {
            uwot::umap_transform(X = NULL, nn_method = object, 
                model = model, n_threads = nbrOfWorkers(), n_epochs = n.epochs, 
                verbose = verbose)
        } else {
            umap_transform(X = object, model = model, n_threads = nbrOfWorkers(), 
                n_epochs = n.epochs, verbose = verbose)
        }
    }, stop("Unknown umap method: ", umap.method, call. = FALSE))
    if (return.model) {
#         umap.output$nn_index <- NULL
        umap.model <- umap.output
        umap.output <- umap.output$embedding
    }
    colnames(x = umap.output) <- paste0(reduction.key, 1:ncol(x = umap.output))
    if (inherits(x = object, what = "dist")) {
        rownames(x = umap.output) <- attr(x = object, "Labels")
    }
    else if (is.list(x = object)) {
        rownames(x = umap.output) <- rownames(x = object$idx)
    }
    else {
        rownames(x = umap.output) <- rownames(x = object)
    }
    umap.reduction <- CreateDimReducObject(embeddings = umap.output, 
        key = reduction.key, assay = assay, global = TRUE)
    if (return.model) {
        Misc(umap.reduction, slot = "model") <- umap.model
    }
    return(umap.reduction)
}


environment(RunUMAP2) <- environment(Seurat:::RunUMAP.default)

mapQuery <- function (exp_query, metadata_query, ref_obj, vars = NULL, verbose = TRUE, 
    do_normalize = TRUE, do_umap = TRUE, sigma = 0.1, return_type = c('symphony', 'Seurat')) 
{
    if (return_type == 'Seurat') {
        que <- Seurat::CreateSeuratObject(
            counts=exp_query,
            meta.data=metadata_query,
            assay='SymphonyQuery'
        )        
    }
    
    if (do_normalize) {
        if (verbose) 
            message("Normalizing")
        exp_query = normalizeData(exp_query, 10000, "log")
    }
    if (verbose) 
        message("Scaling and synchronizing query gene expression")
    idx_shared_genes = which(ref_obj$vargenes$symbol %in% rownames(exp_query))
    shared_genes = ref_obj$vargenes$symbol[idx_shared_genes]
    if (verbose) 
        message("Found ", length(shared_genes), " reference variable genes in query dataset")
    exp_query_scaled = scaleDataWithStats(exp_query[shared_genes, 
        ], ref_obj$vargenes$mean[idx_shared_genes], ref_obj$vargenes$stddev[idx_shared_genes], 
        1)
    exp_query_scaled_sync = matrix(0, nrow = length(ref_obj$vargenes$symbol), 
        ncol = ncol(exp_query))
    exp_query_scaled_sync[idx_shared_genes, ] = exp_query_scaled
    rownames(exp_query_scaled_sync) = ref_obj$vargenes$symbol
    colnames(exp_query_scaled_sync) = colnames(exp_query)
    if (verbose) 
        message("Project query cells using reference gene loadings")
    Z_pca_query = t(ref_obj$loadings) %*% exp_query_scaled_sync
    if (verbose) 
        message("Clustering query cells to reference centroids")
    Z_pca_query_cos = cosine_normalize_cpp(Z_pca_query, 2)
    R_query = soft_cluster(ref_obj$centroids, Z_pca_query_cos, 
        sigma)
    if (verbose) 
        message("Correcting query batch effects")
    if (!is.null(vars)) {
        design = droplevels(metadata_query)[, vars] %>% as.data.frame()
        onehot = design %>% purrr::map(function(.x) {
            if (length(unique(.x)) == 1) {
                rep(1, length(.x))
            }
            else {
                stats::model.matrix(~0 + .x)
            }
        }) %>% purrr::reduce(cbind)
        Xq = cbind(1, intercept = onehot) %>% t()
    }
    else {
        Xq = Matrix(rbind(rep(1, ncol(Z_pca_query)), rep(1, ncol(Z_pca_query))), 
            sparse = TRUE)
    }
    Zq_corr = moe_correct_ref(as.matrix(Z_pca_query), as.matrix(Xq), 
        as.matrix(R_query), as.matrix(ref_obj$cache[[1]]), as.matrix(ref_obj$cache[[2]]))
    colnames(Z_pca_query) = row.names(metadata_query)
    rownames(Z_pca_query) = paste0("PC_", seq_len(nrow(Zq_corr)))
    colnames(Zq_corr) = row.names(metadata_query)
    rownames(Zq_corr) = paste0("harmony_", seq_len(nrow(Zq_corr)))
    umap_query = NULL
    if (do_umap & !is.null(ref_obj$save_uwot_path)) {
        if (verbose) 
            message("UMAP")
        ref_umap_model = uwot::load_uwot(ref_obj$save_uwot_path, 
            verbose = FALSE)
        
        ## UMAP may have been learned on subset of columns
        umap_query = uwot::umap_transform(t(Zq_corr)[, 1:ref_umap_model$norig_col], ref_umap_model)
#         umap_query = uwot::umap_transform(t(Zq_corr), ref_umap_model)
        colnames(umap_query) = c("UMAP1", "UMAP2")
        rownames(umap_query) <- row.names(metadata_query)
    }
    if (verbose) 
        message("All done!")
    
    if (return_type == 'Seurat') {
        que@assays$SymphonyQuery@data <- exp_query
        que@assays$SymphonyQuery@scale.data <- exp_query_scaled_sync
        que[['pca']] <- Seurat::CreateDimReducObject(
            embeddings = t(Z_pca_query),
            loadings = ref_obj$loadings, 
            stdev = as.numeric(apply(Z_pca_query, 1, stats::sd)),
            assay = 'SymphonyQuery',
            key = 'pca_'
        )
        que[['harmony']] <- Seurat::CreateDimReducObject(
            embeddings = t(Zq_corr),
            stdev = as.numeric(apply(Zq_corr, 1, stats::sd)),
            assay = 'SymphonyQuery',
            key = 'harmony_',
            misc=list(R=R_query)
        )
        que <- Seurat::ProjectDim(que, reduction = 'harmony', overwrite = TRUE, verbose = FALSE)
        if (do_umap) {
            que[['umap']] <- Seurat::CreateDimReducObject(
                embeddings = umap_query,
                assay = 'SymphonyQuery',
                key = 'umap_'
            )            
        }
        return(que)
    } else if (return_type == 'symphony') {
        return(list(Z = Zq_corr, Zq_pca = Z_pca_query, R = R_query, 
            Xq = Xq, umap = umap_query, meta_data = metadata_query))
    } else {
        stop(glue('The return type = \"{return_type}\" is not available.'))
    }
    
}

environment(mapQuery) <- environment(symphony::mapQuery)

knnPredict.Seurat <- function(query_obj, ref_obj, label_transfer, k = 5, confidence = TRUE, seed = 0) 
{
    set.seed(seed)
    if (!label_transfer %in% colnames(ref_obj$meta_data)) {
        stop('Label \"{label_transfer}\" is not available in the reference metadata.')
    }
    
    if (confidence) {
        knn_pred <- class::knn(t(ref_obj$Z_corr), Embeddings(query_obj, 'harmony'), 
            ref_obj$meta_data[[label_transfer]], k = k, prob = TRUE)
        knn_prob = attributes(knn_pred)$prob
        query_obj@meta.data[[label_transfer]] <- knn_pred
        query_obj@meta.data[paste0(label_transfer, '_prob')] = knn_prob
    } else {
        knn_pred <- class::knn(t(ref_obj$Z_corr), Embeddings(query_obj, 'harmony'), 
            ref_obj$meta_data[[label_transfer]], k = k, prob = FALSE)
        query_obj@meta.data[[label_transfer]] <- knn_pred
    }
    return(query_obj)
}

symphony <- function(CELLTYPEref, 
                     CELLTYPEquery, 
                     TanIntegrationVars, 
                     THETA_REF, 
                     ENDOIntegrationVars,
                     dir_out) {
  
  ########################################################
  ##### 00. Create dir_out
  ######################################################## 
  
  #create directory
  dir_out <- paste0("/home/common/data/output/projects/ENDO/E044/A047/prolif_sec_",CELLTYPEref,"ref_",CELLTYPEquery,"query",str_c(ENDOIntegrationVars, collapse = ""),"/")
  dir.create(dir_out, recursive = TRUE)
  setwd(dir_out)
  
  ########################################################
  ##### 01. Load and prepare ref and query datasets
  ########################################################
  
  #Load ref data set
  cells_ref <- readRDS(paste0("/home/common/data/output/projects/ENDO/E039/A001/",CELLTYPEref,".rds"))
  cells_ref@meta.data$dataset <- "ref"
  cells_ref <- NormalizeData(cells_ref, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
  ref_exp_full <- cells_ref@assays$RNA@data
  ref_metadata <- cells_ref@meta.data
  ref_metadata <- ref_metadata[, c("sample","stage","subtypes","sample_type_rename","phase")]
  ref_metadata[,"phase"] <- NA
  ref_metadata <- setnames(ref_metadata, old = c("subtypes", "phase"), new = c("cell_type", "batch"))
  
  #Load query data set
  cells_query <- readRDS(paste0(dir_out,"/prolif_sec_",CELLTYPEquery,".rds"))
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
                   vars = ENDOIntegrationVars,      
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
  
  saveRDS(query$meta_data, paste0("prolif_secretory_",CELLTYPEquery,"_HarmonySymphony.RDS"))
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