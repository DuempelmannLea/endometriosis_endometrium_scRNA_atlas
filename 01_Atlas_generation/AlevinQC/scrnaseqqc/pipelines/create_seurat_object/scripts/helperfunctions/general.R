# Helper Functions for core pipe

###########################################################
# Function for Normalisation, HVG selection and scaling
###########################################################
# I: data, snakemake config table
# O: normalised data as seurat object
normalisation_scaling <- function(data, config ) {
  switch((snakemake@config$regr.vars == "NA") + 1, 
         regr.vars <- strsplit(snakemake@config$regr.vars, ",")[[1]], 
         regr.vars <- NULL )
  # load the vars to regress for
  if (paste0(config$norm.method) == c("standart")) {
    n.features <- config$nfeatures
    # data normalisation
    ###########################################################
    # data normalization by the gene expression measurements for each cell by the total expression (lib size normalisation), log transformation
    data <- NormalizeData(object = data, normalization.method = "LogNormalize", 
                          scale.factor = config$scale.factor)
    # Variable feature selection
    ###########################################################
    # Detection of variable genes across the single cells by genewise mean and variance, fitting a model to the relationship of log(variance) and log(mean) using local polynomial regression (loess).The feature values are standartized by the observed mean and expected variance (given by the fitted line).# Feature variance is then calculated on the standardized values.
    data <- FindVariableFeatures(object = data, selection.method = "vst",
                                 nfeatures = n.features )
    # Scale Data
    ###########################################################
    # Scaling the data and remove sources of unwanted variation
    data <- ScaleData(object =   data,
                      features = VariableFeatures(object = data),
                      vars.to.regress = regr.vars,
                      do.scale = TRUE,
                      do.center = TRUE)
    DefaultAssay(data) <- "RNA"
    
  } else if ( paste0(config$norm.method) == "sctransform") { 
    # data normalization by SCtansform using residuals, log transformation
    data <- SCTransform(object = data,  vars.to.regress = regr.vars,
                        conserve.memory=FALSE,
                        seed.use = 123,
                        return.only.var.genes=TRUE)
    data  <-  ScaleData(data, assay = "RNA")
    DefaultAssay(data) <- "SCT"
    
  } else {
    stop("no valid normalisation method selected") 
  }
  return(data)
}


######################################
# Import Alevin counts
######################################
# Parts of the function is taken from Seurat's Read10x parsing function. 
# https://combine-lab.github.io/alevin-tutorial/2018/alevin-seurat/
# I: file paths to count matrix, gene and barcode names
# O: Count matrix in sparse matrix format

ReadAlevin <- function(base.path){
  if (!dir.exists(base.path)) {
    stop("Directory provided does not exist")
  }
  barcode.loc <- file.path(base.path, "quants_mat_rows.txt")
  gene.loc <- file.path(base.path, "quants_mat_cols.txt")
  matrix.loc <- file.path(base.path,  "quants_mat.gz")
  if (!file.exists(barcode.loc)) {
    stop("Barcode file missing")
  }
  if (!file.exists(gene.loc)) {
    stop("Gene name file missing")
  }
  if (!file.exists(matrix.loc)) {
    stop("Expression matrix file missing")
  }
  cell.names <- readLines(con = barcode.loc)
  gene.names <- readLines(con = gene.loc)
  num.cells <- length(x = cell.names)
  num.genes <- length(x = gene.names)
  out.matrix <- matrix(data = NA, nrow = num.genes, ncol = num.cells)
  con <- gzcon(con = file(description = matrix.loc, open = "rb"))
  total.molecules <- 0
  for (n in seq_len(length.out = num.cells)) {
    try(out.matrix[,n ] <- readBin(con = con, what = double(), 
                               endian = "little", n = num.genes),
        silent = TRUE)

  }
  colnames(x = out.matrix) <- cell.names
  rownames(x = out.matrix) <- gene.names
  return(out.matrix)
}


######################################
# rename the gene names to symbol 
######################################
# I: Seurat object , snakemake config file, lookup table containg gene entrez ids
# O: Seurat object with gene names in data and rawdata slots changed to external gene names
rename_features <- function(matrix, lookup, config){
  m <- matrix
  # query
  ## acces mart
  #mart  <-  useMart(biomart = "ensembl", dataset = paste0(config$ensemble_database),  host="uswest.ensembl.org", ensemblRedirect = FALSE) 
  # format ensemble ids
  x <- data.frame( gene_list= rownames(m) ) 
  x$gene_list <- as.character( x$gene_list )
  
  x <- x %>%  
    dplyr::mutate(orig_name = sapply(strsplit(gene_list, split= ".", fixed =TRUE) , "[[", 1) )
  stopifnot(length(x$orig_name) == length(rownames(m)) )
  # check if gene_list is equal to rownames
  stopifnot( (rownames(m) == x$gene_list) == TRUE)
  # query
  # parameters for serach query
  ifelse(config$gene_format == "ensemble_gene_id", filter <-"ensembl_gene_id", filter <- "hgnc_symbol" )
  # create index
  y <- lookup
  # match by indexing for every obj sepereately
  i <- cbind( match(rownames(m), x$gene_list ),
              match(x$orig_name, y$ensembl_gene_id ) )
  
  rownames(m) <- y$external_gene_name[i[,2]]
  # index vector of duplicates
  dup.ind <- which(duplicated(rownames(m)) | duplicated(rownames(m)[length(rownames(m)):1 ])[length(rownames(m)):1])
  # change genenames of duplicates back to ensemble id 
  rownames(m)[dup.ind] <- x$gene_list[dup.ind]
  #check uniqeness
  stopifnot(length(rownames(m)) == length(unique(rownames((m)))))
  return(m)
}
######################################
# Create lookup table between genenames
######################################
# I: Seurat object , snakemake config file
# O: table containing names for ensembl ID, entrez gene, hgnc symbol and gene biotype.
entrez_lookup <- function(matrix, config){
  m <- matrix
  ## acces mart
  hosts = c("useast.ensembl.org", "uswest.ensembl.org", "asia.ensembl.org")
  for (h in 1:length(hosts)) {
    skip_to_next <- FALSE
    #cat("trying to acces:", hosts[h])
    mart = tryCatch({
      useMart(biomart = "ensembl", dataset = paste0(config$ensemble_database), host=hosts[h])
      }, error = function(e) { 
        skip_to_next <<- TRUE
      })
    
    if(skip_to_next) { 
      next 
      } else {
      break
      }
  }
    # format ensemble ids
  x <- data.frame( gene_list= rownames(m) ) 
  
  x$gene_list <- as.character( x$gene_list )
  
  x <- x %>%  
    dplyr::mutate(orig_name = sapply(strsplit(gene_list, split= ".", fixed =TRUE) , "[[", 1) )
  
  # check if gene_list is equal to rownames
  stopifnot( (rownames(m) == x$gene_list) == TRUE)
  # query
  # parameters for serach query
  ifelse(config$gene_format == "ensemble_gene_id", filter <-"ensembl_gene_id", filter <- "hgnc_symbol" )
  # create index
  y <- getBM(
    mart = mart,
    useCache = FALSE,
    values = x$orig_name,
    filter = filter,
    attributes = c("ensembl_gene_id" ,"entrezgene_id", "external_gene_name",  "external_gene_source" ,"gene_biotype")
  )
  # what are the ensembl_gene_id returning NAs in entrez
  unique.ensembles <-  y[is.na(y$entrezgene),]
  
  return(y)
}

