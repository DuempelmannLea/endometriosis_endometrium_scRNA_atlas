## ----VROOM_CONNECTION_SIZE, eval = TRUE, echo = FALSE-------------------------

Sys.setenv(VROOM_CONNECTION_SIZE='512000')


## ----getPackage, eval = FALSE-------------------------------------------------
#  
#    if (!requireNamespace('BiocManager', quietly = TRUE))
#      install.packages('BiocManager')
#  
#    BiocManager::install('PCAtools')
#  

## ----getPackageDevel, eval = FALSE--------------------------------------------
#  
#    if (!requireNamespace('remotes', quietly = TRUE))
#      install.packages('remotes')
#  
#    remotes::install_github('kevinblighe/PCAtools')
#  

## ----Load, message = FALSE----------------------------------------------------

library(PCAtools)


## ----message = FALSE----------------------------------------------------------

library(airway)
library(magrittr)

data('airway')
airway$dex %<>% relevel('untrt')


## ----message = FALSE----------------------------------------------------------

ens <- rownames(airway)

library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(airway), names(symbols))]
rownames(airway) <- symbols
keep <- !is.na(rownames(airway))
airway <- airway[keep,]


## ----message = FALSE, warning = FALSE-----------------------------------------

library('DESeq2')

dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq(dds)
vst <- assay(vst(dds))


## -----------------------------------------------------------------------------

p <- pca(vst, metadata = colData(airway), removeVar = 0.1)


## ----ex1, warning = FALSE, fig.height = 7, fig.width = 6, fig.cap = 'Figure 1: A scree plot'----

screeplot(p, axisLabSize = 18, titleLabSize = 22)


## ----ex2a, eval = FALSE-------------------------------------------------------
#  
#    biplot(p)
#  

## ----ex2b, fig.height = 7, fig.width = 7, fig.cap = 'Figure 2: A bi-plot'-----

biplot(p, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)


## ----message = FALSE----------------------------------------------------------

library(Biobase)
library(GEOquery)

# load series and platform data from GEO
gset <- getGEO('GSE2990', GSEMatrix = TRUE, getGPL = FALSE)
mat <- exprs(gset[[1]])

# remove Affymetrix control probes
mat <- mat[-grep('^AFFX', rownames(mat)),]

# extract information of interest from the phenotype data (pdata)
idx <- which(colnames(pData(gset[[1]])) %in%
               c('relation', 'age:ch1', 'distant rfs:ch1', 'er:ch1',
                 'ggi:ch1', 'grade:ch1', 'size:ch1',
                 'time rfs:ch1'))
metadata <- data.frame(pData(gset[[1]])[,idx],
                       row.names = rownames(pData(gset[[1]])))

# tidy column names
colnames(metadata) <- c('Study', 'Age', 'Distant.RFS', 'ER', 'GGI', 'Grade',
                        'Size', 'Time.RFS')

# prepare certain phenotypes of interest
metadata$Study <- gsub('Reanalyzed by: ', '', as.character(metadata$Study))
metadata$Age <- as.numeric(gsub('^KJ', NA, as.character(metadata$Age)))
metadata$Distant.RFS <- factor(metadata$Distant.RFS,
                               levels = c(0,1))
metadata$ER <- factor(gsub('\\?', NA, as.character(metadata$ER)),
                      levels = c(0,1))
metadata$ER <- factor(ifelse(metadata$ER == 1, 'ER+', 'ER-'),
                      levels = c('ER-', 'ER+'))
metadata$GGI <- as.numeric(as.character(metadata$GGI))
metadata$Grade <- factor(gsub('\\?', NA, as.character(metadata$Grade)),
                         levels = c(1,2,3))
metadata$Grade <- gsub(1, 'Grade 1', gsub(2, 'Grade 2', gsub(3, 'Grade 3', metadata$Grade)))
metadata$Grade <- factor(metadata$Grade, levels = c('Grade 1', 'Grade 2', 'Grade 3'))
metadata$Size <- as.numeric(as.character(metadata$Size))
metadata$Time.RFS <- as.numeric(gsub('^KJX|^KJ', NA, metadata$Time.RFS))

# remove samples from the pdata that have any NA value
discard <- apply(metadata, 1, function(x) any(is.na(x)))
metadata <- metadata[!discard,]

# filter the expression data to match the samples in our pdata
mat <- mat[,which(colnames(mat) %in% rownames(metadata))]

# check that sample names match exactly between pdata and expression data 
all(colnames(mat) == rownames(metadata))


## -----------------------------------------------------------------------------

p <- pca(mat, metadata = metadata, removeVar = 0.1)


## ----ex3a, eval = FALSE-------------------------------------------------------
#  
#    biplot(p)
#  

## ----ex3b, fig.height = 6, fig.width = 7, fig.cap = 'Figure 3: A bi-plot'-----

biplot(p, showLoadings = TRUE, lab = NULL)


## ----ex4, message = FALSE, fig.height = 10, fig.width = 10, fig.cap = 'Figure 4: A pairs plot'----

pairsplot(p)


## ----ex5, fig.height = 6, fig.width = 8, fig.cap = 'Figure 5: A loadings plot'----

plotloadings(p, labSize = 3)


## ----ex6, warning = FALSE, fig.height = 4, fig.width = 8, fig.cap = 'Figure 6: An eigencor plot'----

eigencorplot(p,
             metavars = c('Study','Age','Distant.RFS','ER',
                          'GGI','Grade','Size','Time.RFS'))


## -----------------------------------------------------------------------------

p$rotated[1:5,1:5]

p$loadings[1:5,1:5]


## ----messages = FALSE---------------------------------------------------------

suppressMessages(require(hgu133a.db))
newnames <- mapIds(hgu133a.db,
                   keys = rownames(p$loadings),
                   column = c('SYMBOL'),
                   keytype = 'PROBEID')

# tidy up for NULL mappings and duplicated gene symbols
newnames <- ifelse(is.na(newnames) | duplicated(newnames),
                   names(newnames), newnames)
rownames(p$loadings) <- newnames


## ----warning = FALSE----------------------------------------------------------

horn <- parallelPCA(mat)
horn$n


## -----------------------------------------------------------------------------

elbow <- findElbowPoint(p$variance)
elbow


## ----ex7, fig.height = 7, fig.width = 9, fig.cap = 'Figure 7: Advanced scree plot illustrating optimum number of PCs'----

library(ggplot2)

screeplot(p,
          components = getComponents(p, 1:20),
          vline = c(horn$n, elbow)) +
  
  geom_label(aes(x = horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8))


## -----------------------------------------------------------------------------

which(cumsum(p$variance) > 80)[1]


## ----ex8, fig.height = 7, fig.width = 7.5, fig.cap = 'Figure 8: Colour by a metadata factor, use a custom label, add lines through origin, and add legend'----

biplot(p,
       lab = paste0(p$metadata$Age, ' aÃ±os'),
       colby = 'ER',
       hline = 0, vline = 0,
       legendPosition = 'right')


## ----ex9, message = FALSE, fig.height = 7.5, fig.width = 7, fig.cap = 'Figure 9: Supply custom colours and encircle variables by group'----

biplot(p,
       colby = 'ER', colkey = c('ER+' = 'forestgreen', 'ER-' = 'purple'),
       colLegendTitle = 'ER-\nstatus',
       # encircle config
       encircle = TRUE,
       encircleFill = TRUE,
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)

biplot(p,
       colby = 'ER', colkey = c('ER+' = 'forestgreen', 'ER-' = 'purple'),
       colLegendTitle = 'ER-\nstatus',
       # encircle config
       encircle = TRUE, encircleFill = FALSE,
       encircleAlpha = 1, encircleLineSize = 5,
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)


## ----ex10, message = FALSE, fig.height = 7.5, fig.width = 7, fig.cap = 'Figure 10: Stat ellipses'----

biplot(p,
       colby = 'ER', colkey = c('ER+' = 'forestgreen', 'ER-' = 'purple'),
       # ellipse config
       ellipse = TRUE,
       ellipseType = 't',
       ellipseLevel = 0.95,
       ellipseFill = TRUE,
       ellipseAlpha = 1/4,
       ellipseLineSize = 1.0,
       xlim = c(-125,125), ylim = c(-50, 80),
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)

biplot(p,
       colby = 'ER', colkey = c('ER+' = 'forestgreen', 'ER-' = 'purple'),
       # ellipse config
       ellipse = TRUE,
       ellipseType = 't',
       ellipseLevel = 0.95,
       ellipseFill = TRUE,
       ellipseAlpha = 1/4,
       ellipseLineSize = 0,
       ellipseFillKey = c('ER+' = 'yellow', 'ER-' = 'pink'),
       xlim = c(-125,125), ylim = c(-50, 80),
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)


## ----ex11, message = FALSE, eval = FALSE--------------------------------------
#  
#    biplot(p,
#      colby = 'ER', colkey = c('ER+' = 'forestgreen', 'ER-' = 'purple'),
#      hline = c(-25, 0, 25), vline = c(-25, 0, 25),
#      legendPosition = 'top', legendLabSize = 13, legendIconSize = 8.0,
#      shape = 'Grade', shapekey = c('Grade 1' = 15, 'Grade 2' = 17, 'Grade 3' = 8),
#      drawConnectors = FALSE,
#      title = 'PCA bi-plot',
#      subtitle = 'PC1 versus PC2',
#      caption = '27 PCs â‰ˆ 80%')
#  

## ----ex12a--------------------------------------------------------------------

biplot(p,
       lab = NULL,
       colby = 'ER', colkey = c('ER+'='royalblue', 'ER-'='red3'),
       hline = c(-25, 0, 25), vline = c(-25, 0, 25),
       vlineType = c('dotdash', 'solid', 'dashed'),
       gridlines.major = FALSE, gridlines.minor = FALSE,
       pointSize = 5,
       legendPosition = 'left', legendLabSize = 14, legendIconSize = 8.0,
       shape = 'Grade', shapekey = c('Grade 1'=15, 'Grade 2'=17, 'Grade 3'=8),
       drawConnectors = FALSE,
       title = 'PCA bi-plot',
       subtitle = 'PC1 versus PC2',
       caption = '27 PCs â‰ˆ 80%')


## ----ex12b, message = FALSE, fig.height = 5.5, fig.width = 11, fig.cap = 'Figure 11: Modify line types, remove gridlines, and increase point size'----

biplot(p,
       # loadings parameters
       showLoadings = TRUE,
       lengthLoadingsArrowsFactor = 1.5,
       sizeLoadingsNames = 4,
       colLoadingsNames = 'red4',
       # other parameters
       lab = NULL,
       colby = 'ER', colkey = c('ER+'='royalblue', 'ER-'='red3'),
       hline = 0, vline = c(-25, 0, 25),
       vlineType = c('dotdash', 'solid', 'dashed'),
       gridlines.major = FALSE, gridlines.minor = FALSE,
       pointSize = 5,
       legendPosition = 'left', legendLabSize = 14, legendIconSize = 8.0,
       shape = 'Grade', shapekey = c('Grade 1'=15, 'Grade 2'=17, 'Grade 3'=8),
       drawConnectors = FALSE,
       title = 'PCA bi-plot',
       subtitle = 'PC1 versus PC2',
       caption = '27 PCs â‰ˆ 80%')


## ----ex13a, fig.height = 5.5, fig.width = 6, fig.cap = 'Figure 12: Colour by a continuous variable and plot other PCs'----

# add ESR1 gene expression to the metadata
p$metadata$ESR1 <- mat['205225_at',]

biplot(p,
       x = 'PC2', y = 'PC3',
       lab = NULL,
       colby = 'ESR1',
       shape = 'ER',
       hline = 0, vline = 0,
       legendPosition = 'right') +
  
  scale_colour_gradient(low = 'gold', high = 'red2')


## ----ex13b, eval = FALSE------------------------------------------------------
#  
#    # was always eval = FALSE
#    biplot(p, x = 'PC10', y = 'PC50',
#      lab = NULL,
#      colby = 'Age',
#      hline = 0, vline = 0,
#      hlineWidth = 1.0, vlineWidth = 1.0,
#      gridlines.major = FALSE, gridlines.minor = TRUE,
#      pointSize = 5,
#      legendPosition = 'left', legendLabSize = 16, legendIconSize = 8.0,
#      shape = 'Grade', shapekey = c('Grade 1'=15, 'Grade 2'=17, 'Grade 3'=8),
#      drawConnectors = FALSE,
#      title = 'PCA bi-plot',
#      subtitle = 'PC10 versus PC50',
#      caption = '27 PCs â‰ˆ 80%')
#  

## ----ex14, message = FALSE, fig.height = 8, fig.width = 7, fig.cap = 'Figure 13: Quickly explore potentially informative PCs via a pairs plot'----

pairsplot(p,
          components = getComponents(p, c(1:10)),
          triangle = TRUE, trianglelabSize = 12,
          hline = 0, vline = 0,
          pointSize = 0.4,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'Grade',
          title = 'Pairs plot', plotaxes = FALSE,
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))


## ----ex15, fig.height = 5.5, fig.width = 6, fig.cap = 'Figure 14: arranging a pairs plot horizontally'----

pairsplot(p,
          components = getComponents(p, c(4,33,11,1)),
          triangle = FALSE,
          hline = 0, vline = 0,
          pointSize = 0.8,
          gridlines.major = FALSE, gridlines.minor = FALSE,
          colby = 'ER',
          title = 'Pairs plot', titleLabSize = 22,
          axisLabSize = 14, plotaxes = TRUE,
          margingaps = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))


## ----ex16, fig.height = 5.5, fig.width = 7, fig.cap = 'Figure 15: Determine the variables that drive variation among each PC'----

plotloadings(p,
             rangeRetain = 0.01,
             labSize = 4.0,
             title = 'Loadings plot',
             subtitle = 'PC1, PC2, PC3, PC4, PC5',
             caption = 'Top 1% variables',
             shape = 24,
             col = c('limegreen', 'black', 'red3'),
             drawConnectors = TRUE)


## ----ex17a, fig.height = 7, fig.width = 9, fig.cap = 'Figure 16: plotting absolute component loadings'----

plotloadings(p,
             components = getComponents(p, c(4,33,11,1)),
             rangeRetain = 0.1,
             labSize = 4.0,
             absolute = FALSE,
             title = 'Loadings plot',
             subtitle = 'Misc PCs',
             caption = 'Top 10% variables',
             shape = 23, shapeSizeRange = c(1, 16),
             col = c('white', 'pink'),
             drawConnectors = FALSE)


## ----ex17b, fig.height = 4.5, fig.width = 9, fig.cap = 'Figure 17: plotting absolute component loadings'----

plotloadings(p,
             components = getComponents(p, c(2)),
             rangeRetain = 0.12, absolute = TRUE,
             col = c('black', 'pink', 'red4'),
             drawConnectors = TRUE, labSize = 4) + coord_flip()


## ----ex18a, warning = FALSE, fig.height = 4.25, fig.width = 12, fig.cap = 'Figure 18: Correlate the principal components back to the clinical data'----

eigencorplot(p,
             components = getComponents(p, 1:27),
             metavars = c('Study','Age','Distant.RFS','ER',
                          'GGI','Grade','Size','Time.RFS'),
             col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
             cexCorval = 0.7,
             colCorval = 'white',
             fontCorval = 2,
             posLab = 'bottomleft',
             rotLabX = 45,
             posColKey = 'top',
             cexLabColKey = 1.5,
             scale = TRUE,
             main = 'PC1-27 clinical correlations',
             colFrame = 'white',
             plotRsquared = FALSE)


## ----ex18b, warning = FALSE, fig.height = 5, fig.width = 12, fig.cap = 'Figure 19: Correlate the principal components back to the clinical data'----

eigencorplot(p,
             components = getComponents(p, 1:horn$n),
             metavars = c('Study','Age','Distant.RFS','ER','GGI',
                          'Grade','Size','Time.RFS'),
             col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
             cexCorval = 1.2,
             fontCorval = 2,
             posLab = 'all',
             rotLabX = 45,
             scale = TRUE,
             main = bquote(Principal ~ component ~ Pearson ~ r^2 ~ clinical ~ correlates),
             plotRsquared = TRUE,
             corFUN = 'pearson',
             corUSE = 'pairwise.complete.obs',
             corMultipleTestCorrection = 'BH',
             signifSymbols = c('****', '***', '**', '*', ''),
             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))


## ----ex19, message = FALSE, warning = FALSE, fig.height = 10, fig.width = 15, fig.cap = 'Figure 20: a merged panel of all PCAtools plots'----

pscree <- screeplot(p, components = getComponents(p, 1:30),
                    hline = 80, vline = 27, axisLabSize = 14, titleLabSize = 20,
                    returnPlot = FALSE) +
  geom_label(aes(20, 80, label = '80% explained variation', vjust = -1, size = 8))

ppairs <- pairsplot(p, components = getComponents(p, c(1:3)),
                    triangle = TRUE, trianglelabSize = 12,
                    hline = 0, vline = 0,
                    pointSize = 0.8, gridlines.major = FALSE, gridlines.minor = FALSE,
                    colby = 'Grade',
                    title = '', plotaxes = FALSE,
                    margingaps = unit(c(0.01, 0.01, 0.01, 0.01), 'cm'),
                    returnPlot = FALSE)

pbiplot <- biplot(p,
                  # loadings parameters
                  showLoadings = TRUE,
                  lengthLoadingsArrowsFactor = 1.5,
                  sizeLoadingsNames = 4,
                  colLoadingsNames = 'red4',
                  # other parameters
                  lab = NULL,
                  colby = 'ER', colkey = c('ER+'='royalblue', 'ER-'='red3'),
                  hline = 0, vline = c(-25, 0, 25),
                  vlineType = c('dotdash', 'solid', 'dashed'),
                  gridlines.major = FALSE, gridlines.minor = FALSE,
                  pointSize = 5,
                  legendPosition = 'none', legendLabSize = 16, legendIconSize = 8.0,
                  shape = 'Grade', shapekey = c('Grade 1'=15, 'Grade 2'=17, 'Grade 3'=8),
                  drawConnectors = FALSE,
                  title = 'PCA bi-plot',
                  subtitle = 'PC1 versus PC2',
                  caption = '27 PCs â‰ˆ 80%',
                  returnPlot = FALSE)

ploadings <- plotloadings(p, rangeRetain = 0.01, labSize = 4,
                          title = 'Loadings plot', axisLabSize = 12,
                          subtitle = 'PC1, PC2, PC3, PC4, PC5',
                          caption = 'Top 1% variables',
                          shape = 24, shapeSizeRange = c(4, 8),
                          col = c('limegreen', 'black', 'red3'),
                          legendPosition = 'none',
                          drawConnectors = FALSE,
                          returnPlot = FALSE)

peigencor <- eigencorplot(p,
                          components = getComponents(p, 1:10),
                          metavars = c('Study','Age','Distant.RFS','ER',
                                       'GGI','Grade','Size','Time.RFS'),
                          cexCorval = 1.0,
                          fontCorval = 2,
                          posLab = 'all', 
                          rotLabX = 45,
                          scale = TRUE,
                          main = "PC clinical correlates",
                          cexMain = 1.5,
                          plotRsquared = FALSE,
                          corFUN = 'pearson',
                          corUSE = 'pairwise.complete.obs',
                          signifSymbols = c('****', '***', '**', '*', ''),
                          signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                          returnPlot = FALSE)

library(cowplot)
library(ggplotify)

top_row <- plot_grid(pscree, ppairs, pbiplot,
                     ncol = 3,
                     labels = c('A', 'B  Pairs plot', 'C'),
                     label_fontfamily = 'serif',
                     label_fontface = 'bold',
                     label_size = 22,
                     align = 'h',
                     rel_widths = c(1.10, 0.80, 1.10))

bottom_row <- plot_grid(ploadings,
                        as.grob(peigencor),
                        ncol = 2,
                        labels = c('D', 'E'),
                        label_fontfamily = 'serif',
                        label_fontface = 'bold',
                        label_size = 22,
                        align = 'h',
                        rel_widths = c(0.8, 1.2))

plot_grid(top_row, bottom_row, ncol = 1,
          rel_heights = c(1.1, 0.9))


## -----------------------------------------------------------------------------

p <- pca(mat, metadata = metadata, removeVar = 0.1)
p.prcomp <- list(sdev = p$sdev,
                 rotation = data.matrix(p$loadings),
                 x = data.matrix(p$rotated),
                 center = TRUE, scale = TRUE)

class(p.prcomp) <- 'prcomp'

# for this simple example, just use a chunk of
# the original data for the prediction
newdata <- t(mat[,seq(1,20)])
predict(p.prcomp, newdata = newdata)[,1:5]


## -----------------------------------------------------------------------------

sessionInfo()

