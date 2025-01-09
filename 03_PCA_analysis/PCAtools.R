
###############################################
### Load Libraries and define dir_out
###############################################

library(PCAtools)
library(Seurat)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggbreak) 
library(patchwork)
library(ComplexUpset)
library(ggrepel)

dir_out <- ""
dir.create(dir_out, recursive = TRUE)



###############################################
### Prepare Data and run PCA
###############################################

##Aggregate Expression per sample 
data <- readRDS("/home/common/data/output/projects/ENDO/E044/A033/dataA003A033.rds")
bulk <- AggregateExpression(data, group.by = "sample", assay = "integrated", return.seurat = TRUE)

# Extract counts matrix
mat <- bulk[["integrated"]]@data

##Conduct principal component analysis (PCA):
p <- pca(mat, metadata = metadata, removeVar = 0.1)



###############################################
### Plotting with PCAtools
###############################################

##Bi-plots
biplot1 <-    biplot(p,
                     colby = 'FinalPhaseRefined',
                     labSize = 0,drawConnectors = FALSE, 
                     hline = 0, vline = 0,
                     legendPosition = 'right')
ggsave(
  filename = paste0(dir_out,"FinalPhaseRefined.pdf"),
  plot = biplot1,
  width = 7,
  height = 5)

biplot2 <-   biplot(p,
                    colby = 'Wang_epi_pseudot',
                    labSize = 0,drawConnectors = FALSE, 
                    hline = 0, vline = 0,
                    legendPosition = 'right')
ggsave(
  filename = paste0(dir_out,"Wang_epi_pseudot.pdf"),
  plot = biplot2,
  width = 7,
  height = 5)

biplot3 <-   biplot(p,
                    colby = 'treatment',
                    labSize = 0,drawConnectors = FALSE, 
                    hline = 0, vline = 0,
                    legendPosition = 'right')
ggsave(
  filename = paste0(dir_out,"treatment.pdf"),
  plot = biplot3,
  width = 7,
  height = 5)

##eigencor plot
pdf(paste0(dir_out,"eigencorplot.pdf"), width = 10, height = 5)
eigencorplot(p,
             metavars = c("sample","batch", "grade", "treatment", "ProgesteronePeritoneal", "ProgesteroneSerum", "CycleDay", "CyclePhaseMain","CyclePhase","CyclePhaseRefined"),
             col = c('darkorange3','darkorange', 'white','cadetblue2','cadetblue')
)
dev.off()

sessionInfo()
date()