## ---------------------------------- ##
##  Setup
## ---------------------------------- ##

# load librairies
library(ggplot2)
library(rlang)
library(stringr)
library(dplyr)
library(tidyr)
library(cowplot)
library(gridExtra)
library(patchwork)
library(tibble)
library(tinter)
library(UpSetR)
library(forcats)
library(reshape2)
library(ComplexHeatmap)

# input path 
input_path <- '../Tables/' 

#output path
dir_out <- "../_Data/Figures/"

# source common functions 
source('./common_functions.R')

## ---------------------------------- ##
##  Supplementary Data Figure 6.a
## ---------------------------------- ##

#Load data
ScaiVision_meta <- read.csv(file.path(input_path,'ScaiVision_metadata_prolif_folds.csv')) 

# reshape df
melted_meta <- melt(ScaiVision_meta, id.vars = c("lane", "endo_grade"), measure.vars = c("fold1", "fold2", "fold3", "fold4", "fold5"))

# barplot and save
SFig6a <- ggplot(melted_meta, aes(x = variable, fill = endo_grade)) +
  geom_bar(position = "stack", width = 0.7) +
  facet_wrap(~ fct_rev(value), ncol = 2) +
  labs(title = "Bar Plot of Samples by Endometriosis Grade",
       x = "Fold",
       y = "Number of Samples",fill="endometriosis status") +
  scale_fill_manual(values = c("severe" = "red", "mild" = "orange", "CTL" = "blue")) +
  theme_minimal() +
  coord_flip()
ggsave(filename = paste0(dir_out, 'SFig6a.pdf'),
       plot = SFig6a,
       width =12, height=7)

## ---------------------------------- ##
##  Supplementary Data Figure 6.b
## ---------------------------------- ##
#Load data
ScaiVision_meta <- read.csv(file.path(input_path,'ScaiVision_metadata_prolif_folds.csv')) 

# reshape df
melted_meta <- melt(ScaiVision_meta, id.vars = c("lane", "endo_grade"), measure.vars = c("fold1", "fold2", "fold3", "fold4", "fold5"))

# barplot and save
SFig6b <- ggplot(melted_meta, aes(x = variable, fill = lane)) +
  geom_bar(position = "stack", width = 0.7) +
  facet_wrap(~ fct_rev(value), ncol = 2) +
  labs(title = "Bar Plot of Samples by sequencing batch",
       x = "Fold",
       y = "Number of Samples",fill= "Sequencing") +
  scale_fill_manual(labels=c('1st batch','2nd batch'),values=c('darkorange','darkblue')) +
  theme_minimal() +
  coord_flip()
ggsave(filename = paste0(dir_out, 'SFig6b.pdf'),
       plot = SFig6b,
       width =12, height=7)

## ---------------------------------- ##
##  Supplementary Data Figure 6.c
## ---------------------------------- ##

#Load data
ScaiVision_meta <- read.csv(file.path(input_path,'ScaiVision_metadata_prolif_folds.csv')) 

# Remove 'names' column as it will be used as row labels
rownames(ScaiVision_meta) <- ScaiVision_meta$names

# Convert non-numeric columns to factors
ScaiVision_meta[, c("treatment", "lane",  "endo_grade",  "fold1", "fold2", "fold3", "fold4", "fold5")] <- lapply(ScaiVision_meta[, c("treatment", "lane",  "endo_grade",  "fold1", "fold2", "fold3", "fold4", "fold5")], as.factor)

# Prepare metdata for plotting 
ScaiVision_meta_plot <- ScaiVision_meta %>%
  dplyr::select(c( "lane",  "endo_grade",  "fold1", "fold2", "fold3", "fold4", "fold5")) %>%
  dplyr::rename(endometriosis_status ='endo_grade',
               sequencing_batch = 'lane'
                 )
  

# Create heatmap using ComplexHeatmap 
# Colors were modified using Illustrator
pdf(paste0(dir_out, "SFig6c_ScaiVision_meta.pdf"), width = 7, height = 9)
Heatmap(ScaiVision_meta_plot,
        name = "Metadata",
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 8),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = TRUE,
        show_row_names = TRUE,
        column_title = "",
        row_title = "",
        heatmap_legend_param = list(title = "Legend")
)
dev.off()

## ---------------------------------- ##
##  Supplementary Data Figure 6.d
## ---------------------------------- ##

output_proba_best_learners_corrgpca <- read.csv(file.path(input_path,'output_probabilities_sample_corrgpca.csv'))
output_proba_best_learners_corrgpca <- output_proba_best_learners_corrgpca %>%
  mutate(grade = recode(grade, "CTL" = "non-ENDO", 
                        "severe" = "severe-ENDO", 
                        "mild" = "mild-ENDO"))
plot_violin(data = output_proba_best_learners_corrgpca,
            y_variable='median_output_probability')
## ---------------------------------- ##
##  Supplementary Data Figure 6.e
## ---------------------------------- ##

pca_signature <- read.csv(file.path(input_path,'SupplementaryTable6_PCAsignature.csv')) %>%
  dplyr::pull(gene)
corrgpca_signature <- read.csv(file.path(input_path,'SupplementaryTable6_CorrgPCA_signature.csv')) %>%
  dplyr::pull(gene)

list_signatures <- list(CorrgPCA = corrgpca_signature,
                        PCA = pca_signature
                        )
upset(fromList(list_signatures), order.by = "freq",text.scale = 2)

## ---------------------------------- ##
##  Supplementary Data Figure 6.f
## ---------------------------------- ##

## Predictive score for each of the genes of the 11 gene signature 
predictive_scores_df <- read.csv(file.path(input_path,'summary_predictive_score_captum_11_genes.csv')) %>%
  dplyr::arrange(predictive_score)

plt <- ggplot(predictive_scores_df) +
  geom_col(aes(predictive_score, fct_reorder(feature,predictive_score)), fill = 'darkblue', width = 0.6)+
  labs(x='Predictive contribution',y='')
plt
