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

# input path 
input_path<- '../Tables' 
  
# source common functions 
source('./common_functions.R')

## ---------------------------------- ##
##  Figure 6.a
## ---------------------------------- ##

# @Lea to check 
#Supplementary Data Figure 6a
# reshape df
melted_meta <- melt(meta, id.vars = c("lane", "endo_grade"), measure.vars = c("fold1", "fold2", "fold3", "fold4", "fold5"))
## create the bar plot for endo_EndometriosisGrade
SFig6a <- ggplot(melted_meta, aes(x = variable, fill = endo_grade)) +
  geom_bar(position = "stack", width = 0.7) +
  facet_wrap(~ fct_rev(value), ncol = 2) +
  labs(title = "Bar Plot of Samples by Endo EndometriosisGrade",
       x = "Fold",
       y = "Number of Samples") +
  scale_fill_manual(values = c("severe" = "red", "mild" = "orange", "CTL" = "blue")) +
  theme_minimal() +
  coord_flip()

## ---------------------------------- ##
##  Figure 6.b
## ---------------------------------- ##

# @lea to check 
#Supplementary Data Figure 6b
SFig6b <- ggplot(melted_meta, aes(x = variable, fill = lane)) +
  geom_bar(position = "stack", width = 0.7) +
  facet_wrap(~ fct_rev(value), ncol = 2) +
  labs(title = "Bar Plot of Samples by lane",
       x = "Fold",
       y = "Number of Samples") +
  #scale_fill_manual(values = c("severe" = "red", "mild" = "orange", "CTL" = "blue")) +
  theme_minimal() +
  coord_flip()
SFig6b

## ---------------------------------- ##
##  Figure 6.c
## ---------------------------------- ##

#@lea to check 
#Supplementary Data Figure 6c

metadata <- read.csv('/home/common/data/output/projects/ENDO/E042/A001/summary_plots_cv_splits/metadata_prolif_eval_fold_all.csv')
metadata$new_column <- NULL
metadata$batch <- NULL


# Remove 'names' column as it will be used as row labels
row_names <- metadata$names
metadata <- metadata[, -1]

# Convert non-numeric columns to factors
metadata[, c("EndometriosisStatus", "lane",  "endo_EndometriosisGrade",  "fold1", "fold2", "fold3", "fold4", "fold5")] <- lapply(metadata[, c("EndometriosisStatus", "lane", "endo_EndometriosisGrade", "fold1", "fold2", "fold3", "fold4", "fold5")], as.factor)



# Create heatmap using ComplexHeatmap
pdf(paste0(dir_out, "metadata_prolif_eval_fold_all_heatmap.pdf"), width = 7, height = 9)
Heatmap(metadata,
        name = "Metadata",
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 8),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = TRUE,
        show_row_names = TRUE,
        column_title = "Columns",
        row_title = "Row Names",
        heatmap_legend_param = list(title = "Legend")
)
dev.off()

## ---------------------------------- ##
##  Figure 6.d
## ---------------------------------- ##

output_proba_best_learners_corrgpca <- read.csv(file.path(input_path,'output_probabilities_sample_corrgpca.csv'))
output_proba_best_learners_corrgpca <- output_proba_best_learners_corrgpca %>%
  mutate(grade = recode(grade, "CTL" = "non-ENDO", 
                        "severe" = "severe-ENDO", 
                        "mild" = "mild-ENDO"))
plot_violin(data = output_proba_best_learners_corrgpca,
            y_variable='median_output_probability')
## ---------------------------------- ##
##  Figure 6.e
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
##  Figure 6.f
## ---------------------------------- ##

## Predictive score for each of the genes of the 11 gene signature 
predictive_scores_df <- read.csv(file.path(input_path,'summary_predictive_score_captum_11_genes.csv')) %>%
  dplyr::arrange(predictive_score)

plt <- ggplot(predictive_scores_df) +
  geom_col(aes(predictive_score, fct_reorder(feature,predictive_score)), fill = 'darkblue', width = 0.6)+
  labs(x='Predictive contribution',y='')
plt
