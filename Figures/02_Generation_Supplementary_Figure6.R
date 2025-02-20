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

## to do lea 

## ---------------------------------- ##
##  Figure 6.b
## ---------------------------------- ##

## to do lea 

## ---------------------------------- ##
##  Figure 6.c
## ---------------------------------- ##

## to do lea 

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
