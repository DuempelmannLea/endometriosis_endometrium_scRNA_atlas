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
library(Seurat)

# input path 
input_path <- '../Tables'

## ---------------------------------- ##
##  Generate summary tables 
## ---------------------------------- ##

read_in_tables <- function(file_path,learner_name='learner_uid'){ 
  data <- read.csv(file_path)
  data %>%
    mutate(cvsplit = gsub(".*(CVsplit\\d+).*", "\\1", !!sym(learner_name)),
           learner_id = gsub(".*(learner\\d+).*", "\\1", !!sym(learner_name))) %>%
    select(!c(!!sym(learner_name)))
}


performance_11_genes <- read_in_tables(file.path(input_path,'Performance_11_gene_signature.csv')) %>%
  dplyr::mutate(dimred = '11_genes_signature',
         cell_selection = 'all_cells',
         Model_training_phase = 'validation')

performance_allcells_signature_confirmation <- read_in_tables(file.path(input_path,'Performance_all_cells_signature_confirmation.csv')) %>%
  dplyr::mutate(dimred = paste0('signature_confirmation_',dimred),
         Model_training_phase = 'validation',
         cell_selection = 'all_cells')

performance_cf_signature_confirmation <- read_in_tables(file.path(input_path,'Performance_cycling_fibroblasts_signature_confirmation.csv')) %>%
  dplyr::mutate(cell_selection = 'cycling_fibroblasts',
         Model_training_phase = 'validation',
         dimred = paste0('signature_confirmation_',dimred))

performance_all_cells_validation <- read.csv(file.path(input_path,'Performance_all_cells_validation.csv')) %>%
  dplyr::mutate(cell_selection = 'all_cells',
         Model_training_phase = 'validation') %>%
  dplyr::rename('learner_id'='learner')

performance_all_cells_evaluation <- read_in_tables(file.path(input_path,'Performance_all_cells_evaluation.csv'),learner_name='learner') %>%
  dplyr::mutate(cell_selection = 'all_cells',
         Model_training_phase = 'evaluation')

performance_cf_evaluation <- read_in_tables(file.path(input_path,'Performance_cycling_fibroblasts_evaluation.csv')) %>%
  dplyr::mutate(cell_selection = 'cycling_fibroblasts',
                Model_training_phase = 'evaluation',
                cvsplit = gsub(".*(cvsplit\\d+).*", "\\1", cvsplit))

performance_cf_validation <- read.csv(file.path(input_path,'Performance_cycling_fibroblasts_validation.csv')) %>%
  dplyr::mutate(cell_selection = 'cycling_fibroblasts',
         Model_training_phase = 'validation') %>%
  dplyr::rename('learner_id'='learner')

combined_df <- bind_rows(performance_all_cells_validation, 
                         performance_all_cells_evaluation, 
                         performance_cf_evaluation,
                         performance_cf_validation,
                         performance_11_genes,
                         performance_allcells_signature_confirmation,
                         performance_cf_signature_confirmation)
write.csv(combined_df,file.path(input_path,'Summary_learners.csv'))

table(combined_df$cvsplit,combined_df$fold)
