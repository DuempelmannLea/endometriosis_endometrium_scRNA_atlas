########################################################
##### 00. Load packages and utils
#######################################################

# load packages
suppressPackageStartupMessages({
  source('/endometriosis_endometrium_scRNA_atlas/02_Cell_annotation/01_libs.R')# imports
}) 

#source utils
source('/endometriosis_endometrium_scRNA_atlas/02_Cell_annotation/01_utils.R') # color definitions, symphony function and plotting functions



########################################################
##### 01. Run DEG analysis with muscat 
#######################################################

##Proliferative Phase
#Proliferative phase samples with strict exclusion criteria, AnnotationMain
DEGanalysis_muscat(MenstrualCyclePhase = "Proliferative", Annotation = "AnnotationMain")
  
#Proliferative phase samples with strict exclusion criteria, AnnotationRefined
DEGanalysis_muscat(MenstrualCyclePhase = "Proliferative", Annotation = "AnnotationRefined")

#Proliferative phase samples with strict exclusion criteria, AnnotationUnited
DEGanalysis_muscat(MenstrualCyclePhase = "Proliferative", Annotation = "AnnotationUnited")

##Secretory Phase
#Secretory phase samples with strict exclusion criteria, AnnotationMain
DEGanalysis_muscat(MenstrualCyclePhase = "Secretory", Annotation = "AnnotationMain")

#Secretory phase samples with strict exclusion criteria, AnnotationRefined
DEGanalysis_muscat(MenstrualCyclePhase = "Secretory", Annotation = "AnnotationRefined")

#Secretory phase samples with strict exclusion criteria, AnnotationUnited
DEGanalysis_muscat(MenstrualCyclePhase = "Secretory", Annotation = "AnnotationUnited")
