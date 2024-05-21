########################################################
##### 00. Load packages and utils
#######################################################

# load packages
suppressPackageStartupMessages({
  source('./libs.R')# imports
}) 

#source utils
source('./utils.R') # color definitions, symphony function and plotting functions

########################################################
##### 01. Run Symphony
#######################################################

#00 global
symphony(
  CELLTYPEref = "Tan_global",
  CELLTYPEquery = "ENDO_global"
)

#01 lymphocyte
symphony(
  CELLTYPEref = "Tan_lymphocyte",
  CELLTYPEquery = "ENDO_lymphocyte"
)

#02 myeloid
symphony(
  CELLTYPEref = "Tan_myeloid",
  CELLTYPEquery = "ENDO_myeloid"
)

#03 epithelial
symphony(
  CELLTYPEref = "Tan_epithelial",
  CELLTYPEquery = "ENDO_epithelial"
)

#04 endothelial
symphony(
  CELLTYPEref = "Tan_endothelial",
  CELLTYPEquery = "ENDO_endothelial"
)

#05 stromal fibroblasts
symphony(
  CELLTYPEref = "Tan_stromal",
  CELLTYPEquery = "ENDO_stromal"
)
