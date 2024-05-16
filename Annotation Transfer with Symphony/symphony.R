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
  CELLTYPEquery = "ENDO_global",
  TanIntegrationVars = c("sample","stage"),
  THETA_REF = c(2,1),
  ENDOIntegrationVars = c("batch", "sample")
)

#01 lymphocyte
symphony(
  CELLTYPEref = "Tan_lymphocyte",
  CELLTYPEquery = "ENDO_lymphocyte",
  TanIntegrationVars = c("sample","stage"),
  THETA_REF = c(2,1),
  ENDOIntegrationVars = c("batch", "sample")
)

#02 myeloid
symphony(
  CELLTYPEref = "Tan_myeloid",
  CELLTYPEquery = "ENDO_myeloid",
  TanIntegrationVars = c("sample","stage"),
  THETA_REF = c(2,1),
  ENDOIntegrationVars = c("batch", "sample")
)

#03 epithelial
symphony(
  CELLTYPEref = "Tan_epithelial",
  CELLTYPEquery = "ENDO_epithelial",
  TanIntegrationVars = c("sample"),
  THETA_REF = c(1.5),
  ENDOIntegrationVars = c("batch", "sample")
)

#04 endothelial
symphony(
  CELLTYPEref = "Tan_endothelial",
  CELLTYPEquery = "ENDO_endothelial",
  TanIntegrationVars = c("sample","stage"),
  THETA_REF = c(2,1),
  ENDOIntegrationVars = c("sample","batch")
)

#05 stromal fibroblasts
symphony(
  CELLTYPEref = "Tan_stromal",
  CELLTYPEquery = "ENDO_stromal",
  TanIntegrationVars = c("sample","stage"),
  THETA_REF = c(2,1),
  ENDOIntegrationVars = c("batch", "sample")
)
