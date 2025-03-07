## ---------------------------------------- ##
## Setup
## ---------------------------------------- ##

# load librairies
library(gprofiler2)
library(plyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(cowplot)
library(ggrepel)

# define paths
dir_out <- "../_Data/05b_Functional_enrichment_analysis/"
setwd(dir_out)

## ---------------------------------------- ##
## Identify inflection point to filter out
## lowly expressed genes 
## ---------------------------------------- ##

##Combine muscat output into one dataframe
PATHS <- list.files(path = paste0("/endometriosis_endometrium_scRNA_atlas/_Data/05a_DEG_analysis_muscat/Proliferative*/", pattern = "mm_dream.rds", recursive = TRUE, full.names = TRUE)
RNAlistRDS <- lapply(PATHS, readRDS)
names(RNAlistRDS) <- print(PATHS)
RNAlistRDS <- unlist(RNAlistRDS, recursive = FALSE) #combine sub-lists
RNAlistRDS <- lapply(RNAlistRDS,plyr::rbind.fill)
RNAlistRDSdf <- as.data.frame(do.call(rbind, RNAlistRDS))

##exclude ribosomal genes (comment: does not change the inflection point)
rb.genes_DEGs <- unique(grep(pattern = "^RP[SL]", RNAlistRDSdf$gene, value=TRUE)) #around 150 unique rb genes
RNAlistRDSdf <- RNAlistRDSdf[!(RNAlistRDSdf$gene %in% rb.genes_DEGs), ]

##Determine inflection point
#Determine Inflection Point, No log2 scale, bw = 0.02
plot(density(RNAlistRDSdf$AveExpr))
# Create the density plot and store it in an object
density_plot <- density(RNAlistRDSdf$AveExpr, bw = 0.02) #good values 0.025 or 0.02, 0.01 is too small
plot(density_plot)
# Calculate the second derivative of the density curve
second_derivative <- diff(diff(density_plot$y)) / diff(density_plot$x)^2
# Find where the second derivative crosses zero (inflection point)
inflection_points <- density_plot$x[which(diff(sign(second_derivative)) != 0)]
# Print or plot the inflection point(s)
cat("Inflection Point(s):", inflection_points, "\n")
# You can also plot the density curve with the inflection point(s) marked
plot(density_plot)
abline(v = inflection_points, col = "red", lty = 2)
abline(v = 18.53202, col = "green", lty = 2)
text(x=18.45, y=1.5, srt=90, '18.53202')

#Best result: no log2 scale, bw = 0.02 with resulting inflection point 18.53202

## ---------------------------------------- ##
## Process deg table 
## ---------------------------------------- ##

# 1. Read in DEG table & filter out lowly expressed and non significant genes 
degs_df <- read.csv(path_degs)
degs_df <- degs_df %>%
  dplyr::filter(p_adj.loc < 0.05) %>%
  dplyr::filter(AveExpr > 18.53202)

##--------------------------------##
## Perform go enrichment analysis 
##--------------------------------##

##1. Perform GO enrichment analysis 
run_go_enrichment <- function(query){
  query <- query %>%
    dplyr::arrange(p_adj.loc)
  gostres <- gost(query = query$gene, 
                  organism = "hsapiens", 
                  ordered_query = TRUE, #can use TRUE, p_adj.loc has been arranged 
                  multi_query = FALSE, 
                  significant = TRUE, #only statistically significant results should be returned
                  exclude_iea = FALSE, 
                  measure_underrepresentation = FALSE, 
                  evcodes = TRUE,
                  user_threshold = 0.05, 
                  correction_method = "g_SCS", #g_SCS is default
                  domain_scope = "annotated", #if custom_bg is set, I believe this shouldn't matter too much
                  custom_bg = NULL, #custom background, use!
                  numeric_ns = "", 
                  sources = NULL, #by default all sources are analyzed
                  as_short_link = FALSE)
  df_output <- gostres$result %>%
    dplyr::select(c('p_value','term_size','query_size','intersection_size','term_name'))
}

go_enrichment_df <- lapply(degs_list,run_go_enrichment(degs_df))
write.csv(go_enrichment_df,file = file.path(dir_out,'GO_enrichment_summary.csv'),quote=F,row.names=F)

##--------------------------------##
## Generate summary plot 
##--------------------------------##

go_enrichment_df %>% 
  dplyr::select(p_value, term_size, query_size, intersection_size, term_id, term_name) %>% 
  dplyr::filter(term_size > 25 & term_size < 750) %>% #larger terms are likely to be excessively broad and prone to false positives
  dplyr::arrange(p_adj.loc) %>% 
  ggplot2::ggplot() + 
  ggplot2::aes(y = reorder(term_name, -log(p_value)), x = -log(p_value), size = intersection_size/term_size * 100) + 
  ggplot2::geom_point() + 
  ggplot2::scale_size_continuous("Intersection size / Term size * 100 (%)") + 
  ggplot2::theme_classic() + 
  ggplot2::xlab("-Log(adj. p-value)") + 
  ggplot2::ylab("enriched terms") + 
  ggplot2::ggtitle("all cells")

ggsave(paste0(dir_out, "GO_enrichment_plot.pdf"),
       plot = last_plot(),
       width = 14, height = 14)

