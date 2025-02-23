plot_scatterplot <- function(data,y_variable,output_path){
  dimred <- data %>% pull(dimred) %>% unique
  data <- data %>% 
    select(c('auc','fold')) %>%
    mutate(
           auc = as.numeric(auc))
  
  data.summary <- aggregate(. ~ fold, median, data=data)
  p <- ggplot(data, aes(x=fold, y= !!rlang::sym(y_variable)))
  p <- p + geom_jitter(width=0.2) 
  p <- p + geom_crossbar(data=data.summary, aes(ymin = !!rlang::sym(y_variable), ymax =!!rlang::sym(y_variable)),
                         size=0.2,col="red", width = .2)
  p <- p + labs(x='',shape='Folds')
  p <- p + ylim(c(0,1.2))
  p <- p + annotate(geom="text",
                    x=2, y=0.05, 
                    label= paste0("Median = ",round(median(data[,y_variable]),digits = 2)),
                    size=3)

  p
}
plot_boxplot <- function(data,y_variable,output_path){
  dimred <- data %>% pull(dimred) %>% unique
  data <- data %>% select(c('auc','fold'))
  data.summary <- aggregate(. ~ fold, median, data=data)
  p <- ggplot(data, aes(x=fold, y=!!rlang::sym(y_variable)))
  p <- p + geom_boxplot()
  p <- p + labs(x='')
  p <- p + ylim(c(0,1.2))
  p <- p + ggplot2::geom_hline(yintercept=0.78, linetype = "dotted", col = "salmon")
  p <- p + annotate(geom="text",
                    x=2, y=0.05, 
                    label= paste0("Median = ",round(median(data[,y_variable]),digits = 2)),
                    size=3)
  p
}
plot_boxplot_scailyte_custom <- function(data,accross_dimreds=T,learner_column='learner_uid',yintercept_value=0.5){
  data <- data %>% dplyr::mutate(cvsplit = str_extract(!!sym(learner_column), "(?i)CVsplit\\d+"))
  plot <- data %>%
    ggplot2::ggplot() + 
    ggplot2::aes(x = cvsplit, y = as.numeric(auc), color = cvsplit, fill = cvsplit) + 
    ggplot2::geom_boxplot(alpha = .8, outlier.shape = NA) + 
    ggplot2::geom_jitter(alpha = .2) +
    ggplot2::theme_classic() + 
    ggplot2::ylim(0,1.2) + 
    ggplot2::scale_fill_manual(breaks = c("CVsplit1", "CVsplit2", "CVsplit3"), values = c("#7fc7d3", "#ec6e33", "#076382")) + 
    ggplot2::scale_color_manual(breaks = c("CVsplit1", "CVsplit2", "CVsplit3"), values = darken(c("#7fc7d3", "#ec6e33", "#076382"), amount = 0.2)) + 
    ggplot2::xlab("") + 
    ggplot2::xlab("") +
    ggplot2::geom_hline(yintercept=yintercept_value, linetype = "dotted", col = "salmon")+
    ggplot2::ylab("Confirmation AUC") + 
    ggplot2::theme(strip.background = element_blank(),
                   strip.text = element_text(size = 12),
                   axis.text = element_text(size = 10),
                   legend.position = "none") 
  if (accross_dimreds){
    median_data <- data %>%
      group_by(dimred) %>%
      summarise(median_auc = round(median(auc), digits = 2))
    plot <- plot + ggplot2::facet_wrap(~dimred) 
    plot +  ggplot2::geom_text(data = median_data, aes(x=2, y=0.05, 
                                                                    label = paste0("Median = ", median_auc)), 
                                            size = 3, inherit.aes = FALSE)

  } else {
    plot +  annotate(geom="text",
                     x=2, y=0.05, 
                     label= paste0("Median = ",round(median(data[,'auc']),digits = 2)),
                     size=3)
  
  }
}
go_enrichment_plot <- function(go,top_gene_threshold=50){
  go$result %>% 
    dplyr::select(p_value, term_size, query_size, intersection_size, term_id, term_name) %>% 
    dplyr::filter(term_size > 10 & term_size < 500) %>% #larger terms are likely to be excessively broad and prone to false positives
    dplyr::arrange(p_value) %>% 
    filter(grepl(pattern = "GO:", x = term_id)) %>%
    dplyr::mutate(padj = p.adjust(p = p_value, method = "fdr")) %>% 
    dplyr::filter(padj < 0.05) %>% 
    dplyr::arrange(padj) %>%  
    dplyr::slice_head(n = top_gene_threshold) %>% 
    ggplot2::ggplot() + 
    ggplot2::aes(y = reorder(term_name, -log(padj)), x = -log(padj), color = -log(padj), size = intersection_size/term_size * 100) + 
    ggplot2::geom_point() + 
    ggplot2::scale_color_continuous(name = "-Log(adj. p-value)") + 
    ggplot2::scale_size_continuous("Intersection (%)") + 
    ggplot2::theme_classic() + 
    ggplot2::xlab("-Log(adjusted. p-value)") + 
    ggplot2::ylab("") + 
    ggplot2::ggtitle("GO enrichment")
}
plot_violin <- function(data,y_variable){
  data$grade <- factor(data$grade,levels = c('severe-ENDO','mild-ENDO','non-ENDO'))
  plot <- ggplot(data,aes(x=grade,y=!!rlang::sym(y_variable),fill=grade))
  plot <- plot + geom_violin() 
  plot <- plot + scale_fill_manual(values=c("severe-ENDO" = "#dc7443",
                                            "mild-ENDO" = "#E68E30",
                                            "non-ENDO" = "#7FC7D3"),
                                   breaks=c("severe-ENDO", "mild-ENDO", "non-ENDO"))
  plot <- custom_violin(plot)
  plot <- plot + ylab('Median ENDO probability accross folds')
  plot
}
