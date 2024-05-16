####################################
#### Load libs and set paths
####################################

library(CellChat)
library(paFibroblastshwork)
options(stringsAsFactors = FALSE)

dir_out <- "A076/"


####################################
#### Prepare combine cellchat object for Endo and Control
####################################

#desired label order in plots
labels.levels <- c("VSMC","Prv-STEAP4","Prv-CCL19","Prv-MYH11","fib C7","fib C7-SFRP2","eF4-CXCL14", 
                   "dS1-myofibroblast", "dS2","eF1","eF3","eF2", 
                   "mesothelial","glandular", "TP63+/KRT5+", "lumenal", "lumenal 1",  "lumenal 2","MUC5B+", "ciliated", 
                   "EC-aPCV", "EC-PCV","EC-capillary","EC-tip", "EC-HEV", "EC-artery", "LEC", 
                   "monocytes-CD16+", "monocytes-CD16-", "mast cells","pDC", "mDC","cDC1","pre-cDC2","cDC2","DC3", 
                   "M$\\Phi$1-LYVE1", "M$\\Phi$2-peritoneal", "M$\\Phi$3-APOE", "M$\\Phi$4-infiltrated", "M$\\Phi$5-activated",
                   "B cell", "plasma",  "pNK","NK1","NK2", "NK3", "T$_{Reg}$","T$_N$/T$_{CM}$", "CD4 T$_{RM}$","CD8 T$_{RM}$",  "T$_{EM}$", "CTL",  "ILC")

#Load objects and combine
cellchat.Endo <- readRDS(paste0(dir_out,"Endo_CellChat_Symphony_Refined_Final.rds"))
cellchat.Endo <- updateClusterLabels(cellchat.Endo, new.order = labels.levels)
table(cellchat.Endo@meta[["labels"]])
cellchat.Cntrl <- readRDS(paste0(dir_out,"Control_CellChat_Symphony_Refined_Final.rds"))
cellchat.Cntrl <- updateClusterLabels(cellchat.Cntrl, new.order = labels.levels)
table(cellchat.Cntrl@meta[["labels"]])
object.list <- list(Endo = cellchat.Endo, Cntrl = cellchat.Cntrl)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
saveRDS(cellchat, paste0(dir_out,"cellchat_CntrlvsEndo.rds"))


####################################
#### Plot overall interactions 
####################################

#Comparison of interactions Endo vs Control
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'output_comparison_Refined/compareInteractions.pdf'), width =5,height=5)

#netVisual_heatmap
gd1 <- netVisual_heatmap(cellchat)
gd2 <- netVisual_heatmap(cellchat, measure = "weight")
gd1 + gd2
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'output_comparison_Refined/netVisual_heatmap.pdf'), width =10,height=10)

#netVisual_diffInteraction
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'output_comparison_Refined/netVisual_diffInteraction.pdf'), width =10,height=10)

#Comparison of signalling Endo vs Control
gk1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gk2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gk1 + gk2
ggsave(plot = gk1 + gk2,
       filename = paste0(dir_out, 'output_comparison_Refined/rankNet.pdf'), width =10,height=12)

####################################
#### Plot candidate pathways
####################################

dir.create(paste0(dir_out,"output_comparison_Refined/rankNet/"))
pathways_list <- intersect(cellchat.Cntrl@netP$pathway, cellchat.Endo@netP$pathway) 
for (PATHWAY in pathways_list) {
  pathways.show <- PATHWAY
  par(mfrow = c(1, 2), xpd = TRUE)
  ht <- list()
  for (i in 1:length(object.list)) {
    ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds", title.name = paste(pathways.show, "signaling ", names(object.list)[i]))
  }
  pdf(paste0(dir_out, 'output_comparison_Refined/rankNet/rankNet', PATHWAY, '.pdf'), width = 15, height = 9)
  ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
  dev.off()
}


pathways_list <- intersect(cellchat.Cntrl@netP$pathway, cellchat.Endo@netP$pathway)
for (PATHWAY in pathways_list) {
  pdf(paste0(dir_out,"/netAnalysis_signalingRole_network/netAnalysis_signalingRole_network_",PATHWAY,"_CTL.pdf"), width = 12, height = 3)
  netAnalysis_signalingRole_network(cellchat.Cntrl, signaling = PATHWAY, width = 24, height = 2.5, font.size = 10)
  dev.off()
  pdf(paste0(dir_out,"/netAnalysis_signalingRole_network/netAnalysis_signalingRole_network_",PATHWAY,"_Endo.pdf"), width = 12, height = 3)
  netAnalysis_signalingRole_network(cellchat.Endo, signaling = PATHWAY, width = 24, height = 2.5, font.size = 10)
  dev.off()
}


####################################
#### Identify the upgulated and down-regulated signaling ligand-receptor pairs
####################################

### show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("CXCL"), remove.isolate = FALSE) 
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'output_comparison_Refined/netVisual_bubble_CXCL.pdf'), width =25,height=7)
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("CCL"), remove.isolate = FALSE) 
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'output_comparison_Refined/netVisual_bubble_CCL.pdf'), width =15,height=5)
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("IGFBP"), remove.isolate = FALSE,targets.use = c(5:11)) 
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'output_comparison_Refined/netVisual_bubble_IGFBP.pdf'), width =15,height=5)
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("VEGF"), remove.isolate = FALSE) # sources.use = 4, targets.use = c(5:11),
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'output_comparison_Refined/netVisual_bubble_VEGF.pdf'), width =15,height=5)
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("TNF"), remove.isolate = FALSE) 
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'output_comparison_Refined/netVisual_bubble_TNF.pdf'), width =25,height=5)
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("COMPLEMENT"), remove.isolate = FALSE) 
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("ICAM"), remove.isolate = FALSE) 
ggsave(plot = last_plot(),
       filename = paste0(dir_out, 'output_comparison_Refined/netVisual_bubble_ICAM.pdf'), width =25,height=5)
netVisual_bubble(cellchat, comparison = c(1, 2),signaling = c("IGF"), remove.isolate = FALSE) 


###Circle Plots
setwd(paste0(dir_out,"output_comparison_Refined"))
PATHWAYS <- c( "TNF", "ICAM", "VEGF", "WNT")
for (PATHWAY in PATHWAYS) {
  pathways.show <- PATHWAY
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) 
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
  }
}
dev.off()

# Save plot
pdf(file = paste0("EndoCtl_",PATHWAY,"_circle_plot_edge10.pdf"))
par(mfrow = c(1, length(object.list)), xpd = TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


