FAM_full <- FindAllMarkers(merged_filtered,
                       logfc.threshold = 0.25,
                       group.by = "short_annotation",
                       min.pct = 0.1,
                       only.pos = TRUE,
                       test.use = 'wilcox',
                       slot = 'data')
FAM_cancer <- FindAllMarkers(cluster.cancer,
                           logfc.threshold = 0.25, 
                           min.pct = 0.1,
                           only.pos = TRUE,
                           test.use = 'wilcox',
                           slot = 'data')


write_xlsx(FAM, "DEG_initial.xlsx")
write_xlsx(FAM_full, "DEG_final.xlsx") 
 
cluster.names <- c(
  'CAFs',
  'Endothelial cells',
  'T lymphocytes',
  'Macrophages/Microglia',
  'Oligodendrocytes',
  'GBM')

cluster.names <- c("CAFs",                
                   "Endothelial cells",
                   "T lymphocytes",
                   "Macrophages/Microglia",
                   "Oligodendrocytes",
                   "MES1-like",         
                   "MES2-like",          
                   "NPC1-like",               
                   "NPC2-like",
                   "AC-like",
                   "OPC-like",
                   "Glioblastoma Stem Cells")


cluster.names <- c( 'Glioblastoma Stem Cells',
                     "AC-like",
                     "OPC-like",
                     "MES1-like",
                     "MES2-like",
                     "NPC1-like",
                     "NPC2-like")

# GET top DEG genes:
top.genes <- character()    
for (cl in cluster.names) {
  cluster_genes <- FAM$gene[FAM$cluster == cl]   # genes of current cluster 
  
  # add 5 unique genes
  while (added < 5 && idx <= length(cluster_genes)) {
    g <- cluster_genes[idx]
    if (!(g %in% top.genes)) {                   # get only unique gene
      top.genes <- c(top.genes, g)
      added <- added + 1
    }
    idx <- idx + 1                           
  }
}


#############################______DOTPLOT_____###################################### 

dp <- DotPlot(merged_filtered, features = top.genes, group.by='short_annotation',
              col.min = -2.5, col.max=2.5, scale.min = 0, scale.max = 100) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 1)) + # Поворот подписей на оси X
  ggtitle("") 

svg("vln.svg", width = 11, height = 5) # for DEG 
dp + scale_color_gradientn(
  colours = c("black", "blue1", "deeppink1"),  # Четыре цвета
  values = c(0, 0.1, .55,  1),  # Точки привязки (нормализованные от 0 до 1)
  limits = c(-2.5, 2.5)  # Соответствует col.min и col.max в DotPlot
)
dev.off()

#############################______HEATMAP_____###################################### 

heatmap <- DoHeatmap(
  object = merged_filtered, 
  features = top.genes,  
  group.colors = cluster_colors,
  label = F, 
  disp.min = -2.5, 
  disp.max = 2.5)

heatmap + scale_fill_gradientn(
  colours = c("black", "blue1", "deeppink1"),  
 values = c(0, 0.1, .75,  1),  
  limits = c(-2.5, 2.5)  
)

jpeg("hm.jpg", width = 5500, height = 3100, res=350)
heatmap + scale_fill_gradientn(
  colours = c("black", 'blue', "deeppink1"), 
  limits = c(-2.5, 2.5)  
)
dev.off()