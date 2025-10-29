library(ComplexHeatmap)
library(circlize)

########## COLORS
cluster_colors <- c(
  "CAFs" = 'red',                
  "Endothelial cells" = 'blue',
  "T lymphocytes" = 'forestgreen',
  "Macrophages/Microglia" = 'purple3',
  "Oligodendrocytes" = 'dodgerblue',
  "MES1-like" ='chocolate1',         
  "MES2-like" = 'chocolate4',          
  "NPC1-like" = 'cyan2',               
  "NPC2-like" = 'slateblue1',
  "AC-like" = 'black',
  "OPC-like" = 'orchid1',
  "Glioblastoma Stem Cells" = 'gold1' 
)
cluster_colors_6 <- c('red','blue','forestgreen','purple3', 'deepskyblue','gold1')

########## BLEND UMAP EXPRESSION  
svg("vln.svg", width = 15, height = 6)
FeaturePlot(merged_filtered, features = c("MDK", "PTPRZ1"), blend = TRUE, order = T, blend.threshold = 0.5,
            reduction = "umap", cols = c("slategray4", "springgreen2", "blue"))
dev.off()

########## UMAP EXPRESSION  

svg('vln.svg', width = 8, height = 7)
 FeaturePlot(merged_filtered, "MKI67", order=T, min.cutoff = 0, max.cutoff = 4, slot = "data")  +
   
   scale_color_gradientn(
     colours = c("gray", "blue1", "deeppink1"),   
     values = c(0, 0.2, .65,  1),   
     limits = c(0, 4)   
   )
dev.off()


########## VIOLIN EXPRESSION 
svg('vln.svg', width = 9, height = 7)
VlnPlot(merged_filtered, features = 'SPP1',
        cols  = cluster_colors, group.by = 'final_annot',
        log = T, pt.size = 0)
dev.off()
?VlnPlot

########## CO_EXPRESSION 
prepare.data <- function(obj, 
                         vars, 
                         breaks = c(-Inf, 1, 2, 3, Inf), 
                         labels = c("0-1", "1-2", "2-3", '3+')){
  plot_data <- FetchData(obj, vars = vars, layer = "data")
  for (i in 3:length(vars)){
    plot_data[, paste0(vars[i], '_group')] <- cut(
      plot_data[, vars[i]],
      breaks = breaks,
      labels = labels
    )
  }
  return(plot_data)
}
config <- c("umap_1", "umap_2", "PTN", 'MDK', 
                "PTPRZ1", 'NCL')
plot_data = prepare.data(merged_filtered, config)
 
## sort values in order to get highly expressed dots upper
plot_data_sorted <- plot_data[order(plot_data$EGFR), ]
pack <- c('HBEGF', 'EGFR')

svg("vln.svg", width = 5, height = 5)
ggplot(plot_data_sorted, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(size = HBEGF_group, color = EGFR), alpha = 0.8) + 
  scale_color_gradientn(
    colours = c("darkgray", "blue1", "deeppink1"),   
    values = scales::rescale(c(0, 0.3, 0.7, 1))
  ) +
  theme_minimal() +
  ggtitle(paste0('Co-expression: ', pack[1], ' (size) vs ', pack[2], ' (color)')) +
  scale_size_manual(
    name = paste0(pack[1], " Expression"),
    values = c("0-1" = 0, "1-2" = 0.5, "2-3" = 1, '3+' = 3),
    breaks = c("0-1", "1-2", "2-3", '3+')
  )
dev.off()


########## HEATMAP / DOTPLOT 


II.section <- c('TGFB1', 'TGFBR2', 'ACVRL1', 'ENG', 'NRP1', 'ID1', 'ID3','GRN', 'SORT1','JAG1', 'NOTCH1','NOTCH2','NOTCH3', 'NOTCH4',
                'IGFBP7', 'CD93','SPP1', 'ITGAV', 'ITGA5', 'ITGB1', 'ITGA6','SEMA3F', 'NPR2', 'PLXNA2', 'PTN', 'SDC3', 'NCL',
                'PRSS3', 'PPIA', 'BSG', 'PGF', 'FLT1', 'NAMPT', 'MDK', 'LGALS9', 'P4HB', 'FGF1', 'FGFR1',
                'CXCL12', 'CXCR4', 'ANGPT2', 'ANGPTL2', 'ADM', 'CALCRL','PGF', 'ANGPT2', 'BGN', 'HIF1A','VIM', 'FN1', 'CDH2',
                'HES1', 'HES2', 'HES3', 'HES4', 'HEY1', 'HEY2', 'HEYL', 'ZEB1', 'ZEB2')

II.section.dot <- c('ANGPT2', 'ITGA5', 'ITGB1', 
                    'GRN', 'SORT1', 'JAG1', 'NOTCH4', 'IGFBP7', 'CD93')

III.section <- c('TNFRSF1B', 'TNF', 'SPP1', 'CD44', 'PTPRC', 'PTN', 'NCL', 'PROS1', 'MERTK', 'AXL',
                 'MIF', 'CD74', 'MDK', 'LRP1', 'LGALS9', 'P4HB', 'HAVCR2', 'CXCL12', 'CXCR4', 'CCL5',
                 'CCL3L1', 'CCR1', 'CCL3', 'CCL2', 'C3', 'C3AR1', 'ITGAX', 'ITGB2', 'ANXA1', 'FPR1', 'PDCD1', 'CD274', 'PDCD1LG2')

III.section.dot <- c('SPP1', 'CD44', 'PTN', 'NCL', 'MDK', 'LRP1', 'ANXA1', 'FPR1')

validation_2 <- c('PTN', 'MDK', 'PTPRZ1', 'PPIA', 'BSG', 'LRP1', 'ITGA6', 'ITGB1', 'ITGA5', 'JAG1',
                        'NOTCH4', 'IGFBP7', 'CD93', 'FGF1', 'FGFR1', 'ANXA1', 'FPR1', 'ANGPT2', 'PGF',
                        'FLT1', 'TGFB1', 'PDGFB', 'PDGFRB', 'HBEGF','EGFR', 'ADM', 'CALCRL')

validation_1 <- c('SPP1', 'ITGA5', 'ITGA4', 'ITGA6', 'ITGB5', 'ITGB5', 'MDK', 'PTN',
                   'PTPRZ1', 'NCL', 'JAG1', 'NOTCH4', 'IGFBP7', 'CD93', 'ANXA1', 'FPR1',
                   'ANGPT2', 'PPIA', 'BSG', 'LRP1', 'GRN', 'SORT1', 'HBEGF', 'EGFR', 'PDGFA', 'PDGFRA')
caf.markers <- c('ACTA2', 'VIM', 'PDGFRB', 
                 'COL1A1', 'COL3A1', 'COL4A2', 
                 'NR2F2', 'BGN', 'FAP', 
                 'THY1', 'RGS5', 'FN1')
dp <- DotPlot(merged_filtered, features = II.section.dot, group.by='short_annotation',
              col.min = -2.5, col.max=2.5, scale.min = 0, scale.max = 100) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 1)) + # Поворот подписей на оси X
  ggtitle("") 

svg("vln.svg", width = 6.5, height = 5.5) # 
dp + scale_color_gradientn(
  colours = c("black", "blue1", "deeppink1"),  # Четыре цвета
  values = c(0, 0.1, .55,  1),  # Точки привязки (нормализованные от 0 до 1)
  limits = c(-2.5, 2.5)  # Соответствует col.min и col.max в DotPlot
)
dev.off()




#### complex HEATMAP
# two supplimentary functions 
get_data_for_hm <- function(seurat.obj, genes){
  # Извлеките матрицу из "data" слоя по группам
  expr_mat <- GetAssayData(seurat.obj, assay = "RNA", layer = "data")[genes, ]
  # Агрегируйте по аннотации (усреднение по кластерам)
  annotation <- seurat.obj$short_annotation
  expr_agg <- aggregate(t(as.matrix(expr_mat)), by = list(annotation), FUN = mean)
  rownames(expr_agg) <- expr_agg$Group.1
  expr_agg <- expr_agg[, -1]  # Удалите столбец групп
  return(expr_agg)
}
get_hm <- function(obj, genes, width=4, height=3){
  expr <- get_data_for_hm(obj, genes)
  
  
  expr_numeric <- as.matrix(expr)
  expr_numeric <- apply(expr_numeric, 2, as.numeric)
  
  
  col_fun <- colorRamp2(c(min(expr_numeric, na.rm = TRUE), 
                          median(expr_numeric, na.rm = TRUE), 
                          max(expr_numeric, na.rm = TRUE)), 
                        c("black", "blue1", "deeppink1"))
  
  
  hm <- Heatmap(
    t(expr),  # Транспонируйте для генов в строках
    name = "Expression",  # Название шкалы
    #col = colorRamp2(c(min(expr_agg), max(expr_agg)), c("black", 'cyan', "blue", "magenta")),  # Цвета
    col = col_fun,
    column_names_rot = 90,
    cluster_rows = F,
    cluster_columns = F,
    column_names_side = "top",  # Метки столбцов сверху!
    row_names_side = "left",    # Метки строк слева
    row_names_gp = gpar(fontsize = 7),         
    column_names_gp = gpar(fontsize = 7), 
    heatmap_legend_param = list(
      title = "Expr",                   # Заголовок шкалы
      title_position = "topcenter",     # По центру сверху
      legend_direction = "vertical",  # Горизонтальная шкала
      legend_width = unit(5, "cm")      # Ширина шкалы
      #title_gp = gpar(fontsize = 10),           # Шрифт заголовка шкалы
      #labels_gp = gpar(fontsize = 8)            # Шрифт меток на шкале
    ),  # Шкала по центру (см. ниже)
    border = 'black',
    rect_gp = gpar(col = "gray", lwd = .2),
    #width = unit(10, "cm"),             # Общая ширина
    #height = unit(10, "cm"),  
    
    # Размер ячеек: фиксированный
    width = ncol(t(expr)) * unit(width, "mm"),   # 6 (ANGIO)
    height = nrow(t(expr)) * unit(height, "mm")  # 2.5 (ANGIO)
  )
  return(hm)
}

hm <- get_hm(val, validation_1)
 

pdf("heatmap_60_genes.pdf", width = 4, height = 5)
draw(hm)
dev.off()
