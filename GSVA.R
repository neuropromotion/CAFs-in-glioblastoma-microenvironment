library(GSVA)
library(readxl)
library(dplyr)

add_subtypes <- function(obj, path="path_to_dir/GBM_subtypes_DEG.xlsx"){
  cluster.cancer <- WhichCells(obj, idents = 'GBM')
  cluster.cancer <- subset(obj, cells = cluster.cancer)
  
  cluster.cancer <- RunPCA(cluster.cancer, features = VariableFeatures(object = cluster.cancer))
  cluster.cancer <- FindNeighbors(cluster.cancer, dims = 1:5)
  cluster.cancer <- FindClusters(cluster.cancer, resolution = 0.02)
  cluster.cancer <- RunUMAP(object = cluster.cancer, dims = 1:5)
  
  genes <- read_excel(path)
  
  cluster.cancer <- JoinLayers(object = cluster.cancer)
  expr_matrix <- as.matrix(GetAssayData(object = cluster.cancer, assay = "RNA", slot = "data"))
  
  gene_sets <- list(
    "AC-like" = genes$AC,
    "MES1-like" = genes$MES1,
    "MES2-like" = genes$MES2,
    "OPC-like" = genes$OPC,
    "NPC1-like" = genes$NPC1,
    "NPC2-like" = genes$NPC2,
    'Glioblastoma Stem Cells' = genes$GSC
  )
  expr_matrix_filtered <- expr_matrix[apply(expr_matrix, 1, var) > 0, ]
  gsvapar <- gsvaParam(expr_matrix_filtered, gene_sets)
  gsva_results <- gsva(gsvapar)
  
  cell_annotations <- apply(gsva_results, 2, function(x) names(which.max(x)))
  
  cluster.cancer$GSVA_annotation <- cell_annotations
  #return(cluster.cancer)
  
  # merge
  if (!all(names(cell_annotations) %in% colnames(obj))) {
    stop("Ошибка: не все клетки из cell_annotations найдены в merged_filtered!")
  }
  # 1. Преобразуем фактор аннотаций в data.frame с баркодами
  annot_df <- data.frame(
    barcode = names(Idents(obj)),  # получаем баркоды
    cluster = as.character(Idents(obj)),  # текущие кластеры
    stringsAsFactors = FALSE
  )
  # 2. Создаем data.frame с новыми аннотациями
  new_annot_df <- data.frame(
    barcode = names(cell_annotations),
    new_cluster = cell_annotations,
    stringsAsFactors = FALSE
  )
  
  annot_df <- annot_df %>% 
    left_join(new_annot_df, by = "barcode") %>% 
    mutate(
      final_cluster = ifelse(is.na(new_cluster), cluster, new_cluster)
    )
  # 4. Создаем фактор с обновленными уровнями
  all_levels <- unique(c(levels(Idents(obj)), unique(cell_annotations)))
  annot_df$final_cluster <- factor(annot_df$final_cluster, levels = all_levels)
  
  # 5. Проверяем соответствие порядка баркодов
  identical(annot_df$barcode, colnames(obj))  # должно быть TRUE
  
  # 6. Если порядок совпадает, добавляем аннотацию в объект
  obj$final_annot <- annot_df$final_cluster
  Idents(obj) <- "final_annot"
  
  return(obj)
}
merged_filtered <- add_subtypes(merged_filtered)

