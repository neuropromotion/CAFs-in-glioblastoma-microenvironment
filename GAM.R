library(tradeSeq)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(scales)


# fitGAM2 <- function(seu_obj, DC1){
#   
#   subset_obj <- FindVariableFeatures(seu_obj)
#   
#   hvg_expr <- GetAssayData(subset_obj, layer = "data")[VariableFeatures(subset_obj), ]
#   
#   cells_order <- order(DC1)
#   dc1_sorted <- DC1[cells_order]
#   expr_sorted <- hvg_expr[, cells_order]
#   
#   pseudotime <- matrix(DC1, ncol = 1)
#   cellWeights <- matrix(1, nrow = ncol(expr_sorted), ncol = 1)
#   
#   sce <- fitGAM(expr_sorted, pseudotime = pseudotime, cellWeights = cellWeights)
#   
#   assocRes <- associationTest(sce)
#   
#   gc()
#   return(list(gam=sce, assocres=assocRes))
# }
# 

fitGAM2 <- function(seu_obj, DC1, n_knots = 6){
  library(tradeSeq)
  library(BiocParallel)
  
  # 1. Используем уже существующие VariableFeatures
  hvg <- VariableFeatures(seu_obj)
  if (length(hvg) == 0) stop("VariableFeatures не найдены!")
  
  # 2. Извлекаем сырые каунты (tradeSeq любит именно их)
  # Если у тебя v5 и JoinLayers сделан, это будет одна матрица
  expr_matrix <- GetAssayData(seu_obj, assay = "RNA", layer = "counts")[hvg, ]
  
  # 3. Подготовка метаданных (сохраняем оригинальный порядок клеток)
  pseudotime <- matrix(DC1, ncol = 1)
  rownames(pseudotime) <- colnames(seu_obj)
  
  cellWeights <- matrix(1, nrow = ncol(seu_obj), ncol = 1)
  rownames(cellWeights) <- colnames(seu_obj)
  
  # 4. Настройка параллелизации (используем доступные ядра)
  # Для Linux серверов MulticoreParam подходит лучше всего
  BPPARAM <- BiocParallel::MulticoreParam(workers = parallel::detectCores() - 2)
  
  message("Запуск fitGAM для ", length(hvg), " генов на ", BPPARAM$workers, " ядрах...")
  
  # 5. Запуск модели
  sce <- fitGAM(counts = as.matrix(expr_matrix), 
                pseudotime = pseudotime, 
                cellWeights = cellWeights,
                nknots = n_knots,
                parallel = TRUE,
                BPPARAM = BPPARAM)
  
  # 6. Статистический тест
  message("Выполнение associationTest...")
  assocRes <- associationTest(sce)
  
  gc()
  return(list(gam = sce, assocres = assocRes))
}






plot_hm_gam_v2 <- function(subset_obj, pseudotime, assocRes, colors, 
                           n_genes = 50,
                           title = "Smooth expression variation",
                           f_size = 8,
                           cluster_rows = TRUE){
  
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  
  # 1. Подготовка псевдовремени
  if (is.null(names(pseudotime))) {
    names(pseudotime) <- colnames(subset_obj)
  }
  pt <- pseudotime[colnames(subset_obj)]
  
  # 2. Отбор генов
  sig_res <- assocRes %>%
    as.data.frame() %>%
    filter(pvalue < 0.05) %>%
    arrange(desc(waldStat)) %>% 
    head(n_genes)
  
  sig_genes <- intersect(rownames(sig_res), rownames(subset_obj))
  
  # 3. Сортировка клеток
  cells_order <- order(pt)
  pt_sorted <- pt[cells_order]
  
  # Извлекаем данные и СРАЗУ фиксируем матрицу
  expr_raw <- as.matrix(GetAssayData(subset_obj, layer = "data")[sig_genes, cells_order])
  
  # 4. Сглаживание
  message("Smoothing expression curves...")
  # Применяем сглаживание и ПРИНУДИТЕЛЬНО возвращаем матрицу того же размера
  expr_smooth <- t(apply(expr_raw, 1, function(x) {
    fit <- smooth.spline(seq_along(x), x, spar = 0.7)
    return(fit$y)
  }))
  
  # Восстанавливаем имена строк и колонок (ЭТО ВАЖНО)
  rownames(expr_smooth) <- sig_genes
  colnames(expr_smooth) <- colnames(expr_raw)
  
  # 5. Масштабирование
  mat <- t(scale(t(expr_smooth)))
  mat[is.na(mat)] <- 0
  mat[mat > 2] <- 2
  mat[mat < -2] <- -2
  
  # 6. Аннотация (берем строго по колонкам mat)
  final_cells <- colnames(mat)
  final_clusters <- Idents(subset_obj)[final_cells]
  
  # ПРОВЕРКА ДЛЯ ТЕБЯ (выведется в консоль)
  message(paste("Размер матрицы:", ncol(mat)))
  message(paste("Длина аннотации:", length(final_clusters)))
  
  top_anno <- HeatmapAnnotation(
    Cluster = final_clusters,
    col = list(Cluster = colors),
    show_annotation_name = FALSE
  )
  
  col_fun = colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b"))
  
  # 7. Хитмап БЕЗ column_split (чтобы точно не упало)
  p <- Heatmap(
    mat,
    name = "Z-score",
    column_title = title,
    col = col_fun,
    cluster_columns = FALSE, 
    cluster_rows = cluster_rows,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = f_size),
    top_annotation = top_anno,
    use_raster = TRUE,
    raster_quality = 3
  )
  
  return(p)
}










library(ComplexHeatmap)
library(circlize)
library(viridis)

plot_hm_gam <- function(subset_obj, DC1, assocRes, colors, 
                        cluster_columns = FALSE, 
                        n_genes=50,
                        title="Smooth expression variation",
                        f_size=10,
                        cluster_rows=TRUE
){
  
  
  all_expr <- GetAssayData(subset_obj, layer = "data")
  
  # 1. Находим пересечение: значимые гены, которые РЕАЛЬНО есть в матрице
  present_genes <- intersect(rownames(assocRes)[assocRes$pvalue < 0.05], rownames(all_expr))
  
  if (length(present_genes) == 0) {
    stop("Нет значимых генов (p < 0.05), присутствующих в матрице объекта.")
  }
  
  # 2. Отбираем топ N (используем именно present_genes)
  if (length(present_genes) < n_genes) {
    warning("Requested n_genes exceeds number of significant genes. Using all available ", length(present_genes))
    n_genes <- length(present_genes)
  }
  sig_genes <- present_genes[1:n_genes]
  
  # 3. Сортируем клетки по DC1
  cells_order <- order(DC1)
  dc1_sorted <- DC1[cells_order]
  
  # 4. Вырезаем матрицу (только нужные гены и клетки)
  expr_sorted <- all_expr[sig_genes, cells_order]
  
  # 5. Аннотация кластеров
  clusters <- Idents(subset_obj)
  clusters_sorted <- clusters[cells_order]  
  
  top_anno <- HeatmapAnnotation(
    Cluster = clusters_sorted,
    col = list(Cluster = colors),
    show_legend = TRUE
  )
  
  # 6. Сглаживание LOESS
  # Важно: работаем уже с expr_sorted, где точно нет лишних генов
  expr_smooth <- t(apply(expr_sorted, 1, function(x) {
    loess(x ~ dc1_sorted, span = 0.15)$fitted
  }))
  
  # 7. Исправляем фактор уровней (используем переданный объект)
  clusters_sorted <- factor(clusters_sorted, levels = levels(Idents(subset_obj)))
  
  gc()
  
  mat <- t(scale(t(expr_smooth)))
  z_values <- as.vector(mat)
  df_plot <- data.frame(z = z_values)
  
  z_score_plot <- ggplot(df_plot, aes(x = z)) +
    geom_histogram(aes(y = ..density..), bins = 100, fill = "steelblue", alpha = 0.5) +
    geom_density(color = "red", size = 1) +
    # Добавим линии стандартных границ Z-score (-2, 0, 2)
    geom_vline(xintercept = c(-3, 3), linetype = "dashed", color = "black") +
    # Добавим линии твоих реальных хвостов (99-й квантиль)
    geom_vline(xintercept = quantile(z_values, c(0.001, 0.999)), color = "darkred") +
    labs(title = "Z-score distribution in matrix for heatmap",
         subtitle = "Dashed black lines: -3, 3 | Red lines: 0.1% и 99.9% quantiles",
         x = "Z-score (scaled LOESS values)",
         y = "Density") +
    theme_minimal()
  
  col_fun = colorRamp2(
    c(-3, 0, 3), 
    c("black", "blue2", "#F06543")
  )
  # 8. Строим хитмап
  p <- Heatmap(
    mat,
    cluster_columns = cluster_columns,
    cluster_rows = cluster_rows,
    name = "Z-score",
    column_title = title, #c("black","blue2","#F06543")
    #col = colorRamp2(c(-2,0,1,2), c("black","blue2",'yellow',"#F06543")),
    #col = colorRamp2(c(-1, 0, 2), c("black", "blue2", "#F06543")),
    col = col_fun,
    use_raster = TRUE,
    raster_quality = 5,
    show_row_names = TRUE,
    top_annotation = top_anno,
    row_names_gp = gpar(fontsize = f_size),
    column_gap = unit(0, "mm"),
    border = FALSE,
    rect_gp = gpar(col = NA)
  )
  
  # Отрисовка для получения порядка строк
  p_drawn <- draw(p)
  row_indices <- row_order(p_drawn)
  
  if (is.list(row_indices)) {
    row_indices <- unlist(row_indices)
  }
  
  ordered_genes <- rownames(expr_smooth)[row_indices]
  
  return(list(plot = p, z_plot = z_score_plot, genes = ordered_genes))
}





 

plot_hm_gam_robust <- function(subset_obj, DC1, assocRes, colors, 
                               n_genes=50,
                               title="Adaptive Smooth Expression",
                               f_size=10,
                               span = 0.2) { # Добавил span для контроля гладкости
  library(ComplexHeatmap)
  library(circlize)
  
  # 1. Отбор генов по статистике GAM (WaldStat)
  # Это гарантирует, что мы берем биологически важные гены
  sig_df <- assocRes[!is.na(assocRes$pvalue) & assocRes$pvalue < 0.05, ]
  if(nrow(sig_df) == 0) stop("Нет значимых генов!")
  
  sig_genes <- rownames(sig_df[order(sig_df$waldStat, decreasing = TRUE), ])
  sig_genes <- head(sig_genes, n_genes)
  
  # 2. Сортируем клетки по их реальному положению в пространстве (DC1)
  cells_order <- order(DC1)
  dc1_sorted <- DC1[cells_order]
  
  # 3. Извлекаем данные (data layer - логнормализация)
  all_expr <- GetAssayData(subset_obj, layer = "data")
  expr_sorted <- all_expr[sig_genes, cells_order]
  
  # 4. Адаптивное сглаживание (LOESS)
  # Используем loess прямо по клеткам. 
  # Если переход в датасете резкий - уменьши span до 0.1
  # Если хочешь супер-гладко (как в прошлом датасете) - оставь 0.2
  message("Smoothing...")
  expr_smooth <- t(apply(expr_sorted, 1, function(x) {
    # Метод 'lowess' или 'loess' по реальным координатам DC1
    loess(x ~ dc1_sorted, span = span)$fitted
  }))
  
  # 5. Аннотация кластеров (реальная, по клеткам)
  clusters_sorted <- Idents(subset_obj)[cells_order]
  top_anno <- HeatmapAnnotation(
    Cluster = clusters_sorted,
    col = list(Cluster = colors),
    show_legend = TRUE
  )
  
  # 6. Z-score с защитой от выбросов
  mat <- t(scale(t(expr_smooth)))
  mat[mat > 3] <- 3; mat[mat < -3] <- -3
  
  # 7. Цвета (твоя классика)
  col_fun = colorRamp2(c(-3, 0, 3), c("black", "blue2", "#F06543"))
  
  # 8. Хитмап
  p <- Heatmap(
    mat,
    cluster_columns = FALSE, 
    cluster_rows = TRUE, 
    name = "Z-score",
    column_title = title,
    col = col_fun,
    top_annotation = top_anno,
    use_raster = TRUE,
    raster_quality = 5,
    show_row_names = TRUE, # Показываем гены, если их не слишком много
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = f_size),
    border = FALSE,
    rect_gp = gpar(col = NA)
  )
  
  return(list(plot = p, genes = rownames(mat)))
}

# plot_hm_gam <- function(seu_obj, DC1, assocRes, colors, 
#                         cluster_columns = FALSE, 
#                         n_genes=50,
#                         title="Smooth expression variation",
#                         f_size=10
#                         ){
#   
#   subset_obj <- FindVariableFeatures(seu_obj)
#   
#   hvg_expr <- GetAssayData(subset_obj, layer = "data")[VariableFeatures(subset_obj), ]
#   
#   cells_order <- order(DC1)
#   dc1_sorted <- DC1[cells_order]
#   expr_sorted <- hvg_expr[, cells_order]
#   
#   
#   clusters <- Idents(seu_obj)
#   
#   clusters_sorted <- clusters[cells_order]  
#   
#   top_anno <- HeatmapAnnotation(
#     Cluster = clusters_sorted,
#     col = list(Cluster = colors),
#     show_legend = TRUE
#   )
#   
#    
#   sig_genes <- rownames(assocRes)[assocRes$pvalue < 0.05]
#   if (length(sig_genes) < n_genes) {
#     warning("Requested n_genes exceeds number of significant genes. Using all available ", length(sig_genes), " genes.")
#     n_genes <- length(sig_genes)
#   }
#   sig_genes <- sig_genes[1:n_genes]
#   
#   
#   expr_smooth <- t(apply(expr_sorted[sig_genes, ], 1, function(x) {
#     loess(x ~ dc1_sorted, span = 0.15)$fitted
#   }))
#   
#   clusters_sorted <- factor(clusters_sorted, levels = levels(Idents(seu_obj)))
#   
#   gc()
#   
#   p <- Heatmap(
#     t(scale(t(expr_smooth))),
#     cluster_columns = cluster_columns,   # сохраняем порядок клеток по DC1
#     cluster_rows = TRUE,       # можно кластеризовать гены по форме
#     name = "Z-score",
#     column_title = title,
#     col = colorRamp2(c(-2,0,2), c("black","blue2","#F06543")),
#     use_raster = TRUE,
#     raster_quality = 5,
#     show_row_names = TRUE,
#     top_annotation = top_anno,
#     #column_split = clusters_sorted,
#     
#     # row_sizes
#     #height = unit(row_h * n_genes, "mm"), # Общая высота хитмапа зависит от кол-ва генов
#     row_names_gp = gpar(fontsize = f_size), # Размер шрифта генов
#     
#     
#     column_gap = unit(0, "mm"),
#     border = FALSE,
#     rect_gp = gpar(col = NA)
#   )
#   
#   # --- ДОБАВЬТЕ ЭТИ СТРОКИ ---
#   p_drawn <- draw(p) # заставляем хитмап рассчитать кластеризацию
#   row_indices <- row_order(p_drawn) # достаем индексы строк
#   
#   if (is.list(row_indices)) {
#     row_indices <- unlist(row_indices)
#   }
#   
#   # Имена генов в порядке отображения
#   ordered_genes <- rownames(expr_smooth)[row_indices]
#   
#   
#   # Возвращаем и график, и список генов
#   return(list(plot = p, genes = ordered_genes))
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
