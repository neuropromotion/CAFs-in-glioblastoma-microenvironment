library(destiny)
library(Seurat)
library(dplyr)
library(plotly)
library(scatterplot3d)
library(viridis)
library(ComplexHeatmap)

# make_dm <- function(seurat_object, seed=42){
#   set.seed(seed)
#   subset_obj <- FindVariableFeatures(seurat_object)
#   hvg_expr <- GetAssayData(subset_obj, layer = "data")[VariableFeatures(subset_obj), ]
#   dm_hvg <- DiffusionMap(as.matrix(t(hvg_expr)), n_pcs = 50)
#   
#   return(dm_hvg)
# }

make_dm <- function(seurat_object, n_pcs = 50, seed = 42){
  set.seed(seed)
  
  # 1. Проверяем, есть ли уже вариабельные гены в объекте
  hvg <- VariableFeatures(seurat_object)
  
  if (length(hvg) == 0) {
    stop("В объекте не найдены VariableFeatures. Сначала запустите FindVariableFeatures или передайте их в объект.")
  }
  
  # 2. Извлекаем данные только для этих генов
  # Используем слой "data" (логарифмированные данные)
  hvg_expr <- GetAssayData(seurat_object, layer = "data")[hvg, ]
  
  # 3. Строим DiffusionMap
  # Транспонируем, так как DiffusionMap ждет клетки в строках, а гены в столбцах
  dm_hvg <- destiny::DiffusionMap(as.matrix(t(hvg_expr)), n_pcs = n_pcs)
  
  return(dm_hvg)
}


run_dpt <- function(dm_object) {
  # 1. Находим индекс клетки с минимальным значением первой диффузионной компоненты (DC1)
  # Это и будет наш "root"
  root_cell_index <- which.min(eigenvectors(dm_object)[, 1])
  
  message(paste("Root cell selected at index:", root_cell_index))
  
  # 2. Расчет DPT с указанием этой точки как tip
  dpt_res <- destiny::DPT(dm_object, tips = root_cell_index)
  
  # 3. Извлечение значений псевдовремени
  # По умолчанию dpt_res содержит колонку dpt, рассчитанную от выбранного tip
  pseudotime_values <- dpt_res$dpt
  
  return(list(dpt_obj = dpt_res, time = pseudotime_values, root = root_cell_index))
}



get_2d_plot_dm <- function(obj, DC1, DC2, colors = NULL,
                           title = 'DiffusionMap (HVG)',
                           xlab = "DC1",
                           ylab = "DC2",
                           size = 1,
                           stroke = 0.5,
                           alpha = 1) {
  
  # Получаем уровни кластеров
  cluster_ids <- as.character(Idents(obj))
  unique_clusters <- levels(Idents(obj))
  
  # Если цвета не заданы, генерируем стандартную палитру (как в Seurat)
  if (is.null(colors)) {
    # Используем палитру из пакета scales, которая имитирует дефолтный ggplot2/Seurat
    library(scales)
    cols_vector <- hue_pal()(length(unique_clusters))
    names(cols_vector) <- unique_clusters
    fill_cols <- cols_vector[cluster_ids]
  } else {
    fill_cols <- colors[cluster_ids]
  }
  
  # Отрисовка
  plot(DC1, DC2,
       pch = 21,                                     
       bg = adjustcolor(fill_cols, alpha.f = alpha),  
       col = "black",                                
       lwd = stroke,                                 
       cex = size,
       main = title,
       xlab = xlab,
       ylab = ylab)
}

get_2d_feature_plot_dm <- function(obj, DC1, DC2, 
                                   feature, 
                                   title = NULL,
                                   size = 1, 
                                   stroke = 0.2, 
                                   alpha = 1,
                                   min_cutoff = "q5", # Обрезка шума снизу
                                   max_cutoff = "q95") { # Обрезка выбросов сверху
  
  
  # 1. Извлекаем данные экспрессии
  # FetchData удобен тем, что сам ищет и в генах, и в метаданных
  data_to_plot <- FetchData(obj, vars = feature, layer = "data")
  expr_values <- data_to_plot[[1]]
  
  # 2. Обработка выбросов (как в FeaturePlot)
  if (is.character(min_cutoff)) {
    q_min <- quantile(expr_values, as.numeric(gsub("q", "", min_cutoff))/100)
    expr_values[expr_values < q_min] <- q_min
  }
  if (is.character(max_cutoff)) {
    q_max <- quantile(expr_values, as.numeric(gsub("q", "", max_cutoff))/100)
    expr_values[expr_values > q_max] <- q_max
  }
  
  # 3. Собираем таблицу для ggplot
  df <- data.frame(
    dc1 = DC1,
    dc2 = DC2,
    expression = expr_values
  )
  
  # 4. Отрисовка
  p <- ggplot(df, aes(x = dc1, y = dc2, fill = expression)) +
    # Используем shape 21 для точек с обводкой
    geom_point(shape = 21, size = size, stroke = stroke, color = "black", alpha = alpha) +
    # Настраиваем градиент (классический Seurat-style или Viridis)
    scale_fill_gradientn(colors = inferno(100), name = "Expr") +
    theme_classic() +
    labs(
      title = if (is.null(title)) feature else title,
      x = "DC1",
      y = "DC2"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "italic"),
      legend.position = "right"
    )
  
  return(p)
}

# get_3d_plot_dm <- function(dc, cluster_ids, colors, title='Diffusion Map DC1-DC2-DC3'){
# 
#   p <- plot_ly(x = dc[[1]], y = dc[[2]], z = dc[[3]],
#           color = cluster_ids, colors = colors,
#           type = 'scatter3d', mode = 'markers', marker = list(size = 3)) %>%
#     layout(title = title)
# 
#   return(p)
# }

get_3d_plot_dm <- function(dc, cluster_ids, colors = NULL, title = 'Diffusion Map DC1-DC2-DC3') {
  library(plotly)
  library(scales)
  
  # Превращаем в фактор, чтобы сохранить порядок уровней
  cluster_ids <- as.factor(cluster_ids)
  unique_clusters <- levels(cluster_ids)
  
  # Если цвета не заданы, генерируем стандартную палитру Seurat (hue_pal)
  if (is.null(colors)) {
    cols_vector <- hue_pal()(length(unique_clusters))
    # Plotly удобнее принимать именованный вектор или просто список цветов
    colors_to_use <- cols_vector
  } else {
    colors_to_use <- colors
  }
  
  p <- plot_ly(
    x = dc[[1]], 
    y = dc[[2]], 
    z = dc[[3]],
    color = cluster_ids, 
    colors = colors_to_use, # Применит цвета по порядку уровней фактора
    type = 'scatter3d', 
    mode = 'markers', 
    marker = list(
      size = 5
    )
  ) %>%
    layout(
      title = title,
      scene = list(
        xaxis = list(title = "DC1"),
        yaxis = list(title = "DC2"),
        zaxis = list(title = "DC3")
      )
    )
  
  return(p)
}
# get_3d_plot_expr <- function(seu_obj, genes, dc, cluster_colors,
#                              title="3D DiffusionMap (Clusters / Gene expression)"){
#   # seu_obj - seurat object 
#   # genes - list of genes to show expression
#   # dc - list of diffusion map's tree components: list(dm$DC1, dm$DC2, dm$DC3)
#   # cluster_colors: vector of colors for each cluster: c('red', 'blue')
#   
#   dc_df <- data.frame(
#     DC1 = dc[[1]],
#     DC2 = dc[[2]],
#     DC3 = dc[[3]],
#     cluster = Idents(seu_obj)
#   )
#   
#   # -------- Clusters --------
#   clusters <- unique(as.character(dc_df$cluster))
#    
#   expr_mat <- GetAssayData(seu_obj, layer = "data")[genes, ]
#   expr_mat <- log1p(expr_mat)
#   
#   vmin <- quantile(expr_mat, 0.01, na.rm = TRUE)
#   vmax <- quantile(expr_mat, 0.99, na.rm = TRUE)
#   
#   inferno_scale <- lapply(
#     seq(0, 1, length.out = 256),
#     function(x) list(x, viridis::inferno(256)[round(x * 255) + 1])
#   )
#   
#   # -------- Base plot (clusters broken into traces) --------
#   p <- plot_ly()
#   
#   for (k in seq_along(clusters)) {
#     idx <- dc_df$cluster == clusters[k]
#     
#     p <- p %>%
#       add_trace(
#         x = dc_df$DC1[idx],
#         y = dc_df$DC2[idx],
#         z = dc_df$DC3[idx],
#         type = "scatter3d",
#         mode = "markers",
#         marker = list(
#           size = 4,
#           color = cluster_colors[clusters[k]]
#         ),
#         name = clusters[k],
#         showlegend = TRUE,
#         visible = TRUE
#       )
#   }
#   
#   # -------- Gene traces (hidden) --------
#   for (i in seq_along(genes)) {
#     p <- p %>%
#       add_trace(
#         x = dc_df$DC1,
#         y = dc_df$DC2,
#         z = dc_df$DC3,
#         type = "scatter3d",
#         mode = "markers",
#         marker = list(
#           size = 3,
#           color = expr_mat[i, ],
#           colorscale = inferno_scale,
#           cmin = vmin,
#           cmax = vmax,
#           showscale = TRUE
#         ),
#         name = genes[i],
#         showlegend = FALSE,
#         visible = FALSE
#       )
#   }
#   
#   # -------- Dropdown buttons --------
#   buttons <- list()
#   
#   buttons[[1]] <- list(
#     label = "Clusters",
#     method = "update",
#     args = list(
#       list(visible = c(rep(TRUE, length(clusters)), rep(FALSE, length(genes))))
#     )
#   )
#   
#   for (i in seq_along(genes)) {
#     vis <- c(rep(FALSE, length(clusters)), rep(FALSE, length(genes)))
#     vis[length(clusters) + i] <- TRUE
#     
#     buttons[[i + 1]] <- list(
#       label = genes[i],
#       method = "update",
#       args = list(list(visible = vis))
#     )
#   }
#   
#   # -------- Layout --------
#   p <- p %>%
#     layout(
#       title = title,
#       updatemenus = list(
#         list(
#           type = "dropdown",
#           x = 0.02,
#           y = 1,
#           buttons = buttons
#         )
#       )
#     )
#   
#   return(p)
# }


get_smooth_drivers <- function(sce, assoc_res, top_n = 20, cutoff_cor = 0.8) {
  # 1. Берем значимые гены (p-value < 0.05)
  sig_genes <- rownames(assoc_res[assoc_res$pvalue < 0.05, ])
  
  # 2. Получаем предсказанные значения (smoothers) для всех генов
  # Это "чистые" линии без шума клеток
  yhat <- predictSmooth(sce, gene = sig_genes, nPoints = 100)
  
  # 3. Считаем монотонность для каждого гена
  # (Корреляция между временем и предсказанием модели)
  time_points <- 1:100
  monotonicity <- apply(yhat, 1, function(x) cor(x, time_points, method = "spearman"))
  
  # 4. Собираем таблицу
  df_drivers <- data.frame(
    gene = sig_genes,
    wald_stat = assoc_res[sig_genes, "waldStat"],
    cor_val = monotonicity,
    trend = ifelse(monotonicity > 0, "Up", "Down")
  )
  
  # 5. Фильтруем только те, что выше порога плавности
  smooth_drivers <- df_drivers[abs(df_drivers$cor_val) >= cutoff_cor, ]
  smooth_drivers <- smooth_drivers[order(-smooth_drivers$wald_stat), ]
  
  return(head(smooth_drivers, top_n))
}

get_3d_plot_expr <- function(seu_obj, genes, dc, cluster_colors = NULL,
                             title = "3D DiffusionMap",
                             pt_size = 4,      # <-- Новый аргумент для размера
                             pt_opacity = 1) { # <-- Новый аргумент для прозрачности
  library(plotly)
  library(viridis)
  library(scales)
  
  # 1. Подготовка координат
  dc_df <- data.frame(
    DC1 = dc[[1]], DC2 = dc[[2]], DC3 = dc[[3]],
    cluster = Idents(seu_obj)
  )
  
  clusters <- levels(dc_df$cluster)
  if (is.null(cluster_colors)) {
    cluster_colors <- setNames(hue_pal()(length(clusters)), clusters)
  }
  
  expr_data <- FetchData(seu_obj, vars = genes, layer = "scale.data")
  
  inferno_scale <- lapply(seq(0, 1, length.out = 256), function(x) {
    list(x, viridis::inferno(256)[round(x * 255) + 1])
  })
  
  p <- plot_ly()
  
  # -------- Слой Кластеров --------
  for (k in seq_along(clusters)) {
    idx <- which(dc_df$cluster == clusters[k])
    p <- p %>% add_trace(
      x = dc_df$DC1[idx], y = dc_df$DC2[idx], z = dc_df$DC3[idx],
      type = "scatter3d", mode = "markers",
      marker = list(
        size = pt_size, 
        color = cluster_colors[clusters[k]],
        opacity = pt_opacity
      ),
      name = clusters[k], showlegend = TRUE, visible = TRUE
    )
  }
  
  # -------- Слои Генов --------
  for (gene in genes) {
    vec <- expr_data[[gene]]
    vmin <- quantile(vec, 0.05, na.rm = TRUE)
    vmax <- quantile(vec, 0.95, na.rm = TRUE)
    
    p <- p %>% add_trace(
      x = dc_df$DC1, y = dc_df$DC2, z = dc_df$DC3,
      type = "scatter3d", mode = "markers",
      marker = list(
        size = pt_size, # Используем тот же размер
        color = vec, 
        colorscale = inferno_scale, 
        cmin = vmin, cmax = vmax,
        showscale = TRUE,
        opacity = pt_opacity
      ),
      name = gene, showlegend = FALSE, visible = FALSE
    )
  }
  
  # -------- Кнопки и Layout --------
  # (Логика кнопок остается прежней)
  n_clusters <- length(clusters)
  n_genes <- length(genes)
  buttons <- list(list(label = "Clusters", method = "update", 
                       args = list(list(visible = c(rep(TRUE, n_clusters), rep(FALSE, n_genes))))))
  
  for (i in seq_along(genes)) {
    vis <- rep(FALSE, n_clusters + n_genes)
    vis[n_clusters + i] <- TRUE
    buttons[[i + 1]] <- list(label = genes[i], method = "update", args = list(list(visible = vis)))
  }
  
  p <- p %>% layout(
    title = title,
    updatemenus = list(list(type = "dropdown", x = 0.05, y = 1, buttons = buttons)),
    scene = list(xaxis = list(title = "DC1"), yaxis = list(title = "DC2"), zaxis = list(title = "DC3"))
  )
  
  return(p)
}



# get_3d_plot_expr <- function(seu_obj, genes, dc, cluster_colors = NULL,
#                              title = "3D DiffusionMap (Clusters / Gene expression)") {
#   library(plotly)
#   library(viridis)
#   library(scales)
#   
#   # 1. Подготовка данных координат и кластеров
#   dc_df <- data.frame(
#     DC1 = dc[[1]],
#     DC2 = dc[[2]],
#     DC3 = dc[[3]],
#     cluster = Idents(seu_obj)
#   )
#   
#   # 2. Логика цветов для кластеров (стандартная палитра, если NULL)
#   clusters <- levels(dc_df$cluster) # сохраняем порядок уровней
#   if (is.null(cluster_colors)) {
#     cluster_colors <- hue_pal()(length(clusters))
#     names(cluster_colors) <- clusters
#   }
#   
#   # 3. Извлечение данных экспрессии
#   # Используем FetchData — это надежнее, так как он сам найдет гены и в SCT, и в RNA
#   expr_data <- FetchData(seu_obj, vars = genes, layer = "data")
#   
#   # Если это один ген, превращаем в матрицу для совместимости с циклом
#   if (length(genes) == 1) {
#     expr_mat <- as.matrix(expr_data)
#     colnames(expr_mat) <- genes
#   } else {
#     expr_mat <- t(as.matrix(expr_data))
#   }
#   
#   # Расчет лимитов для шкалы (чтобы выбросы не портили градиент)
#   vmin <- quantile(expr_mat, 0.01, na.rm = TRUE)
#   vmax <- quantile(expr_mat, 0.99, na.rm = TRUE)
#   
#   # Создаем шкалу Inferno
#   inferno_scale <- lapply(
#     seq(0, 1, length.out = 256),
#     function(x) list(x, viridis::inferno(256)[round(x * 255) + 1])
#   )
#   
#   p <- plot_ly()
#   
#   # -------- Слой Кластеров --------
#   for (k in seq_along(clusters)) {
#     idx <- dc_df$cluster == clusters[k]
#     p <- p %>%
#       add_trace(
#         x = dc_df$DC1[idx], y = dc_df$DC2[idx], z = dc_df$DC3[idx],
#         type = "scatter3d", mode = "markers",
#         marker = list(size = 4, color = cluster_colors[clusters[k]]),
#         name = clusters[k], showlegend = TRUE, visible = TRUE
#       )
#   }
#   
#   # -------- Слои Генов (скрыты по умолчанию) --------
#   for (i in seq_along(genes)) {
#     p <- p %>%
#       add_trace(
#         x = dc_df$DC1, y = dc_df$DC2, z = dc_df$DC3,
#         type = "scatter3d", mode = "markers",
#         marker = list(
#           size = 3, color = expr_mat[genes[i], ],
#           colorscale = inferno_scale, cmin = vmin, cmax = vmax,
#           showscale = TRUE
#         ),
#         name = genes[i], showlegend = FALSE, visible = FALSE
#       )
#   }
#   
#   # -------- Кнопки переключения --------
#   buttons <- list()
#   
#   # Кнопка для возврата к кластерам
#   buttons[[1]] <- list(
#     label = "Clusters",
#     method = "update",
#     args = list(list(visible = c(rep(TRUE, length(clusters)), rep(FALSE, length(genes)))))
#   )
#   
#   # Кнопки для каждого гена
#   for (i in seq_along(genes)) {
#     vis <- c(rep(FALSE, length(clusters)), rep(FALSE, length(genes)))
#     vis[length(clusters) + i] <- TRUE
#     
#     buttons[[i + 1]] <- list(
#       label = genes[i],
#       method = "update",
#       args = list(list(visible = vis))
#     )
#   }
#   
#   p <- p %>%
#     layout(
#       title = title,
#       updatemenus = list(
#         list(type = "dropdown", x = 0.02, y = 1, buttons = buttons)
#       ),
#       scene = list(xaxis = list(title = "DC1"), yaxis = list(title = "DC2"), zaxis = list(title = "DC3"))
#     )
#   
#   return(p)
# }
library(ComplexHeatmap)
library(viridis)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(Matrix)

get_hm_transition_prob <- function(seu_obj, dm, colors, group_by = NULL, eps = 1e-6) {
  # dm - diffusion map object, where DM coordinates dm$DC1, dm$DC2 ... are presented 
  # colors - map of colors, like: 
  dc1_coords <- data.frame(
    cluster = Idents(seu_obj),
    DC1 = dm$DC1,
    DC2 = dm$DC2
  )
  
  
  M <- dm@transitions    
  
  M_mat <- as.matrix(M)
  
  dc1 <- dm$DC1
  
  tapply(dc1, Idents(seu_obj), median)
  
  dc1 <- -dc1
  ord <- order(dc1)
  
  M_ord <- M_mat[ord, ord] # sort cells along DC1
  

  M_log <- log10(M_ord + eps) # add eps to stop crushing 
   
  clusters <- Idents(seu_obj)[ord]
  clusters <- factor(clusters, levels = unique(clusters))
   
  clusters <- factor(clusters, levels = unique(clusters))
  if (!all(levels(clusters) %in% names(colors))) {
    stop("Not all clusters have assigned colors!")
  }
  
  # --- НАСТРОЙКИ  ЛЕГЕНДЫ ---
  leg_attr <- list(
    title_gp = gpar(fontsize = 14, fontface = "bold"),  
    labels_gp = gpar(fontsize = 12),                  
    grid_height = unit(8, "mm"),                       
    grid_width = unit(8, "mm")                        
  )
  # --- НАСТРОЙКИ  ЛЕГЕНДЫ кластеров---
  bar_attr <- list(
    title_gp = gpar(fontsize = 14, fontface = "bold"),
    labels_gp = gpar(fontsize = 12),
    legend_height = unit(60, "mm"),                   
    grid_width = unit(8, "mm")                       
  )
  
  ha_row <- rowAnnotation(
    Cluster = clusters, 
    col = list(Cluster = colors), 
    annotation_height = unit(3, "mm"), 
    show_annotation_name = FALSE,
    show_legend = FALSE # Чтобы не дублировать легенду слева и сверху
  )
  
  ha_col <- HeatmapAnnotation(
    Cluster = clusters, 
    col = list(Cluster = colors), 
    annotation_height = unit(3, "mm"), 
    show_annotation_name = FALSE,
    annotation_legend_param = list(Cluster = leg_attr) # Применяем размер
  )
  
  ########  ########  ########  ########  ########
  ########  ########  ########  ########  ########
  #col_range <- quantile(M_log, probs = c(0.05, 0.95), na.rm = TRUE)
  breaks <- seq(quantile(M_log, 0.01), quantile(M_log, 0.99), length.out = 100)
  colors_func <- colorRamp2(breaks, inferno(100))
  ht <- Heatmap(
    M_log,
    name = "log10(P)",
    col = colors_func,
    #colorRamp2(seq(col_range[1], col_range[2], length.out = 100), inferno(100)),
    show_row_names = FALSE,
    show_column_names = FALSE,
    left_annotation = ha_row, 
    top_annotation = ha_col,
    heatmap_legend_param = bar_attr,
    cluster_rows = TRUE,          
    cluster_columns = TRUE,        
    row_gap = unit(0, "mm"),
    column_gap = unit(0, "mm"),
    show_parent_dend_line = FALSE, 
    rect_gp = gpar(col = "gray", lwd = 0.1),
    border = TRUE,
    use_raster = TRUE,
    raster_quality = 4 
  )
  
  
  return(ht)
}










