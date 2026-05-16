library(Seurat)
library(Matrix)
library(data.table)
library(DoubletFinder)
library(ggplot2)
library(celda)
library(SingleCellExperiment)
#library(scran)
library(harmony)
library(tidyr)
library(ComplexHeatmap) 
library(circlize)


main_initial_colors <- c("CAFs" = '#00B295', 
                       'Pericytes'='#F56476',
                       'Endothelial cells'='#314CB6',
                       'ECs'='#314CB6',
                       'Oligodendrocytes'='dodgerblue',
                       'TAMs'='#642CA9',
                       'GBM'='gold3'
)

main_final_colors <- c("CAFs" = '#00B295', 
                      'Pericytes'='#F56476',
                      'Endothelial cells'='#314CB6',
                      'ECs'='#314CB6',
                      'Mesenchymal ECs'='gold',
                      'TAMs'='#642CA9',
                      'Oligodendrocytes'='dodgerblue',
                      'MES1-like'='chocolate1',
                      'MES2-like'='chocolate4',          
                      'NPC1-like'='#5AD2F4',               
                      'NPC2-like'='slateblue1',
                      'OPC-like'='orchid3',
                      'AC-like'='black',
                      'GSCs'='gold3'
)

main_gams_color <- c('Pericytes'='#F56476',
                    'ECs'='#314CB6',
                    'Mesenchymal ECs'='gold2',
                    'CAFs'='#00B295'
)


val_initial_colors <- c('Mesenchymal cells'='#F56476',
                      'Endothelial cells'='#314CB6',
                      'T lymphocytes'='forestgreen',
                      'TAMs'='#642CA9',
                      'Oligodendrocytes'='dodgerblue',
                      'GBM'='gold2'
)

val_final_colors <- c('Mesenchymal cells'='#F56476',
                      'CAFs' = '#00B295',
                      'Pericytes' = '#F56476',
                       'Endothelial cells'='#314CB6',
                       'T lymphocytes'='forestgreen',
                       'TAMs'='#642CA9',
                       'Oligodendrocytes'='dodgerblue',
                       'MES1-like'='chocolate1',
                       'MES2-like'='chocolate4',          
                       'NPC1-like'='#5AD2F4',               
                       'NPC2-like'='slateblue1',
                       
                       'OPC-like'='orchid3',
                       'AC-like'='black',
                       'GSCs'='gold3'
)
val_gams_color <- c('Pericytes'='#F56476',
                 'ECs'='#314CB6',
                 'CAFs'='#00B295'
)


save_html <- function(fig, path){
  htmlwidgets::saveWidget(fig, path, selfcontained = TRUE)
}



change_order <- function(seu, new.order, group.by) {
  
  # 1. Проверка наличия колонки
  if(!group.by %in% colnames(seu@meta.data)) {
    stop(paste("Column", group.by, "not found in metadata"))
  }
  
  # 2. Получаем текущие уникальные значения
  current_values <- unique(as.character(seu@meta.data[[group.by]]))
  
  # 3. Проверка: совпадают ли переданные имена с теми, что есть в данных
  if(!setequal(current_values, new.order)) {
    missing <- setdiff(current_values, new.order)
    extra <- setdiff(new.order, current_values)
    stop(paste("Mismatch in cluster names.\nMissing in new.order:", 
               paste(missing, collapse=", "), 
               "\nExtra in new.order:", paste(extra, collapse=", ")))
  }
  
  # 4. Перезаписываем колонку как фактор с НОВЫМ порядком уровней
  # Это гарантирует, что везде (DimPlot, DotPlot, VlnPlot) порядок будет как в new.order
  seu@meta.data[[group.by]] <- factor(seu@meta.data[[group.by]], levels = new.order)
  
  # 5. Устанавливаем эту колонку как активную идентичность (опционально, но удобно)
  Idents(seu) <- group.by
  
  return(seu)
}


add_GBM_subtypes <- function(seu, annot_name, path='/home/amismailov/genes.xlsx'){
  
  seu[[annot_name]] <- as.character(Idents(seu))
  gbm <- subset(seu, idents = c('GBM'))
  if(ncol(gbm) == 0) stop("No GBM cells found!")
  
  modules <- readxl::read_excel(path)
  
  gbm <- AddModuleScore(
    gbm,
    features = list(
      AC   = modules$AC,
      NPC1 = modules$NPC1,
      NPC2 = modules$NPC2,
      MES1 = modules$MES1,
      MES2 = modules$MES2,
      OPC = modules$OPC,
      GSC = modules$GSC
    ),
    name = c("AC", "NPC1", "NPC2", 'MES1', 'MES2', 'OPC', 'GSC')
  )
  
  scores <- gbm@meta.data[, c("AC1", "NPC12", "NPC23", 'MES14', 'MES25', 'OPC6', 'GSC7')]
  
  gbm$GBM_state <- colnames(scores)[
    apply(scores, 1, which.max)
  ]
  
  
  Idents(gbm) <- 'GBM_state'
  new_names <- c('AC1'='AC-like','GSC7'='GSCs','MES14'='MES1-like','MES25'='MES2-like','NPC12'='NPC1-like','NPC23'='NPC2-like','OPC6'='OPC-like')
  gbm <- RenameIdents(gbm, new_names)
  gbm$final <- as.character(Idents(gbm))
  
  
  common_cells <- intersect(colnames(gbm), colnames(seu))
  
  # Создаем новый вектор аннотаций
  current_annots <- as.character(Idents(seu))
  names(current_annots) <- colnames(seu)
  
  # Заменяем GBM на их специфичные стейты
  current_annots[common_cells] <- gbm$final
  
  # Сохраняем в объект
  seu[[annot_name]] <- current_annots
  
  # Опционально: ставим новую аннотацию как активные Idents
  # Idents(seu) <- annot_name
  
  gc()
  return(seu)
  
}


make_2d_plot <- function(seu, 
                         colors = NULL, 
                         comp_1 = 1,
                         comp_2 = 2,
                         group.by = NULL,
                         shape = 21,
                         size = 2,
                         stroke = 0.2,
                         alpha = 0.8) {
  
  library(ggplot2)
  library(Seurat)
  
  # 1. Логика выбора группировки (group.by или Idents)
  if (is.null(group.by)) {
    cell_groups <- Idents(seu)
  } else {
    cell_groups <- seu[[group.by]][, 1]
  }
  
  # 2. Формируем данные для визуализации
  df <- data.frame(
    dim1 = Embeddings(seu, "umap")[, comp_1],
    dim2 = Embeddings(seu, "umap")[, comp_2],
    cell_type = cell_groups
  )
  
  # 3. Базовый объект графика
  p <- ggplot(df, aes(x = dim1, y = dim2, fill = cell_type)) +
    geom_point(shape = shape, size = size, stroke = stroke, color = "black", alpha = alpha) +
    theme_classic() +
    labs(x = paste0("UMAP_", comp_1), y = paste0("UMAP_", comp_2)) +
    theme(
      legend.key.size = unit(1.75, "lines"),
      legend.text = element_text(size = 13),
      legend.title = element_blank()
    ) +
    guides(
      fill = guide_legend(override.aes = list(size = 6))
    )
  
  # 4. Логика раскраски
  if (!is.null(colors)) {
    # Если палитра передана — используем её
    p <- p + scale_fill_manual(values = colors)
  } else {
    # Если цветов нет — используем стандартную палитру ggplot2 (discrete)
    # Это предотвратит ошибку отсутствия аргумента 'values'
    p <- p + scale_fill_discrete()
  }
  
  return(p)
}





make_hm_complex <- function(seu, genes,
                            clusters = NULL, # Новый аргумент для списка кластеров
                            cluster_rows = FALSE,
                            cluster_cols = FALSE,
                            side = "right",
                            layer = 'data',
                            font_size_row = 10,
                            font_size_col = 12,
                            z_limit = 2,
                            lwd = 1,
                            column_title = '',
                            assay = 'RNA') {
  
  # 1. Получаем средние значения для существующих данных
  # Важно: AverageExpression берет текущие Idents(seu)
  avg <- AverageExpression(seu, features = genes, layer = layer)
  mat_raw <- as.matrix(avg[[assay]])
  
  # 2. Определяем итоговый список кластеров
  # Если список не передан, берем те, что есть в объекте
  if (is.null(clusters)) {
    clusters <- colnames(mat_raw)
  }
  
  # 3. Создаем "полную" матрицу (Гены x Кластеры)
  # Заполняем NA по умолчанию
  mat_full <- matrix(NA, 
                     nrow = length(genes), 
                     ncol = length(clusters), 
                     dimnames = list(genes, clusters))
  
  # Находим пересечение по генам и по кластерам
  existing_genes <- intersect(genes, rownames(mat_raw))
  existing_clusters <- intersect(clusters, colnames(mat_raw))
  
  # Заполняем только те ячейки, где данные есть и в списке, и в объекте
  if(length(existing_genes) > 0 && length(existing_clusters) > 0) {
    mat_full[existing_genes, existing_clusters] <- mat_raw[existing_genes, existing_clusters]
  }
  
  # 4. Z-score по строкам
  # scale() корректно обрабатывает NA: он их игнорирует при расчете среднего/SD
  # и возвращает NA на тех же позициях.
  mat_scaled <- t(apply(mat_full, 1, function(x) {
    # Если в строке совсем нет данных или нет вариации — оставляем NA
    if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) {
      return(rep(NA, length(x)))
    }
    # Масштабируем (автоматически сохраняет NA на пустых местах)
    return(as.vector(scale(x)))
  }))
  
  colnames(mat_scaled) <- clusters
  rownames(mat_scaled) <- genes
  
  # Clipping
  mat_scaled[mat_scaled > z_limit] <- z_limit
  mat_scaled[mat_scaled < -z_limit] <- -z_limit
  
  # 5. Цветовая шкала
  col_fun <- colorRamp2(
    seq(-z_limit, z_limit, length.out = 3),
    c('gray', "navy", "firebrick1")#'coral', 
  )
  
  # 6. Heatmap
  ht <- Heatmap(
    mat_scaled,
    name = "Z-score",
    col = col_fun,
    #na_col = "#000", # Черный для отсутствующих генов И кластеров
    
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_title = column_title,
    row_names_side = side,
    rect_gp = gpar(col = "black", lwd = lwd),
    row_names_gp = gpar(fontsize = font_size_row),
    column_names_gp = gpar(fontsize = font_size_col),
    border = TRUE
  )
  
  ht_drawn <- draw(ht)
  
  # 7. Возврат порядка генов
  if (cluster_rows) {
    row_idx <- row_order(ht_drawn)
    if (is.list(row_idx)) row_idx <- unlist(row_idx)
    ordered_genes <- rownames(mat_scaled)[row_idx]
  } else {
    ordered_genes <- genes
  }
  
  
  return(list(
    plot = ht_drawn,
    genes = ordered_genes,
    matrix = mat_scaled  # <-- Добавляем это
  ))
}









make_hm <- function(seu, genes,
                    cluster_rows = FALSE,
                    cluster_cols = FALSE,
                    side = "right",
                    layer = 'data',
                    z_limit = 3) { # Добавляем лимит для Z-score
  
  # 1. Получаем средние значения
  avg <- AverageExpression(seu, features = genes, layer = layer)
  mat <- as.matrix(avg$RNA)
  
  existing_genes <- genes[genes %in% rownames(mat)]
  mat <- mat[existing_genes, , drop = FALSE]
  
  # 2. Ручной Z-score с клиппингом (вместо scale = "row" в pheatmap)
  # Это дает контроль над визуализацией
  mat_scaled <- t(scale(t(mat))) 
  mat_scaled[mat_scaled > z_limit] <- z_limit
  mat_scaled[mat_scaled < -z_limit] <- -z_limit
  
  # 3. Визуализация
  p <- pheatmap(
    mat_scaled, # Используем подготовленную матрицу
    scale = "none", # scale уже сделан вручную
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    #color = colorRampPalette(c("blue", "white", "red"))(100), # Стандарт для Z-score
    border_color = "gray",
    fontsize_row = 8,
    name = "Z-score",
    row_names_side = side
  )
  
  # ... остальная часть кода с draw и возвратом ...

  
  # 3. Для ComplexHeatmap нужно вызвать draw, чтобы объект "ожил"
  p_drawn <- draw(p)
  
  # 4. Извлекаем порядок
  if (cluster_rows) {
    # Если кластеризация включена, берем результат вычислений
    row_idx <- row_order(p_drawn)
    # Если row_order возвращает список (бывает при split), берем первый элемент
    if(is.list(row_idx)) row_idx <- unlist(row_idx)
    ordered_genes <- rownames(mat)[row_idx]
  } else {
    # Если выключена — это будет наш existing_genes
    ordered_genes <- rownames(mat)
  }
  
  return(list(
    plot = p,
    genes = ordered_genes
  ))
}



library(viridis)
make_feature_plot <- function(
    seu,
    feature,
    comp_1 = 1,
    comp_2 = 2,
    reduction = "umap",
    shape = 21,
    size = 2,
    stroke = 0.2,
    alpha = 0.9, 
    order = FALSE,
    title=NULL
) {
  
  library(Seurat)
  library(ggplot2)
  
  df <- data.frame(
    dim1 = Embeddings(seu, reduction)[, comp_1],
    dim2 = Embeddings(seu, reduction)[, comp_2],
    expression = FetchData(seu, vars = feature)[,1]
  )
  
  if (order) {
    df <- df[order(df$expression), ]
  }
  
  # если title не задан — используем имя гена
  if (is.null(title)) {
    title <- feature
  }
  
  p <- ggplot(df, aes(x = dim1, y = dim2, fill = expression)) +
    geom_point(
      shape = shape,
      size = size,
      stroke = stroke,
      color = "black",
      alpha = alpha
    ) +
    scale_fill_viridis_c(option = "inferno") +
    labs(title = title) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.key.size = unit(1.75, "lines"),
      legend.text = element_text(size = 13),
      legend.title = element_blank()
    )
  
  return(p)
}

get_2d_umap_expr <- function(
    seu_obj,
    genes,
    reduction = "umap",
    cluster_colors,
    title = "2D UMAP (Clusters / Gene expression)",
    point_size=6
) {
  
  library(plotly)
  library(Seurat)
  library(viridis)
  
  # --- embeddings ---
  emb <- Embeddings(seu_obj, reduction)
  
  df <- data.frame(
    UMAP1 = emb[,1],
    UMAP2 = emb[,2],
    cluster = as.character(Idents(seu_obj))
  )
  
  #clusters <- unique(as.character(df$cluster))
  clusters <- levels(Idents(seu_obj))
  
  # --- expression matrix ---
  expr_mat <- GetAssayData(seu_obj, layer = "data")[genes, ]
  expr_mat <- log1p(expr_mat)
  
  vmin <- quantile(expr_mat, 0.01, na.rm = TRUE)
  vmax <- quantile(expr_mat, 0.99, na.rm = TRUE)
  
  inferno_scale <- lapply(
    seq(0, 1, length.out = 256),
    function(x) list(x, viridis::inferno(256)[round(x * 255) + 1])
  )
  
  # --- base plot ---
  p <- plot_ly()
  
  # =========================
  # CLUSTER TRACES
  # =========================
  
  for (k in seq_along(clusters)) {
    
    clust_name <- clusters[k]
    
    if (!clust_name %in% names(cluster_colors)) {
      stop(paste("No color provided for cluster:", clust_name))
    }
    
    idx <- df$cluster == clust_name
    
    p <- p %>%
      add_trace(
        x = df$UMAP1[idx],
        y = df$UMAP2[idx],
        type = "scatter",
        mode = "markers",
        marker = list(
          size = point_size,
          color = cluster_colors[[clust_name]]
        ),
        name = clust_name,
        showlegend = TRUE,
        visible = TRUE
      )
  }
  
  # =========================
  # GENE TRACES (hidden)
  # =========================
  
  for (i in seq_along(genes)) {
    
    expr_vec <- expr_mat[i, ]
    expr_vec[expr_vec < vmin] <- vmin
    expr_vec[expr_vec > vmax] <- vmax
    
    p <- p %>%
      add_trace(
        x = df$UMAP1,
        y = df$UMAP2,
        type = "scatter",
        mode = "markers",
        marker = list(
          size = point_size,
          color = expr_vec,
          colorscale = inferno_scale,
          cmin = vmin,
          cmax = vmax,
          showscale = TRUE
        ),
        name = genes[i],
        showlegend = FALSE,
        visible = FALSE
      )
  }
  
  # =========================
  # DROPDOWN
  # =========================
  
  buttons <- list()
  
  # clusters button
  buttons[[1]] <- list(
    label = "Clusters",
    method = "update",
    args = list(
      list(visible = c(rep(TRUE, length(clusters)),
                       rep(FALSE, length(genes))))
    )
  )
  
  # gene buttons
  for (i in seq_along(genes)) {
    
    vis <- c(rep(FALSE, length(clusters)),
             rep(FALSE, length(genes)))
    
    vis[length(clusters) + i] <- TRUE
    
    buttons[[i + 1]] <- list(
      label = genes[i],
      method = "update",
      args = list(list(visible = vis))
    )
  }
  
  # =========================
  # LAYOUT
  # =========================
  
  p <- p %>%
    layout(
      title = title,
      xaxis = list(
        title = "UMAP_1",
        zeroline = FALSE,
        showline = FALSE,
        showgrid = FALSE
      ),
      yaxis = list(
        title = "UMAP_2",
        zeroline = FALSE,
        showline = FALSE,
        showgrid = FALSE
      ),
      updatemenus = list(
        list(
          type = "dropdown",
          x = 0.02,
          y = 1,
          buttons = buttons
        )
      )
    )
  
  return(p)
}

 
decontex_workflow <- function(seu){
  seurat_copy <- seu
  
  counts_mat <- GetAssayData(seurat_copy, layer = "counts", assay = "RNA") # raw counts matrix extraction
  cell_md <- seurat_copy@meta.data # metadata extraction

  sce <- SingleCellExperiment(
    assays = list(counts = counts_mat),
    colData = cell_md
  )
  # Add clusters
  get_groups <- function(sobj){
    sobj <- NormalizeData(sobj, verbose = FALSE)
    sobj <- FindVariableFeatures(object = sobj, nfeatures = 2000, verbose = FALSE, selection.method = 'vst')
    sobj <- ScaleData(sobj, verbose = FALSE)
    sobj <- RunPCA(sobj, npcs = 20, verbose = FALSE)
    sobj <- FindNeighbors(sobj, dims = 1:20, verbose = FALSE)
    sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE)
    return(sobj@meta.data[['seurat_clusters']])
  }
  add_groups <- function(sobj){
    sobj$soup_group <- get_groups(sobj)
    return(sobj)
  }
  seurat_copy <- add_groups(seurat_copy)
  sce <- decontX(sce, z = seurat_copy$soup_group) # run decontX
  # Extraction clear count and metadata
  clean_counts <- decontXcounts(sce)
  clean_md     <- as.data.frame(colData(sce))
  # create new clean Seurat-object
  seu_clean <- CreateSeuratObject(
    counts    = clean_counts,
    project   = "GBM_decontaminated",
    meta.data = clean_md
  )
  seu_clean$contamination <- sce$decontX_contamination
  return(seu_clean)
} # DecontX pipeline



standard_workflow <- function(obj, group.by.vars=NULL, n_components=10, res=0.1){
  seu <- NormalizeData(obj)                             
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)   
  seu <- ScaleData(seu, features = rownames(seu)) 
  seu <- RunPCA(seu, features = VariableFeatures(object = seu))
  if(!is.null(group.by.vars)){
    seu <- RunHarmony(
      object = seu,
      group.by.vars = group.by.vars
    )
    seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:n_components)
    seu <- FindClusters(seu, resolution = res)
    seu <- RunUMAP(seu, reduction = "harmony", dims = 1:n_components)
  }
  else{
    seu <- FindNeighbors(seu, reduction = "pca", dims = 1:n_components)
    seu <- FindClusters(seu, resolution = res)
    seu <- RunUMAP(seu, reduction = "pca", dims = 1:n_components)
  }
  gc()
  return(seu)
} # Norm, Scale, PCA, UMAP

standard_workflow_sub <- function(obj, norm=F, scale=T, harmony.dims=20, theta=2,
                                  nfeatures=2000, group.by.vars=NULL, n_components=5, res=0.1){
  seu <- obj
  if(norm){
    seu <- NormalizeData(seu)    
  }
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = nfeatures) 
  if(scale){
    seu <- ScaleData(seu, features = rownames(seu)) 
  }
  
  seu <- RunPCA(seu, features = VariableFeatures(object = seu))
  
  if(!is.null(group.by.vars)){
    seu <- RunHarmony(
      object = seu,
      group.by.vars = group.by.vars,
      dims.use=1:harmony.dims,
      theta=theta
    )
    seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:n_components)
    seu <- FindClusters(seu, resolution = res)
    seu <- RunUMAP(seu, reduction = "harmony", dims = 1:n_components)
  }
  else{
    seu <- FindNeighbors(seu, reduction = "pca", dims = 1:n_components)
    seu <- FindClusters(seu, resolution = res)
    seu <- RunUMAP(seu, reduction = "pca", dims = 1:n_components)
  }
  gc()
  return(seu)
}


deduplex_v3 <- function(obj, join=FALSE) {
  # 1. Подготовка временного объекта для расчетов
  # Мы делаем копию, чтобы не засорять основной объект слоями PCA и ScaleData
  temp <- obj 
  
  if(join){
    temp <- JoinLayers(object = temp)
  }
  
  temp <- NormalizeData(temp)
  # 3. Стандартный цикл для DoubletFinder
  temp <- FindVariableFeatures(temp, verbose = FALSE)
  temp <- ScaleData(temp, verbose = FALSE)
  temp <- RunPCA(temp, npcs = 30, verbose = FALSE)
  
  # 4. Поиск дублетов
  # nExp: 0.05 - это ожидаемый процент дублетов (5%). 
  # Для 10X это значение обычно зависит от количества загруженных клеток.
  nExp <- round(0.05 * ncol(temp)) 
  
  temp <- doubletFinder(temp, PCs = 1:10, pN = 0.25, pK = 0.09, 
                        nExp = nExp, reuse.pANN = NULL, sct = FALSE)
  
  # 5. Извлекаем колонку с результатами классификации
  df_col <- grep("DF.classifications", colnames(temp@meta.data), value = TRUE)
  
  # 6. Переносим результаты в исходный объект
  # Нам нужны только метаданные, чтобы основной объект остался "чистым"
  res_table <- temp@meta.data[, df_col, drop = FALSE]
  obj <- AddMetaData(obj, metadata = res_table)
  
  cat('Cells Before: ', ncol(obj), ' | ')
  
  # 7. Фильтрация
  # Используем классификацию из только что добавленной колонки
  obj <- obj[, obj@meta.data[[df_col]] == "Singlet"]
  
  cat('Cells After: ', ncol(obj), '\n')
  
  return(obj)
}


deduplex_v4 <- function(obj, join = FALSE, min_cells = 100) {
  # Проверка на размер объекта
  if (ncol(obj) < min_cells) {
    warning("Объект содержит меньше ", min_cells, " клеток. Пропускаем DoubletFinder.")
    return(obj)
  }
  
  # Создаём временный объект
  temp <- obj
  
  if (join) {
    temp <- JoinLayers(object = temp)
  }
  
  # Нормализация и стандартная подготовка для DoubletFinder
  temp <- NormalizeData(temp)
  temp <- FindVariableFeatures(temp, verbose = FALSE)
  temp <- ScaleData(temp, verbose = FALSE)
  temp <- RunPCA(temp, npcs = 30, verbose = FALSE)
  
  nExp <- round(0.05 * ncol(temp))  # ожидаемый процент дублетов
  
  # Попытка запуска DoubletFinder
  temp <- tryCatch({
    doubletFinder(temp, PCs = 1:10, pN = 0.25, pK = 0.09,
                  nExp = nExp, reuse.pANN = NULL, sct = FALSE)
  }, error = function(e) {
    warning("DoubletFinder не сработал: ", e$message)
    return(NULL)
  })
  
  # Если DoubletFinder не сработал — возвращаем исходный объект
  if (is.null(temp)) return(obj)
  
  # Поиск колонок DF.classifications
  df_cols <- grep("DF.classifications", colnames(temp@meta.data), value = TRUE)
  
  if (length(df_cols) == 0) {
    warning("DF.classifications колонка не найдена. Пропускаем фильтр.")
    return(obj)
  }
  
  # Если колонок несколько — берём первую
  df_col <- df_cols[1]
  
  # Переносим результаты в исходный объект
  res_table <- temp@meta.data[, df_col, drop = FALSE]
  obj <- AddMetaData(obj, metadata = res_table)
  
  cat("Cells Before: ", ncol(obj), " | ")
  
  # Фильтрация Singlet
  singlet_idx <- which(obj@meta.data[[df_col]] == "Singlet")
  doublets <- which(obj@meta.data[[df_col]] != "Singlet")
  cat('Doublet size: ', length(doublets))
  if (length(singlet_idx) == 0) {
    warning("Нет клеток с классификацией Singlet. Возвращаем исходный объект.")
    return(obj)
  }
  
  obj <- obj[, singlet_idx]
  
  cat("Cells After: ", ncol(obj), "\n")
  
  return(obj)
}








deduplex_v3_integrated <- function(obj, stlit.by='orig.ident',
                                   project_name='seurat_DF_filtered'){
  # Divide and process 
  seu_list <- SplitObject(obj, split.by = stlit.by)
  seu_list <- lapply(seu_list, function(x) {
    if (ncol(x) < 100) return(x) 
    return(deduplex_v3(x))
  })
  
  # Merge
  seu <- merge(
    seu_list[[1]], 
    y = seu_list[-1], 
    project = project_name
  )
  
  seu <- JoinLayers(object = seu)
  gc()
  return(seu)
}


# LOAD SMARTSEQ2
smartseq_txt_integration <- function(path='/mnt/jack-5/amismailov/CAF_study/data/GSE135045/', 
                                     project_name='seurat_project'){
  path = path
  files <- list.files(path)
  seurat_list <- lapply(files, function(f){
    mat <- read.delim(paste0(path, f), row.names = 1)
    CreateSeuratObject(counts = mat)
  })
  names(seurat_list) <- files
  
  obj <- merge(seurat_list[[1]], y = seurat_list[-1],
                     add.cell.ids  = files,
                     project = project_name)
  
  obj$sample <- rownames(obj@meta.data)
  
  obj@meta.data <- separate(obj@meta.data, col = 'sample', into = c('Sample', 'Barcode'),
                                  sep = '_', extra = "merge")
  
  obj@meta.data$orig.ident <- paste0(project_name, '_', obj@meta.data$Sample)
  obj@meta.data$Sample <- NULL
  
  rm(seurat_list, files)
  gc()
  
  return(obj)
}






get_data_for_hm <- function(seurat.obj, genes){
  # Извлеките матрицу из "data" слоя по группам
  expr_mat <- GetAssayData(seurat.obj, assay = "RNA", layer = "data")[genes, ]
  # Агрегируйте по аннотации (усреднение по кластерам)
  annotation <- seurat.obj$short_annotation
  expr_agg <- aggregate(t(as.matrix(expr_mat)), by = list(annotation), FUN = mean)
  rownames(expr_agg) <- expr_agg$Group.1
  expr_agg <- expr_agg[, -1]  # Удалите столбец групп
  gc()
  return(expr_agg)
}
get_complex_hm <- function(obj, genes, width=4, height=3){
  expr <- get_data_for_hm(obj, genes)
  
  
  expr_numeric <- as.matrix(expr)
  expr_numeric <- apply(expr_numeric, 2, as.numeric)
  
  
  col_fun <- colorRamp2(c(min(expr_numeric, na.rm = TRUE), 
                          median(expr_numeric, na.rm = TRUE), 
                          max(expr_numeric, na.rm = TRUE)), 
                        c("black", "blue1", "#F06543"))
  
  
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
  gc()
  return(hm)
}


# GET INTERSECTION FAM

get_robust_markers <- function(
    seu,
    cluster_name,
    logfc.threshold = 0.25,
    min.pct = 0.1,
    padj.threshold = 0.05,
    test.use = "wilcox",
    only.pos = TRUE
) {
  
  # Проверка
  if (!cluster_name %in% levels(Idents(seu))) {
    stop("Cluster not found in Idents(seu)")
  }
  
  clusters <- levels(Idents(seu))
  other_clusters <- setdiff(clusters, cluster_name)
  
  markers_list_genes <- list()
  markers_list_full  <- list()
  
  for (cl in other_clusters) {
    
    markers <- FindMarkers(
      seu,
      ident.1 = cluster_name,
      ident.2 = cl,
      logfc.threshold = logfc.threshold,
      min.pct = min.pct,
      test.use = test.use,
      only.pos = only.pos
    )
    
    # сохраняем полную таблицу
    markers_list_full[[cl]] <- markers
    
    # фильтрация
    markers_filtered <- markers[markers$p_val_adj < padj.threshold, ]
    
    markers_list_genes[[cl]] <- rownames(markers_filtered)
  }
  
  # пересечение
  common_genes <- Reduce(intersect, markers_list_genes)
  
  return(common_genes)
}



# Differentially expressed genes
FAM2 <- function(seu,
                 n_genes = 50,
                 return_FAM = FALSE,
                 logfc.threshold = 0.25, 
                 min.pct = 0.1,
                 only.pos = TRUE,
                 test.use = 'wilcox',
                 layer = 'data'){
  
  FAM <- FindAllMarkers(seu,
                        logfc.threshold = logfc.threshold, 
                        min.pct = min.pct,
                        only.pos = only.pos,
                        test.use = test.use,
                        layer = layer)
  
  FAM <- FAM[FAM$p_val_adj<0.05,]
  
  clusters <- names(table(Idents(seu)))
  
  top_genes <- list()
  for (cl in clusters) {
    cluster_genes <- FAM$gene[FAM$cluster == cl][1:n_genes]  
    top_genes[[as.character(cl)]] <- cluster_genes
  }
  
  gc()
  if(return_FAM){
    return(c(FAM, top_genes))
    }
  
  return(top_genes)
  
}







library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(dplyr)

filter_and_remove_doublets <- function(seurat_obj) {
  # 1. Сплит объекта по образцам
  sample_list <- SplitObject(seurat_obj, split.by = "orig.ident")
  
  clean_samples <- lapply(names(sample_list), function(s_name) {
    message(paste0("--- Processing sample: ", s_name, " ---"))
    
    # Конвертация в SingleCellExperiment (требуется для scDblFinder)
    sce <- as.SingleCellExperiment(sample_list[[s_name]])
    
    # Запуск поиска дуплетов
    # clusters = TRUE помогает алгоритму лучше различать типы клеток
    sce <- scDblFinder(sce, clusters = TRUE)
    
    # Извлекаем результаты
    db_results <- as.data.frame(colData(sce))
    num_doublets <- sum(db_results$scDblFinder.class == "doublet")
    total_cells <- nrow(db_results)
    perc <- round((num_doublets / total_cells) * 100, 2)
    
    message(paste0("Found ", num_doublets, " doublets out of ", total_cells, " cells (", perc, "%)"))
    
    # Добавляем информацию о дуплетах обратно в метаданные Seurat подсета
    sample_list[[s_name]]$scDblFinder_class <- db_results$scDblFinder.class
    sample_list[[s_name]]$scDblFinder_score <- db_results$scDblFinder.score
    
    # Фильтруем: оставляем только синглеты
    subset_clean <- subset(sample_list[[s_name]], subset = scDblFinder_class == "singlet")
    
    return(subset_clean)
  })
  
  # 2. Собираем всё обратно в один объект
  message("--- Merging clean samples back together ---")
  merged_clean <- merge(clean_samples[[1]], y = clean_samples[-1])
  
  return(merged_clean)
}
 





 

PrettyDotPlot <- function(obj, genes, title = "Z-score Gene Expression", dot_scale = 8) {
  
  # 1. Проверяем наличие генов в объекте
  all_genes <- rownames(obj)
  genes_to_use <- intersect(genes, all_genes)
  
  if(length(genes_to_use) == 0) {
    stop("Ни один из указанных генов не найден в объекте!")
  }
  
  # 2. Принудительное масштабирование (чтобы был Z-score)
  # Используем только нужные гены, чтобы не перегружать память
  obj <- ScaleData(obj, features = genes_to_use, verbose = FALSE)
  
  # 3. Генерируем базовые данные DotPlot
  # col.min/max ограничивают Z-score для лучшего контраста
  p <- DotPlot(obj, features = genes_to_use, col.min = -2.5, col.max = 2.5)
  
  # 4. Полная кастомизация
  p_custom <- p + 
    # Слой с обводкой и заливкой Viridis
    geom_point(aes(size = pct.exp, fill = avg.exp.scaled), 
               shape = 21, 
               color = "black", 
               stroke = 0.6) + 
    
    # Цветовая шкала Magma (Z-score)
    scale_fill_viridis_c(option = "plasma", 
                         name = "Z-score\nExpression",
                         limits = c(-2.5, 2.5), 
                         oob = scales::squish) + 
    
    # Настройка размеров точек (убираем конфликт слоев)
    scale_size(range = c(0, dot_scale), name = "Percent\nExpressed") +
    
    # Оформление
    theme_minimal() +
    theme(
      # Исправлено: убраны лишние аргументы
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                                 size = 10, face = "bold", color = "black"),
      axis.text.y = element_text(size = 10, face = "bold", color = "black"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid.major = element_line(color = "gray95"),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    ) +
    
    # Удаляем стандартную легенду Seurat, оставляем нашу кастомную
    guides(color = "none") +
    
    labs(title = title, x = "", y = "")
  
  return(suppressMessages(suppressWarnings(p_custom)))
}

save_pdf <- function(plot, path, w=5, h=5){
  ggsave(
    filename = path,
    plot = plot,
    width = w,
    height = h,
    device = cairo_pdf   
  )
}

