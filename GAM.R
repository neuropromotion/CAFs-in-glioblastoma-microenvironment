library(tradeSeq)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(scales)
library(BiocParallel)
library(dplyr)
library(viridis)

fitGAM2 <- function(seu_obj, DC1, n_knots = 6){
  # use only variable genes
  hvg <- VariableFeatures(seu_obj)
  if (length(hvg) == 0) stop("VariableFeatures не найдены!")
  
  # get raw count data
  expr_matrix <- GetAssayData(seu_obj, assay = "RNA", layer = "counts")[hvg, ]

  # save original order of cells
  pseudotime <- matrix(DC1, ncol = 1)
  rownames(pseudotime) <- colnames(seu_obj)
  
  cellWeights <- matrix(1, nrow = ncol(seu_obj), ncol = 1)
  rownames(cellWeights) <- colnames(seu_obj)
  
  BPPARAM <- BiocParallel::MulticoreParam(workers = parallel::detectCores() - 2)
  
  message("Запуск fitGAM для ", length(hvg), " генов на ", BPPARAM$workers, " ядрах...")
  
  # launch gam
  sce <- fitGAM(counts = as.matrix(expr_matrix), 
                pseudotime = pseudotime, 
                cellWeights = cellWeights,
                nknots = n_knots,
                parallel = TRUE,
                BPPARAM = BPPARAM)
  
  # statistics
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
  
  
  # pseudotime
  if (is.null(names(pseudotime))) {
    names(pseudotime) <- colnames(subset_obj)
  }
  pt <- pseudotime[colnames(subset_obj)]
  
  # use significant genes
  sig_res <- assocRes %>%
    as.data.frame() %>%
    filter(pvalue < 0.05) %>%
    arrange(desc(waldStat)) %>% 
    head(n_genes)
  
  sig_genes <- intersect(rownames(sig_res), rownames(subset_obj))
  
  # sort cells
  cells_order <- order(pt)
  pt_sorted <- pt[cells_order]
  
  expr_raw <- as.matrix(GetAssayData(subset_obj, layer = "data")[sig_genes, cells_order])
  
  # smoothing
  message("Smoothing expression curves...")

  expr_smooth <- t(apply(expr_raw, 1, function(x) {
    fit <- smooth.spline(seq_along(x), x, spar = 0.7)
    return(fit$y)
  }))
  
  rownames(expr_smooth) <- sig_genes
  colnames(expr_smooth) <- colnames(expr_raw)
  
  mat <- t(scale(t(expr_smooth)))
  mat[is.na(mat)] <- 0
  mat[mat > 2] <- 2
  mat[mat < -2] <- -2
  
  final_cells <- colnames(mat)
  final_clusters <- Idents(subset_obj)[final_cells]
  
  # check
  message(paste("Размер матрицы:", ncol(mat)))
  message(paste("Длина аннотации:", length(final_clusters)))
  
  top_anno <- HeatmapAnnotation(
    Cluster = final_clusters,
    col = list(Cluster = colors),
    show_annotation_name = FALSE
  )
  
  col_fun = colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b"))
  
  # plot hm
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


plot_hm_gam <- function(subset_obj, DC1, assocRes, colors, 
                        cluster_columns = FALSE, 
                        n_genes=50,
                        title="Smooth expression variation",
                        f_size=10,
                        cluster_rows=TRUE
){
  
  all_expr <- GetAssayData(subset_obj, layer = "data")
  present_genes <- intersect(rownames(assocRes)[assocRes$pvalue < 0.05], rownames(all_expr))
  
  if (length(present_genes) == 0) {
    stop("Нет значимых генов (p < 0.05), присутствующих в матрице объекта.")
  }
  if (length(present_genes) < n_genes) {
    warning("Requested n_genes exceeds number of significant genes. Using all available ", length(present_genes))
    n_genes <- length(present_genes)
  }
  sig_genes <- present_genes[1:n_genes]
  
  # sort cells along DC1 (pseudotime)
  cells_order <- order(DC1)
  dc1_sorted <- DC1[cells_order]
  
  expr_sorted <- all_expr[sig_genes, cells_order]
  
  clusters <- Idents(subset_obj)
  clusters_sorted <- clusters[cells_order]  
  
  top_anno <- HeatmapAnnotation(
    Cluster = clusters_sorted,
    col = list(Cluster = colors),
    show_legend = TRUE
  )
  
  expr_smooth <- t(apply(expr_sorted, 1, function(x) {
    loess(x ~ dc1_sorted, span = 0.15)$fitted
  }))
  
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
                               span = 0.2) {  
  
  sig_df <- assocRes[!is.na(assocRes$pvalue) & assocRes$pvalue < 0.05, ]
  if(nrow(sig_df) == 0) stop("Нет значимых генов!")
  
  sig_genes <- rownames(sig_df[order(sig_df$waldStat, decreasing = TRUE), ])
  sig_genes <- head(sig_genes, n_genes)
  
  cells_order <- order(DC1)
  dc1_sorted <- DC1[cells_order]
  
  all_expr <- GetAssayData(subset_obj, layer = "data")
  expr_sorted <- all_expr[sig_genes, cells_order]
  
  message("Smoothing...")
  expr_smooth <- t(apply(expr_sorted, 1, function(x) {
    loess(x ~ dc1_sorted, span = span)$fitted
  }))
  
  clusters_sorted <- Idents(subset_obj)[cells_order]
  top_anno <- HeatmapAnnotation(
    Cluster = clusters_sorted,
    col = list(Cluster = colors),
    show_legend = TRUE
  )
  
  mat <- t(scale(t(expr_smooth)))
  mat[mat > 3] <- 3; mat[mat < -3] <- -3
  
  col_fun = colorRamp2(c(-3, 0, 3), c("black", "blue2", "#F06543"))
  
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
