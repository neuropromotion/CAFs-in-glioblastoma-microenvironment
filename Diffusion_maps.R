library(destiny)
library(Seurat)
library(dplyr)
library(plotly)
library(scatterplot3d)
library(viridis)
library(ComplexHeatmap)
library(scales)
library(circlize)
library(Matrix)

make_dm <- function(seurat_object, n_pcs = 50, seed = 42){
  set.seed(seed)
  
  hvg <- VariableFeatures(seurat_object)
  
  if (length(hvg) == 0) {
    stop("Variable genes are not found")
  }
  
  hvg_expr <- GetAssayData(seurat_object, layer = "data")[hvg, ]
  dm_hvg <- destiny::DiffusionMap(as.matrix(t(hvg_expr)), n_pcs = n_pcs)
  return(dm_hvg)
}


run_dpt <- function(dm_object) {
  # Root is the first (or last - it is not important actually) cell on DC1 component
  root_cell_index <- which.min(eigenvectors(dm_object)[, 1])
  message(paste("Root cell selected at index:", root_cell_index))
  
  dpt_res <- destiny::DPT(dm_object, tips = root_cell_index)
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
  
  cluster_ids <- as.character(Idents(obj))
  unique_clusters <- levels(Idents(obj))
  
  if (is.null(colors)) {
    library(scales)
    cols_vector <- hue_pal()(length(unique_clusters))
    names(cols_vector) <- unique_clusters
    fill_cols <- cols_vector[cluster_ids]
  } else {
    fill_cols <- colors[cluster_ids]
  }
  
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
                                   min_cutoff = "q5", 
                                   max_cutoff = "q95") {  
  
  
  data_to_plot <- FetchData(obj, vars = feature, layer = "data")
  expr_values <- data_to_plot[[1]]
  
  if (is.character(min_cutoff)) {
    q_min <- quantile(expr_values, as.numeric(gsub("q", "", min_cutoff))/100)
    expr_values[expr_values < q_min] <- q_min
  }
  if (is.character(max_cutoff)) {
    q_max <- quantile(expr_values, as.numeric(gsub("q", "", max_cutoff))/100)
    expr_values[expr_values > q_max] <- q_max
  }
  
  df <- data.frame(
    dc1 = DC1,
    dc2 = DC2,
    expression = expr_values
  )
  
  p <- ggplot(df, aes(x = dc1, y = dc2, fill = expression)) +
    geom_point(shape = 21, size = size, stroke = stroke, color = "black", alpha = alpha) +
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

get_3d_plot_dm <- function(dc, cluster_ids, colors = NULL, title = 'Diffusion Map DC1-DC2-DC3') {
  
  cluster_ids <- as.factor(cluster_ids)
  unique_clusters <- levels(cluster_ids)
  
  if (is.null(colors)) {
    cols_vector <- hue_pal()(length(unique_clusters))
    colors_to_use <- cols_vector
  } else {
    colors_to_use <- colors
  }
  
  p <- plot_ly(
    x = dc[[1]], 
    y = dc[[2]], 
    z = dc[[3]],
    color = cluster_ids, 
    colors = colors_to_use,  
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


get_smooth_drivers <- function(sce, assoc_res, top_n = 20, cutoff_cor = 0.8) {
  sig_genes <- rownames(assoc_res[assoc_res$pvalue < 0.05, ])
  yhat <- predictSmooth(sce, gene = sig_genes, nPoints = 100)
  
  time_points <- 1:nPoints
  monotonicity <- apply(yhat, 1, function(x) cor(x, time_points, method = "spearman"))
  
  df_drivers <- data.frame(
    gene = sig_genes,
    wald_stat = assoc_res[sig_genes, "waldStat"],
    cor_val = monotonicity,
    trend = ifelse(monotonicity > 0, "Up", "Down")
  )
  
  smooth_drivers <- df_drivers[abs(df_drivers$cor_val) >= cutoff_cor, ]
  smooth_drivers <- smooth_drivers[order(-smooth_drivers$wald_stat), ]
  
  return(head(smooth_drivers, top_n))
}

get_3d_plot_expr <- function(seu_obj, genes, dc, cluster_colors = NULL,
                             title = "3D DiffusionMap",
                             pt_size = 4,       
                             pt_opacity = 1) {  
  

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
  
  # layer with clusters
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
  
  # layer with gene expression
  for (gene in genes) {
    vec <- expr_data[[gene]]
    vmin <- quantile(vec, 0.05, na.rm = TRUE)
    vmax <- quantile(vec, 0.95, na.rm = TRUE)
    
    p <- p %>% add_trace(
      x = dc_df$DC1, y = dc_df$DC2, z = dc_df$DC3,
      type = "scatter3d", mode = "markers",
      marker = list(
        size = pt_size,  
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




get_hm_transition_prob <- function(seu_obj, dm, colors, group_by = NULL, eps = 1e-6) {

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
  
  M_ord <- M_mat[ord, ord] 
  M_log <- log10(M_ord + eps) 
   
  clusters <- Idents(seu_obj)[ord]
  clusters <- factor(clusters, levels = unique(clusters))
   
  clusters <- factor(clusters, levels = unique(clusters))
  if (!all(levels(clusters) %in% names(colors))) {
    stop("Not all clusters have assigned colors!")
  }
  
  # legend
  leg_attr <- list(
    title_gp = gpar(fontsize = 14, fontface = "bold"),  
    labels_gp = gpar(fontsize = 12),                  
    grid_height = unit(8, "mm"),                       
    grid_width = unit(8, "mm")                        
  )

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
    show_legend = FALSE 
  )
  
  ha_col <- HeatmapAnnotation(
    Cluster = clusters, 
    col = list(Cluster = colors), 
    annotation_height = unit(3, "mm"), 
    show_annotation_name = FALSE,
    annotation_legend_param = list(Cluster = leg_attr) 
  )
  

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


