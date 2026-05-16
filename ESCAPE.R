library(enrichR)
library(org.Hs.eg.db)
library(GO.db)
library(AnnotationDbi)
library(stringr)
library(magrittr)
library(escape)
library(ComplexHeatmap)
library(viridis)

run_enrich <- function(FAM, n_genes=250){
  up_sets <- list()
  clusters <- unique(FAM$cluster)
  
  
  
  for (clust in clusters) {
    deg_clust <- FAM[FAM$cluster == clust, ]
    deg_clust <- deg_clust[deg_clust$p_val_adj < 0.05 & deg_clust$avg_log2FC > 0.5, ]
    
    # Топ-100 (или все, если меньше 100)
    if (nrow(deg_clust) > 0) {
      top_idx <- order(deg_clust$avg_log2FC, decreasing = TRUE)[1:min(n_genes, nrow(deg_clust))]
      up_genes <- rownames(deg_clust)[top_idx]
      up_sets[[clust]] <- up_genes
    }
  }
  
  
  enrich_results = list()
  for (clust in clusters) {
    current_results <- enrichr(up_sets[[clust]], databases = c(
      "GO_Biological_Process_2023",
      "GO_Molecular_Function_2023"#, 
      #"KEGG_2021_Human",
      #"Hallmark_Gene_Sets",
      #"Jensen_DISEASES"
    ))
    enrich_results[[clust]] = current_results
  }
  gc()
  return(enrich_results)
}


run_enrich_list <- function(genes){

  current_results <- enrichr(genes, databases = c(
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023"#, 
    #"KEGG_2021_Human"#,
    #"Hallmark_Gene_Sets",
    #"Jensen_DISEASES"
  ))

  
  gc()
  return(current_results)
}



get_go_genes <- function(
    term_ids,
    library_names = c("GO_Biological_Process_2023",
                      "GO_Molecular_Function_2023")
) {
  
  final_list <- list()
  
  # --- Сначала скачиваем все библиотеки ---
  libraries_data <- list()
  
  for (lib in library_names) {
    message(paste0("Downloading library: ", lib, "..."))
    
    api_url <- paste0(
      "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=",
      lib
    )
    
    con <- url(api_url)
    lines <- tryCatch(readLines(con), error = function(e) NULL)
    close(con)
    
    if (!is.null(lines)) {
      libraries_data[[lib]] <- lines
    } else {
      message(sprintf("  Не удалось скачать библиотеку %s. Пропускаю.", lib))
    }
  }
  
  # --- Теперь идём по term_ids (ВАЖНО: сохраняем порядок) ---
  for (id in term_ids) {
    
    found <- FALSE
    
    for (lib in names(libraries_data)) {
      
      lines <- libraries_data[[lib]]
      pattern <- paste0("\\(", id, "\\)")
      target_line <- grep(pattern, lines, value = TRUE)
      
      if (length(target_line) > 0) {
        
        parts <- strsplit(target_line[1], "\t")[[1]]
        
        term_full_name <- parts[1]
        genes <- parts[3:length(parts)]
        genes <- genes[genes != ""]
        n_genes <- length(genes)
        
        list_name <- paste0(term_full_name, " [", n_genes, " genes]")
        
        final_list[[list_name]] <- genes
        message(sprintf("  Добавлено: %s", list_name))
        
        found <- TRUE
        break  # ← прекращаем поиск в других библиотеках
      }
    }
    
    if (!found) {
      message(sprintf("  ID %s не найден ни в одной библиотеке.", id))
    }
  }
  
  gc()
  
  if (length(final_list) == 0) {
    warning("Ни один из ID не найден.")
  }
  
  return(final_list)
}



run_escape <- function(seu, gene_sets, method = 'UCell', group_by = "ident", groups=1000) {
  
  message(paste0("Running escape with method: ", method))
  
  # 1. Запуск escape
  seu <- runEscape(seu, 
                   method = method,
                   gene.sets = gene_sets, 
                   groups = groups, 
                   min.size = 3,
                   new.assay.name = method)
  
  # 2. Агрегация данных
  # Важно: используем слой 'data', так как UCell пишет туда
  avg_exp <- AverageExpression(seu, assays = method, group.by = group_by, slot = "data")
  pathway_avg <- as.matrix(avg_exp[[method]])
  
  # --- ЛОГИКА СОПОСТАВЛЕНИЯ ИМЕН ---
  # Печатаем для диагностики (потом можно убрать)
  message("Названия путей в расчетах: ", paste(head(rownames(pathway_avg), 2), collapse = ", "))
  message("Ожидаемые названия: ", paste(head(names(gene_sets), 2), collapse = ", "))
  
  # Вместо жесткого сопоставления, попробуем найти пересечение
  common_names <- intersect(rownames(pathway_avg), names(gene_sets))
  
  if (length(common_names) == 0) {
    # Если точного совпадения нет, возможно escape убрал спецсимволы.
    # Попробуем сопоставить по порядку (если количество совпадает)
    if (nrow(pathway_avg) == length(gene_sets)) {
      message("Точных совпадений имен нет, но количество путей совпадает. Использую оригинальные имена.")
      rownames(pathway_avg) <- names(gene_sets)
      existing_order <- names(gene_sets)
    } else {
      stop("Имена путей в ассае не совпадают с оригинальными именами списка.")
    }
  } else {
    existing_order <- names(gene_sets)[names(gene_sets) %in% common_names]
  }
  
  pathway_avg <- pathway_avg[existing_order, , drop = FALSE]
  # ------------------------------
  
  # 3. Z-score (с защитой от NaN)
  pathway_z <- t(scale(t(pathway_avg)))
  pathway_z[is.na(pathway_z)] <- 0 # Убираем NaN если они возникли
  
  # 4. Отрисовка
  ht <- Heatmap(pathway_z,
                name = "z-score",
                cluster_rows = FALSE,      
                cluster_columns = FALSE,   
                row_names_gp = gpar(fontsize = 12),
                column_names_gp = gpar(fontsize = 12, rotation = 90),  
                rect_gp = gpar(col = "black", lwd = .7),
                column_names_side = "top",
                row_names_side = 'left',
                width = ncol(pathway_z)*unit(8, "mm"),
                height = nrow(pathway_z)*unit(8, "mm"),
                
                
                col = viridis::inferno(100))
  
  return(list(obj = seu, plot = ht))
}



plot_enrichment_three <- function(enrich_list,
                                  padj_cutoff = 0.05,
                                  max_terms = 12,
                                  title = "Functional enrichment") {
  
  db_map <- c(
    GO_Biological_Process_2023 = "GO Biological Process",
    GO_Molecular_Function_2023 = "GO Molecular Function"
  )
  
  prepare_df <- function(df, db_name) {
    
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    df$Term <- as.character(df$Term)
    
    df %>%
      filter(Adjusted.P.value < padj_cutoff) %>%
      arrange(Adjusted.P.value) %>%
      slice_head(n = max_terms) %>%
      mutate(Database = db_name,
             logP = -log10(Adjusted.P.value),
             Term = str_trunc(Term, 70))
  }
  
  df_all <- bind_rows(
    prepare_df(enrich_list$GO_Biological_Process_2023, db_map[1]),
    prepare_df(enrich_list$GO_Molecular_Function_2023, db_map[2])
  )
  
  if (is.null(df_all) || nrow(df_all) == 0)
    stop("No significant enrichment terms found.")
  
  # порядок терминов внутри каждой панели
  df_all <- df_all %>%
    group_by(Database) %>%
    mutate(Term = factor(Term, levels = rev(unique(Term)))) %>%
    ungroup()
  
  ggplot(df_all,
         aes(x = logP,
             y = Term)) +
    geom_point(aes(size = Combined.Score,
                   color = logP)) +
    scale_color_viridis_c(option = "magma") +
    facet_wrap(~Database, scales = "free_y", ncol = 1) +
    theme_minimal(base_size = 15) +
    labs(
      x = expression(-log[10]("Adjusted P-value")),
      y = "",
      size = "Combined score",
      color = expression(-log[10]("FDR")),
      title = title
    ) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    )
}