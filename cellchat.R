library(CellChat)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(viridis)
library(stringr)
library(ggrepel)
#-------------------------------------------------------------------------------------

cellchat_workflow <- function(gcf){
  options(stringsAsFactors = FALSE)
  gcf <- JoinLayers(object = gcf)
  data.input <- GetAssayData(gcf, assay = "RNA", layer = "data") # normalized data matrix
  labels <- Idents(gcf)
  meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
  cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
  cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  interaction_input <- read.csv(file = '/cellchat_config/interaction_input_CellChatDB.csv')
  complex_input <- read.csv(file = '/cellchat_config/complex_input_CellChatDB.csv', row.names = 1)
  cofactor_input <- read.csv(file = '/cellchat_config/cofactor_input_CellChatDB.csv', row.names = 1)
  geneInfo <- read.csv(file = '/cellchat_config/geneInfo_CellChatDB.csv', row.names = 1)

  rownames(interaction_input) <- make.unique(as.character(interaction_input$Unnamed..0))
  interaction_input$Unnamed..0 <- NULL
  interaction_input$X <- NULL
  
  CellChatDB <- list()
  CellChatDB$interaction <- interaction_input
  CellChatDB$complex <- complex_input
  CellChatDB$cofactor <- cofactor_input
  CellChatDB$geneInfo <- geneInfo
  
  CellChatDB.ss <- subsetDB(CellChatDB, search = "Secreted Signaling", key='annotation')
  cellchat@DB <- CellChatDB.ss
  
  
  cellchat <- subsetData(cellchat) 
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)  

  cellchat <- filterCommunication(cellchat, min.cells = 1)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat) # Агрегируйте сети коммуникаций
  return(cellchat)
}

bubble_wrapper <- function(object, 
                           sources = NULL, 
                           targets = NULL, 
                           pairLR_use = NULL, 
                           color_option = "plasma", 
                           base_size = 10,
                           stroke_width = 0.5,
                           title = '') {
  
  # 1. Определяем источники (если NULL — берем все)
  if (is.null(sources)) {
    sources.use <- levels(object@idents)
  } else {
    sources.use <- sources
  }
  
  # 2. Извлекаем данные
  suppressMessages({
    res <- netVisual_bubble(object, 
                            sources.use = sources.use, 
                            targets.use = targets, 
                            remove.isolate = TRUE, 
                            return.data = TRUE)
  })
  
  df <- res$communication
  
  if (is.null(df) || nrow(df) == 0) {
    message("Данные не найдены.")
    return(NULL)
  }
  
  df$source_target <- paste0(df$source, " -> ", df$target)
  
  if (!is.null(pairLR_use)) {
    df <- df[df$interaction_name %in% pairLR_use | df$interaction_name_2 %in% pairLR_use, ]
  }
  
  p <- ggplot(df, aes(x = source_target, y = interaction_name_2)) +
    geom_point(aes(size = -log10(pval), fill = prob), 
               shape = 21, 
               color = "black", 
               stroke = stroke_width) +
    
    scale_fill_viridis_c(option = color_option, 
                         name = "Communication\nprobability", 
                         oob = scales::squish) +
    
    scale_size_continuous(range = c(1, 6), name = "-log10(p-value)") +
    
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, h = 1, v = 1, color = "black"), # 45 градусов для читаемости
      axis.text.y = element_text(color = "black"),
      axis.title = element_blank(),
      panel.grid.major = element_line(color = "grey95"),
      legend.key.height = unit(0.8, "cm")
    ) +
    labs(title = title)
  
  return(p)
}


bubble_wrapper_v2 <- function(object, 
                           sources, 
                           targets = NULL, 
                           pairLR_use = NULL,    # Добавляем этот аргумент
                           color_option = "plasma", 
                           base_size = 10,
                           stroke_width = 0.5,
                           title='') {
  
  # 1. Извлекаем данные
  suppressMessages({
    res <- netVisual_bubble(object, 
                            sources.use = sources, 
                            targets.use = targets, 
                            remove.isolate = TRUE, 
                            return.data = TRUE)
  })
  
  df <- res$communication
  
  # 2. Фильтрация по списку взаимодействий
  if (!is.null(pairLR_use)) {
    # Фильтруем по колонке interaction_name (или interaction_name_2)
    df <- df[df$interaction_name %in% pairLR_use | df$interaction_name_2 %in% pairLR_use, ]
  }
  
  # Проверка после фильтрации
  if (is.null(df) || nrow(df) == 0) {
    message("Внимание: После фильтрации данных не осталось.")
    return(NULL)
  }
  
  # 3. Визуализация (остается прежней)
  p <- ggplot(df, aes(x = target, y = interaction_name_2)) +
    geom_point(aes(size = -log10(pval), fill = prob), 
               shape = 21, 
               color = "black", 
               stroke = stroke_width) +
    
    scale_fill_viridis_c(option = color_option, 
                         name = "Communication\nprobability", 
                         oob = scales::squish) +
    
    scale_size_continuous(range = c(1, 6), name = "-log10(p-value)") +
    
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 90, h = 1, v = 0.5, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title = element_blank(),
      panel.grid.major = element_line(color = "grey95"),
      legend.key.height = unit(0.8, "cm")
    ) +
    labs(title = title)
  
  return(p)
}


bubble_wrapper_v3 <- function(object, 
                              sources = NULL, 
                              targets = NULL, 
                              pairLR_use = NULL, 
                              color_option = "plasma", 
                              base_size = 11,
                              stroke_width = 0.5,
                              title = '') {
  
  suppressMessages({
    res <- netVisual_bubble(object, 
                            sources.use = sources, 
                            targets.use = targets, 
                            remove.isolate = TRUE, 
                            return.data = TRUE)
  })
  
  df <- res$communication
  if (is.null(df) || nrow(df) == 0) return(NULL)
  

  df$pval <- mapply(function(s, t, i) object@net$pval[s, t, i], 
                    as.character(df$source), 
                    as.character(df$target), 
                    as.character(df$interaction_name))
  
  if (!is.null(pairLR_use)) {
    df <- df[df$interaction_name %in% pairLR_use | df$interaction_name_2 %in% pairLR_use, ]
  }
  
  df$p_cat <- "p > 0.05"
  df$p_cat[df$pval <= 0.05] <- "0.01 < p < 0.05"
  df$p_cat[df$pval <= 0.01] <- "p < 0.01" 
  
  existing_levels <- intersect(c("p < 0.01", "0.01 < p < 0.05", "p > 0.05"), unique(df$p_cat))
  df$p_cat <- factor(df$p_cat, levels = existing_levels)
  
  df$interaction_pretty <- gsub(" - ", " \u2192 ", df$interaction_name_2)
  
  p <- ggplot(df, aes(x = target, y = interaction_pretty)) +
    geom_point(aes(size = p_cat, fill = prob), 
               shape = 21, color = "black", stroke = stroke_width) +
    
    scale_size_manual(values = c("p < 0.01" = 6, 
                                 "0.01 < p < 0.05" = 3.5, 
                                 "p > 0.05" = 1.5),
                      name = "Significance") +
    
    scale_fill_viridis_c(option = color_option, 
                         name = "Communication\nProbability") +
    
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, h = 1, v = 1, color = "black", face = "bold"),
      axis.text.y = element_text(color = "black", face = "bold"), 
      axis.title = element_blank(),
      panel.grid.major = element_line(color = "grey95"),
      legend.title = element_text(face = "bold")
    ) +
    labs(title = title)
  
  return(p)
}


bubble_wrapper_HR <- function(object, 
                              sources = NULL, 
                              targets = NULL, 
                              pairLR_use = NULL, 
                              hr_data = NULL,  
                              color_option = "viridis", 
                              base_size = 11,
                              stroke_width = 0.5,
                              title = '') {
  
  clean_str <- function(x) {
    toupper(gsub("[^[:alnum:]]", "", x))
  }
  
  suppressMessages({
    res <- netVisual_bubble(object, 
                            sources.use = sources, 
                            targets.use = targets, 
                            remove.isolate = TRUE, 
                            return.data = TRUE)
  })
  
  df <- res$communication
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df$pval <- mapply(function(s, t, i) object@net$pval[s, t, i], 
                    as.character(df$source), 
                    as.character(df$target), 
                    as.character(df$interaction_name))
  
  if (!is.null(hr_data)) {
    hr_lookup <- setNames(as.numeric(hr_data), clean_str(names(hr_data)))
    df$clean_key <- clean_str(df$interaction_name_2)
    df$hr_val <- hr_lookup[df$clean_key]
    df <- df[!is.na(df$hr_val), ]
  }
  
  if (!is.null(pairLR_use)) {
    allowed_keys <- clean_str(pairLR_use)
    df <- df[clean_str(df$interaction_name_2) %in% allowed_keys, ]
  }
  
  if (nrow(df) == 0) {
    message("Ошибка: Ни одна пара не совпала. Проверьте имена в hr_data и pairLR_use.")
    return(NULL)
  }
  
  df <- df %>%
    group_by(target, interaction_name_2) %>%
    summarize(
      pval = min(pval), 
      prob = max(prob), 
      hr_val = hr_val[1], 
      .groups = 'drop'
    )
  
  all_log_hrs <- log(as.numeric(hr_data))
  max_log_val <- max(all_log_hrs, na.rm = TRUE)
  
  df$hr_score <- log(df$hr_val) / max_log_val
  df$p_cat <- cut(df$pval, 
                  breaks = c(-Inf, 0.01, 0.05, Inf), 
                  labels = c("p < 0.01", "0.01 < p < 0.05", "p > 0.05"))
  
  df <- df %>% arrange(hr_val)
  df$interaction_pretty <- gsub(" - ", " \u2192 ", df$interaction_name_2)
  df$interaction_pretty <- factor(df$interaction_pretty, levels = unique(df$interaction_pretty))
  
  p <- ggplot(df, aes(x = target, y = interaction_pretty)) +
    geom_point(aes(size = p_cat, fill = hr_score), 
               shape = 21, color = "black", stroke = stroke_width) +
    
    scale_size_manual(values = c("p < 0.01" = 6, 
                                 "0.01 < p < 0.05" = 3.5, 
                                 "p > 0.05" = 1.5),
                      name = "Significance") +
    
    scale_fill_viridis_c(option = color_option, 
                         name = "Normalized\nlog(HR)",
                         limits = c(0, 1),
                         breaks = c(0, 0.5, 1)) +
    
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, h = 1, v = 1, color = "black", face = "bold"),
      axis.text.y = element_text(color = "black", face = "bold"), 
      axis.title = element_blank(),
      panel.grid.major = element_line(color = "grey95"),
      legend.title = element_text(face = "bold")
    ) +
    labs(title = title)
  
  return(p)
}

comm_prob_boxplot <- function(object, 
                                 sources = NULL, 
                                 targets = NULL, 
                                 pairLR_use = NULL, 
                                 color_option = "plasma", 
                                 base_size = 11,
                                 title = 'Communication Probability Distribution') {
  
  suppressMessages({
    res <- netVisual_bubble(object, 
                            sources.use = sources, 
                            targets.use = targets, 
                            remove.isolate = TRUE, 
                            return.data = TRUE)
  })
  
  df <- res$communication
  if (is.null(df) || nrow(df) == 0) return(NULL)
  if (!is.null(pairLR_use)) {
    df <- df[df$interaction_name %in% pairLR_use | df$interaction_name_2 %in% pairLR_use, ]
  }

  df$interaction_pretty <- gsub(" - ", " \u2192 ", df$interaction_name_2)
  
  library(dplyr)
  df <- df %>%
    group_by(interaction_pretty) %>%
    mutate(median_prob = median(prob)) %>%
    ungroup()
  
  # 5. ВИЗУАЛИЗАЦИЯ
  library(ggplot2)
  p <- ggplot(df, aes(x = interaction_pretty, y = prob)) + 
    geom_boxplot(aes(fill = median_prob), outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1, color = "black") +
    scale_fill_viridis_c(option = color_option, name = "Median\nComm. Prob", limits = c(0, 1)) +
    
    coord_flip() +
    theme_bw(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black", face = "bold"),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 9, face = "bold")
    ) +
    labs(
      title = title,
      x = "Ligand-Receptor Interaction",
      y = "Communication Probability"
    )
  
  return(p)
}



comm_prob_histogram_v3 <- function(object, 
                                        sources = NULL, 
                                        targets = NULL, 
                                        pairLR_use = NULL, 
                                        bins = 30,
                                        fill_color = "steelblue",
                                        base_size = 12,
                                        title = 'Distribution of Communication Probabilities') {
  
  suppressMessages({
    res <- netVisual_bubble(object, 
                            sources.use = sources, 
                            targets.use = targets, 
                            remove.isolate = TRUE, 
                            return.data = TRUE)
  })
  
  df <- res$communication
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  if (!is.null(pairLR_use)) {
    df <- df[df$interaction_name %in% pairLR_use | df$interaction_name_2 %in% pairLR_use, ]
  }
    
  breaks <- seq(0, max(df$prob), length.out = bins + 1)
  df$bin <- cut(df$prob, breaks = breaks, include.lowest = TRUE)
   
  df_summary <- df %>%
    group_by(bin) %>%
    summarise(
      count = n(), 
      labels = paste(head(unique(interaction_name_2), 3), collapse = "\n"),
      .groups = 'drop'
    ) %>%
    filter(count > 0) 
  
  bin_levels <- levels(df$bin)
  df_summary$bin_center <- sapply(as.character(df_summary$bin), function(x) { 
    bounds <- as.numeric(gsub("[^0-9.]", "", unlist(strsplit(x, ","))))
    mean(bounds)
  })
   
  median_val <- median(df$prob)
  p <- ggplot(df, aes(x = prob)) + 
    geom_histogram(breaks = breaks, fill = fill_color, color = "white", alpha = 0.7) +
     
    
    #geom_vline(xintercept = median_val, color = "red", linetype = "dashed") +
    
    theme_bw(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = title,
      x = "Communication Probability",
      y = "Count"
    )+ coord_cartesian(ylim = c(0, max(df_summary$count) * 1.3))
  
  # 4. ТАБЛИЦА (надежный dplyr::select)
  final_table <- df %>%
    dplyr::select(source, target, interaction_name_2, prob) %>%
    dplyr::arrange(desc(prob))
  
  return(list(plot = p, table = final_table))
}
