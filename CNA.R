library(Seurat)
library(readr)
library(dplyr)
library(MASS)
library(dplyr)
source('GBM-CARE-WT_CNA_utils.R')

get_cna_matrix <- function(seu, mapping_path='/home/amismailov/GBM_finder/genes_chr_mapping.csv'){
  expr <- GetAssayData(seu, layer = "data")
  
  mapping <- read_csv(mapping_path)
  mapping <- mapping %>%
    filter(Chromosome %in% c(as.character(1:22), "X"))
  
  common_genes <- intersect(rownames(expr), mapping$Gene_name)
  
  expr <- expr[common_genes, ]
  
  mapping <- mapping %>%
    filter(Gene_name %in% common_genes)
  
  
  # order
  mapping <- mapping[match(rownames(expr), mapping$Gene_name), ]
  
  rownames(expr) <- paste0(
    mapping$Gene_name, "|",
    mapping$Chromosome
  )
  
  cna_matrix <- as.matrix(expr)
  rm(expr)
  gc()
  return(cna_matrix)
}


get_hm <- function(seu, cna_smooth){
  seu@meta.data$CellID <- colnames(seu)
  cells_to_plot <- seu@meta.data %>%
    group_by(final_annotation) %>%
    slice_sample(n = 100) %>%
    pull(CellID) # Убедись, что колонка с баркодами называется CellID, или используй rownames()
  
  rownames(cna_smooth) <- gsub("\\|.*", "", rownames(cna_smooth))
  # 2. Фильтруем матрицу и аннотацию
  cna_matrix_small <- t(cna_smooth[, cells_to_plot])
  row_ann_small <- data.frame(
    CellID = cells_to_plot,
    V1 = seu$final_annotation[cells_to_plot]
  )
   
  p_final <- plot_cna_matrix_safe(
    cna_matrix = cna_matrix_small,
    row_facet = row_ann_small,
    title = "GBM Inferred CNA (Downsampled)"
  )
  
  gc()
  return(p_final)
}

final_prediction <- function(seu, cna_smooth, reference_cells, mapping_path='/home/amismailov/GBM_finder/genes_chr_mapping.csv'){
  mapping <- read_csv(mapping_path)
  mapping <- mapping %>%
    filter(Chromosome %in% c(as.character(1:22), "X"))
  
  # 1. Группируем гены по хромосомам (используем твой mapping)
  # Убедись, что имена генов в cna_smooth и mapping совпадают (мы убирали суффикс |CHR ранее)
  genes_by_chr <- split(mapping$Gene_name, mapping$Chromosome)
  
  # 2. Считаем средний сигнал для каждой хромосомы в каждой клетке
  chr_signals <- sapply(names(genes_by_chr), function(chr) {
    genes <- genes_by_chr[[chr]]
    genes_present <- intersect(genes, rownames(cna_smooth))
    if (length(genes_present) > 0) {
      return(colMeans(cna_smooth[genes_present, , drop = FALSE]))
    } else {
      return(rep(0, ncol(cna_smooth)))
    }
  }) # Результат: матрица [клетки x хромосомы]
  
  # 3. Фитируем нормальное распределение по референсным клеткам (Step 1 из статьи)
  # Для каждой хромосомы берем значения в reference_cells и считаем Mu и SD
  ref_stats <- apply(chr_signals[reference_cells, ], 2, function(x) {
    # fitdistr из MASS для нормального распределения
    fit <- MASS::fitdistr(x, "normal")
    return(c(mu = fit$estimate[["mean"]], sd = fit$estimate[["sd"]]))
  })
  
  # 4. Считаем P-values (Z-test) для каждой клетки на каждой хромосоме
  z_scores <- t(sapply(1:nrow(chr_signals), function(i) {
    (chr_signals[i, ] - ref_stats["mu", ]) / ref_stats["sd", ]
  }))
  rownames(z_scores) <- rownames(chr_signals)
  
  p_values <- 2 * pnorm(-abs(z_scores))
  
  # Коррекция Benjamini-Hochberg по каждой хромосоме
  adj_p_values <- apply(p_values, 2, p.adjust, method = "BH")
  
  # 5. Определяем события (Real и Suspicious)
  # Real: adj_p < 0.05
  # Suspicious: unadjusted p < 0.05
  has_real_event <- adj_p_values < 0.05
  has_suspicious_event <- p_values < 0.05
  
  # Специфические для GBM события (согласно статье)
  is_chr7_gain <- has_real_event[, "7"] & z_scores[, "7"] > 0
  is_chr10_loss <- has_real_event[, "10"] & z_scores[, "10"] < 0
  
  # Базовая классификация Step 1
  # Малигнантные - те, у кого есть хотя бы одно "реальное" событие
  seu$step1_class <- "nonmalignant"
  seu$step1_class[rowSums(has_real_event) >= 1] <- "malignant"
  seu$step1_class[rowSums(has_real_event) == 0 & rowSums(has_suspicious_event) >= 1] <- "unresolved"
  
  # Уточнение Step 2 (Correlation)
  seu$step2_class <- seu$step1_class
  seu$step2_class[seu$step1_class == "malignant" & seu$cna_cor < 0.5] <- "unresolved"
  seu$step2_class[seu$step1_class == "nonmalignant" & seu$cna_cor > 0.35] <- "unresolved"
  
  cor_threshold <- 0.4
  
  # --- ШАГ 2: Малигнантные по "уликам" ---
  # В GBM +7 и -10 настолько специфичны, что их наличия достаточно для уверенности
  is_sure_malignant <- (adj_p_values[, "7"] < 0.05 & z_scores[, "7"] > 0) | 
    (adj_p_values[, "10"] < 0.05 & z_scores[, "10"] < 0)
  
  # Добавим "подозрительные" клетки (unadjusted p < 0.01)
  is_suspected_malignant <- (p_values[, "7"] < 0.01 & z_scores[, "7"] > 0) | 
    (p_values[, "10"] < 0.01 & z_scores[, "10"] < 0)
  
  # --- ШАГ 3: Новая логика финального статуса ---
  seu$final_CNA_status_soft <- "Non-malignant"
  
  # Клетка Malignant, если:
  # ЛИБО у нее высокая корреляция (> cor_threshold)
  # ЛИБО у нее есть железное событие (+7/-10) даже при средней корреляции
  seu$final_CNA_status_soft[seu$cna_cor > cor_threshold | is_sure_malignant] <- "Malignant"
  
  # Клетка остается Unresolved только если она совсем противоречивая:
  # Например, высокая корреляция, но при этом она из "нормального" кластера (5, 8, 10)
  # Или наоборот - низкая корреляция, но сидит в опухолевом кластере.
  # Но для "мягкой" аннотации можно просто сделать бинарно:
  
  # Если хочешь совсем без Unresolved (чисто бинарно):
  seu$final_CNA_status_soft <- ifelse(
    seu$cna_cor > cor_threshold | is_sure_malignant, 
    "Malignant", 
    "Non-malignant"
  )
  
  return(seu)
}










