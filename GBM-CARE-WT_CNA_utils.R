library(ggpubr)
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
# Reference correction
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################

.refrange2 <- function (cna, ...) {
  dots = list(...)
  v = sapply(dots, function(ref) rowMeans(cna[, ref, drop = F]), 
             simplify = F)
  res = list(Min = do.call(pmin, v), Max = do.call(pmax, v))
  res
}

.refcenter2 <- function (v, Min, Max, noise = NULL, zero_inflate = T) {
  if (!is.null(noise)) {
    Min = Min - noise
    Max = Max + noise
  }
  above = v > Max
  below = v < Min
  normal = !above & !below
  
  v[above] <- v[above] - Max
  v[below] <- v[below] - Min
  
  if(zero_inflate)
    v[normal] <- 0
  else
    v[normal] <- v[normal] - mean(c(Max, Min))
  
  v
}

refCorrect2 <- function (cna, noise = NULL, isLog = FALSE, zero_inflate = T, ...) {
  dots = list(...)
  genes = rownames(cna)
  # Args = c(list(cna = cna, isLog = isLog), dots)
  Args = c(list(cna = cna), dots)
  rg = do.call(.refrange2, Args)
  n = nrow(cna)
  cna = t(sapply(1:n, function(i) {
    .refcenter2(cna[i, ], Min = rg$Min[i], Max = rg$Max[i], noise = noise, zero_inflate = zero_inflate)
  }))
  rownames(cna) = genes
  cna
}

#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
# Plotting functions
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################

plot_cna_matrix <- function(cna_matrix, genome = "hg38", title = NULL, row_facet = NULL) {
  
  infercna::useGenome(genome)
  
  cg <- infercna::splitGenes(colnames(cna_matrix))
  cna_matrix <- cna_matrix[, unname(unlist(cg))]
  cl <- (infercna::splitGenes(colnames(cna_matrix)) %>% lengths)
  cl <- cl[cl > 0]
  
  # Превращаем в обычный data.frame, чтобы melt не ругался
  cna_df <- as.data.frame(cna_matrix)
  cna_df$CellID <- rownames(cna_matrix)
  
  # Явно вызываем melt из reshape2
  dm <- reshape2::melt(data = cna_df,
                       id.vars = c("CellID"),
                       variable.name = "Gene")
  
  names(cg) <- paste0("CHR", names(cg))
  gene2chr <- sapply(names(cg), function(chr) setNames(rep(chr, length(cg[[chr]])), cg[[chr]]), USE.NAMES = F)
  gene2chr <- unlist(gene2chr, recursive = F)
  
  dm$CHR <- gene2chr[as.character(dm$Gene)]
  dm$CHR <- factor(dm$CHR, levels = unique(dm$CHR))
  
  if(!is.null(row_facet)) {
    dm <- dm %>%
      left_join(row_facet, by = "CellID")
  }
  
  p <- ggplot(data = dm, mapping = aes(x = Gene, y = CellID)) +
    geom_raster(mapping = aes(fill = value))
  
  if(is.null(row_facet)) {
    p <- p + facet_grid(cols = vars(CHR), scales = "free")
  } else {
    p <- p + facet_grid(rows = vars(V1), cols = vars(CHR), scales = "free")
  }
  
  p <- p +
    scale_fill_gradient2(name = "Inferred CNA\n[log2 ratio]", low = "dodgerblue", mid = "white", high = "red", limits = c(-1, 1), oob = scales::squish) +
    xlab("Genes") +
    ylab("Cells") +
    ggtitle(title) +
    theme_pubr() +
    theme(axis.title = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
          legend.title = element_text(size = 16), legend.text = element_text(size = 12, angle = 90), legend.position = "none",
          plot.title = element_text(size = 20, face = "bold"),
          strip.text = element_text(size = 12)) +
    theme(panel.spacing = unit(.0, "lines"),
          panel.border = element_rect(color = "black", fill = NA), 
          strip.background = element_rect(color = "black"))
  
  return(p)
}

plot_cna_matrix_safe <- function(cna_matrix, genome = "hg38", title = NULL, row_facet = NULL) {
  infercna::useGenome(genome)
  cg <- infercna::splitGenes(colnames(cna_matrix))
  cna_matrix <- cna_matrix[, unname(unlist(cg))]
  
  cna_df <- as.data.frame(cna_matrix)
  cna_df$CellID <- rownames(cna_matrix)
  
  dm <- reshape2::melt(data = cna_df, id.vars = c("CellID"), variable.name = "Gene")
  
  gene2chr <- unlist(sapply(names(cg), function(chr) setNames(rep(paste0("CHR", chr), length(cg[[chr]])), cg[[chr]]), USE.NAMES = F))
  dm$CHR <- gene2chr[as.character(dm$Gene)]
  
  # КРИТИЧЕСКИЙ МОМЕНТ: Убираем гены, которые не привязались к хромосомам
  dm <- dm[!is.na(dm$CHR), ]
  dm$CHR <- factor(dm$CHR, levels = paste0("CHR", c(1:22, "X")))
  
  if(!is.null(row_facet)) {
    dm <- dm %>% left_join(row_facet, by = "CellID")
    dm <- dm[!is.na(dm$V1), ] # Убираем клетки без аннотации
  }
  
  ggplot(data = dm, mapping = aes(x = Gene, y = CellID)) +
    geom_raster(mapping = aes(fill = value)) +
    facet_grid(rows = vars(V1), cols = vars(CHR), scales = "free", space = "free") +
    scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red", limits = c(-0.5, 0.5), oob = scales::squish) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.spacing = unit(0, "lines")) +
    ggtitle(title)
}

plot_cna_sig_per_chr <- function(cna_sig, event_stats, q025, q975, title = NULL) {
  
  ggplot(cna_sig, aes(x = CHR, y = Sig, color = isReference)) +
    geom_boxplot(aes(fill = isReference), outlier.shape = NA, color = "black") +
    geom_text(data = data.frame(CHR = factor(names(event_stats), levels = c(paste0("CHR", 1:22), "CHRX")), Stat = paste0(round(event_stats, 0), "%")), mapping = aes(label = Stat, x = CHR, y = .25), size = 8, inherit.aes = F) +
    facet_grid(cols = vars(CHR), scales = "free_x") +
    scale_color_discrete(name = "") +
    geom_hline(data = data.frame(CHR = factor(names(q975), levels = c(paste0("CHR", 1:22), "CHRX")), Q975 = q975), mapping = aes(yintercept = Q975), linetype = "dashed", size = 1, color = "black") +
    geom_hline(data = data.frame(CHR = factor(names(q025), levels = c(paste0("CHR", 1:22), "CHRX")), Q025 = q025), mapping = aes(yintercept = Q025), linetype = "dashed", size = 1, color = "black") +
    theme_pubr() +
    scale_y_continuous(limits = c(-.5, .5)) +
    xlab("Chromosome") +
    ylab("CNA signal") +
    ggtitle(title) +
    theme(axis.title = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 16),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          legend.title = element_text(size = 20), legend.text = element_text(size = 16),
          plot.title = element_text(size = 20, face = "bold"),
          strip.text = element_text(size = 16)) +
    theme(panel.spacing = unit(.0, "lines"),
          panel.border = element_rect(color = "black", fill = NA),
          strip.background = element_rect(color = "black"))
}

plot_cna_sig_cor_dist <- function(cna_scores, norm_fit_sig, norm_fit_cor) {
  
  ggplot(data = cna_scores, aes(x = CNAsig, fill = isReference, y = after_stat(ndensity))) +
    geom_histogram(color = "black", alpha = .5, bins = 100) +
    # this plots the the data of the normal distribution that we fitted
    geom_density(data = data.frame(x = norm_fit_sig), mapping = aes(x = x, y = after_stat(ndensity)), inherit.aes = F, size = 2) +
    geom_vline(xintercept = quantile(norm_fit_sig, .95), linetype = "dashed", size = 2) +
    xlab("CNA signal") +
    ylab("Density (scales to 1)") + 
    theme_classic() +
    theme(axis.text = element_text(size=20), axis.title = element_text(size=20),
          plot.title = element_text(size=24), legend.text = element_text(size=16), legend.title = element_text(size = 20),
          legend.key=element_blank(), panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"), aspect.ratio = 1, legend.position = "none") +
    ggplot(data = cna_scores, aes(x = CNAcor, fill = isReference)) +
    geom_histogram(mapping = aes(y = after_stat(ndensity)), color = "black", alpha = .5, bins = 100) +
    geom_density(data = data.frame(x = norm_fit_cor), mapping = aes(x = x, y = after_stat(ndensity)), inherit.aes = F, size = 2) +
    geom_vline(xintercept = quantile(norm_fit_cor, .95), linetype = "dashed", size = 2) +
    xlab("CNA correlation") +
    ylab("Density (scales to 1)") + 
    theme_classic() +
    theme(axis.text = element_text(size=20), axis.title = element_text(size=20),
          plot.title = element_text(size=24), legend.text = element_text(size=16), legend.title = element_text(size = 20),
          legend.key=element_blank(), panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"), aspect.ratio = 1, legend.position = "none")
}

plot_step_1_result <- function(cna_scores) {
  
  ggscatter(data = cna_scores, x = "UMAP1", y = "UMAP2", color = "Cluster") +
    scale_color_discrete("Cluster") +
    theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16),
          legend.title = element_text(size = 20), legend.text = element_text(size = 16),
          plot.title = element_text(size = 20, face = "bold"), legend.position = "right",
          strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 16)) +
    ggscatter(data = cna_scores, x = "UMAP1", y = "UMAP2", color = "CellType") +
    scale_color_discrete("Cell Type") +
    theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16),
          legend.title = element_text(size = 20), legend.text = element_text(size = 16),
          plot.title = element_text(size = 20, face = "bold"), legend.position = "right",
          strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 16)) +
    ggscatter(data = cna_scores, x = "CNAcor", y = "CNAsig", color = "CNAclassification", xlab = "CNA correlation", ylab = "CNA signal") +
    scale_color_discrete("Classification") +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16),
          legend.title = element_text(size = 20), legend.text = element_text(size = 16),
          plot.title = element_text(size = 20, face = "bold"), legend.position = "right",
          strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 16)) +
    ggscatter(data = cna_scores, x = "UMAP1", y = "UMAP2", color = "CNAclassification") +
    scale_color_discrete("Classification") +
    theme(axis.title = element_text(size = 20), axis.text = element_text(size = 16),
          legend.title = element_text(size = 20), legend.text = element_text(size = 16),
          plot.title = element_text(size = 20, face = "bold"), legend.position = "right",
          strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 16))
}
