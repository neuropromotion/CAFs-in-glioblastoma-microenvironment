# Dataset URL:
# https://www.10xgenomics.com/datasets/human-glioblastoma-multiforme-3-v-3-whole-transcriptome-analysis-3-standard-4-0-0

library(Seurat)
library(DoubletFinder)
feature_plot <- function(data, feature){
  FeaturePlot(data, feature, order=T, min.cutoff = 0, max.cutoff = 4)  +
    scale_color_gradientn(
      colours = c("gray", "blue1", "deeppink1"),   
      values = c(0, 0.2, .65,  1),  
      limits = c(0, 4)  
    )
}
count_plots <- function(obj, 
                        feature_lower=600,
                        feature_upper=6000,
                        count_lower=200,
                        count_upper=1e5,
                        mito=15){
  
  vln_plot1 <- VlnPlot(obj, features = "nFeature_raw", layer = "counts") + 
    theme(panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),       # Убрать подпись оси x
          axis.text.x = element_blank(),        # Убрать метки оси x
          axis.ticks.x = element_blank()) + 
    geom_hline(yintercept = feature_upper, linetype = "dashed", color = "red", cex=1)+ 
    geom_hline(yintercept = feature_lower, linetype = "dashed", color = "red", cex=1)+ 
    NoLegend()
  
  vln_plot2 <- VlnPlot(obj, features = "nCount_raw", layer = "counts") + 
    theme(panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),       # Убрать подпись оси x
          axis.text.x = element_blank(),        # Убрать метки оси x
          axis.ticks.x = element_blank()) + 
    geom_hline(yintercept = count_upper, linetype = "dashed", color = "red", cex=1)+ 
    geom_hline(yintercept = count_lower, linetype = "dashed", color = "red", cex=1)+ 
    NoLegend()
  
  vln_plot3 <- VlnPlot(obj, features = "percent.mt", layer = "counts") + 
    theme(panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),       # Убрать подпись оси x
          axis.text.x = element_blank(),        # Убрать метки оси x
          axis.ticks.x = element_blank()) + 
    geom_hline(yintercept = mito, linetype = "dashed", color = "red", cex=1) +
    NoLegend()
  combined_plot <- vln_plot1 | vln_plot2 | vln_plot3
  return(combined_plot)
}
deduplex <- function(obj){
  nExp <- round(0.05 * ncol(obj))
  pK <- 0.09 
  obj <- doubletFinder(obj, PCs = 1:10, pN = 0.25, pK = 0.09, 
                       nExp = nExp, reuse.pANN = NULL, sct = FALSE)
  
  df_col <- grep("DF.classifications", colnames(obj@meta.data), value = TRUE)
  cat('Cells Before deduplex: ', ncol(obj), '\n')
  obj <- subset(
    obj,
    subset = !!rlang::sym(df_col) == "Singlet"
  )
  cat('Cells Before deduplex: ', ncol(obj), '\n')
  return(obj)
} # DoubletFinder
get_chromosome_means <- function(counts, path_to_mapped_genes='path_to/mapped_genes.txt'){
  row_sums <- rowSums(counts)
  counts <- counts[row_sums > 0, ]
  norm_counts <- t(t(counts) / colSums(counts) * 1000)  
  log_counts <- log2(norm_counts + 1) 
  gene_filter <- rowSums(log_counts) >= 100
  log_counts_filtered <- log_counts[gene_filter, ]
  # Remove HLA-genes
  hla_genes <- grep("HLA", rownames(log_counts_filtered), value = TRUE)
  log_counts_filtered <- log_counts_filtered[!(rownames(log_counts_filtered) %in% hla_genes), ]
  gene_chromosome_map <- read.table(path_to_mapped_genes, sep = '\t', header = T)  # Файл с генами и их хромосомами
  chromosome_means <- sapply(unique(gene_chromosome_map$Chromosome.scaffold.name), function(chrom) {
    genes <- intersect(
      gene_chromosome_map$HGNC.symbol[gene_chromosome_map$Chromosome.scaffold.name == chrom],
      rownames(log_counts_filtered)  
    )
    
    if (length(genes) > 0) {
      colMeans(log_counts_filtered[genes, , drop = FALSE], na.rm = TRUE)
    } else {
      rep(NA, ncol(log_counts_filtered))  # Если генов нет, возвращаем NA
    }
  })
  chromosome_means <- as.data.frame(chromosome_means)
  colnames(chromosome_means) <- paste0("Chr", rep(1:22))
  rownames(chromosome_means) <- rownames(t(log_counts_filtered))
  return(chromosome_means)
}

# load and save RDS
#saveRDS(val, file = "/mnt/jack-5/amismailov/CAF_study/validation_2/val2.rds")
#val <- readRDS(file = "/mnt/jack-5/amismailov/CAF_study/validation_2/val2.rds")


data.dir <- "path_to_dir"
m <- Read10X(data.dir = data.dir)
val <- CreateSeuratObject(counts = m, project = "GBM", min.cells = 3, min.features = 200)

cts <- val[["RNA"]]$counts
val$nFeature_raw <- Matrix::colSums(cts > 0)   # Кол-во детектированных генов
val$nCount_raw   <- Matrix::colSums(cts)       # Суммарные UMI
val$percent.mt   <- PercentageFeatureSet(val, pattern = "^MT-")  # Митохондрии

count_plots(val, feature_lower=300, feature_upper=8000, count_upper=9e4)

#24381 x 5960
val <- subset(
  val,
  subset = nFeature_raw >= 300 & nFeature_raw <= 8000 & nCount_raw <= 9e4 & percent.mt <= 15
)
#24381 x 5398

val <- NormalizeData(val)                             # LogNormalize
val <- FindVariableFeatures(val, selection.method = "vst", nfeatures = 2000)    # HVG
val <- ScaleData(val, features = rownames(val)) # регрессия по mt
val <- RunPCA(val, features = VariableFeatures(object = val))
ElbowPlot(val)
val <- FindNeighbors(val, dims = 1:15)
val <- FindClusters(val, resolution = 0.1)
val <- RunUMAP(val, dims = 1:15)
val <- deduplex(val)
DimPlot(val, label=T)

# Then GSVA annotation adding (GSVA.R script required)
 
# CELL_CHAT _____________________________
library(CellChat)
library(openxlsx)

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
  interaction_input <- read.csv(file = 'path_to_dir/interaction_input_CellChatDB.csv')
  complex_input <- read.csv(file = 'path_to_dir/complex_input_CellChatDB.csv', row.names = 1)
  cofactor_input <- read.csv(file = 'path_to_dir/cofactor_input_CellChatDB.csv', row.names = 1)
  geneInfo <- read.csv(file = 'path_to_dir/geneInfo_CellChatDB.csv', row.names = 1)
  #View(interaction_input)
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
  
  
  unique(cellchat@idents)
  dim(cellchat@data)
  # извлечение генов, которые есть только в базе данных взаимодействия
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  #future::plan("multisession", workers = 4)
  # Идентифицировать сверхэкспрессию лигандов/рецепторов в одном типе клеток
  #options(future.globals.maxSize = 1024 * 1024^2)  # 1 GiB ЕСЛИ ОШИБКА БЕЗ
  cellchat <- identifyOverExpressedGenes(cellchat)
  
  # Выявление взаимных пар сверхэкспрессии
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Сглаженные значения экспрессии (предназначены для устранения эффекта выпадения, по желанию не используются)
  # We also provide a function to project gene expression data onto protein-protein interaction (PPI) network. Specifically, a diffusion process is used to smooth genes’ expression values based on their neighbors’ defined in a high-confidence experimentally validated protein-protein network. This function is useful when analyzing single-cell data with shallow sequencing depth because the projection reduces the dropout effects of signaling genes, in particular for possible zero expression of subunits of ligands/receptors. 
  cellchat <- projectData(cellchat, PPI.human)
  # Расчет возможностей взаимодействия
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)  
  
  
  # Пары взаимодействия с низким процентом клеток, экспрессирующих фильтр
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat) # Агрегируйте сети коммуникаций
  return(cellchat)
}
cellchat <- cellchat_workflow(val)

print(cellchat@netP$pathways) # все пары взаимодействий

cafs_pairs <- data.frame(interaction_name = c('SPP1_ITGAV_ITGB5',
                                              'SPP1_ITGA5_ITGB1',
                                              'SPP1_ITGA4_ITGB1',
                                              'MDK_ITGA6_ITGB1',
                                              'JAG1_NOTCH4.1',
                                              'IGFBP7_CD93',
                                              'ANXA1_FPR1',
                                              'ANGPT2_ITGA5_ITGB1',
                                              'PPIA_BSG',
                                              'PTN_PTPRZ1', 
                                             'MDK_PTPRZ1', 
                                             'PTN_NCL',
                                             'MDK_NCL',
                                             'MDK_LRP1',
                                             'GRN_SORT1',
                                             'HBEGF_EGFR', 
                                             'PDGFA_PDGFRA'))
svg("vln.svg", width = 3.9, height = 3.9) 
netVisual_bubble(cellchat, 
                 sources.use = c('CAFs'),
                 targets.use = c('T cells'),
                 #pairLR.use = cafs_pairs,
                 remove.isolate = T,
                 font.size=10,
                 n.colors = 11,
                 color.grid = "black",  
                 angle.x = 90) +  
  scale_color_gradientn(
    colors = c("black", 'cyan', "blue", "magenta"),
    limits = c(0, 1),
    values = scales::rescale(c(0, 0.1, 0.4, 1)))
dev.off() 

#### All communications with related information:
comm_df <- subsetCommunication(cellchat)   
comm_df_pathway <- subsetCommunication(cellchat, slot.name = "netP")
write.xlsx(comm_df, file = "CellChat_Communications_validation_2.xlsx", sheetName = "Ligand_Receptor", rowNames = FALSE)
#### #### #### #### #### #### #### #### #### 


#--------------------------CNV----------------------------------------------------------------
counts <- GetAssayData(val, assay = "RNA", layer = "counts")
chromosome_means <- get_chromosome_means(counts)

write.csv(
  data.frame(chromosome_means),
  file = "chromosome_means_validation_v2.csv",
  row.names = TRUE,            # не пишем номер строки
  quote = FALSE                 # без кавычек, если не нужны
)



####### ML RES:
res <- read_csv('result_ML_val2.csv')
tail(res)
res_df <- as.data.frame(res)
head(res_df)

ann <- ifelse(res_df$pred == 1, "Glioblastoma cells", "Stromal cells")
names(ann) <- res_df$label
val$ML_annotation <- ann[Cells(val)]
gcf
#View(merged_filtered@meta.data)



svg('vln.svg', width = 6, height = 5)
DimPlot(val, group.by = 'ML_annotation', cols = c('#FF3E96', '#1E90FF'))
dev.off()

DimPlot(val, group.by = 'init_annot', cols=colors_initial)
 

