# Dataset URL:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173278
library(Seurat)
library(Matrix)
library(data.table) 
library(DoubletFinder)
library(ggplot2)

#saveRDS(gcf, file = "/mnt/jack-5/amismailov/CAF_study/GCF/gcf.rds")
gcf <- readRDS(file = "/mnt/jack-5/amismailov/CAF_study/GCF/gcf.rds")

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
feature_plot <- function(data, feature){
  FeaturePlot(data, feature, order=T, slot = 'data')  +
    scale_color_gradientn(
      colours = c("gray", "blue1", "deeppink1"),  
      values = c(0, 0.3, .65,  1),  
    )
}
count_plots <- function(obj, 
                        feature_lower=600,
                        feature_upper=6000,
                        count_lower=200,
                        count_upper=1e5,
                        mito=15,
                        pt.size=0){
  
  vln_plot1 <- VlnPlot(obj, features = "nFeature_raw", layer = "counts", pt.size=pt.size) + 
    theme(panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),       # Убрать подпись оси x
          axis.text.x = element_blank(),        # Убрать метки оси x
          axis.ticks.x = element_blank()) + 
    geom_hline(yintercept = feature_upper, linetype = "dashed", color = "red", cex=1)+ 
    geom_hline(yintercept = feature_lower, linetype = "dashed", color = "red", cex=1)+ 
    NoLegend()
  
  vln_plot2 <- VlnPlot(obj, features = "nCount_raw", layer = "counts", pt.size=pt.size) + 
    theme(panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),       # Убрать подпись оси x
          axis.text.x = element_blank(),        # Убрать метки оси x
          axis.ticks.x = element_blank()) + 
    geom_hline(yintercept = count_upper, linetype = "dashed", color = "red", cex=1)+ 
    geom_hline(yintercept = count_lower, linetype = "dashed", color = "red", cex=1)+ 
    NoLegend()
  
  vln_plot3 <- VlnPlot(obj, features = "percent.mt", layer = "counts", pt.size=pt.size) + 
    theme(panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),       # Убрать подпись оси x
          axis.text.x = element_blank(),        # Убрать метки оси x
          axis.ticks.x = element_blank()) + 
    geom_hline(yintercept = mito, linetype = "dashed", color = "red", cex=1) +
    NoLegend()
  combined_plot <- vln_plot1 | vln_plot2 | vln_plot3
  return(combined_plot)
}

path_gcf <- '/mnt/jack-5/amismailov/CAF_study/GCF/filtered'

gcf <- Read10X(data.dir = path_gcf)
gcf <- CreateSeuratObject(counts = gcf, project = "GCF", min.cells = 0, min.features = 0)

cts_gcf <- gcf[["RNA"]]$counts
gcf$nFeature_raw <- Matrix::colSums(cts_gcf > 0)   # Кол-во детектированных генов
gcf$nCount_raw   <- Matrix::colSums(cts_gcf)       # Суммарные UMI
gcf$percent.mt   <- PercentageFeatureSet(gcf, pattern = "^MT-")  # Митохондрии

count_plots(gcf, feature_lower=1000, feature_upper=5500, count_upper=5500, count_lower=3500)
# Before = 75075 cells x 23076 genes
gcf <- subset(
  gcf,
  subset = nFeature_raw >= 1000 & nFeature_raw <= 5500 & 
    nCount_raw <= 5500 & nCount_raw >= 3500 &percent.mt <= 15
)
# After = 68585 cells x 23076 genes


#--------------------------AMBIENT RNA REMOVAL (SoupX)----------------------------------------------------------------
library(SoupX)

raw_dir <- 'path_to_dir/raw'
filt_dir <- 'path_to_dir/filtered'
tod <- Read10X(data.dir = raw_dir)      # raw droplets
toc <- Read10X(data.dir = filt_dir)     # filtered cells

# create seurat object for filtered data:
seu <- CreateSeuratObject(counts = toc, assay = "RNA", min.features = 0, min.cells = 0)
cts <- seu[["RNA"]]$counts
seu$nFeature_raw <- Matrix::colSums(cts > 0)   # Кол-во детектированных генов
seu$nCount_raw   <- Matrix::colSums(cts)       # Суммарные UMI
seu$percent.mt   <- PercentageFeatureSet(seu, pattern = "^MT-")  # Митохондрии

seu <- subset(
  seu,
  subset = nFeature_raw >= 1000 & nFeature_raw <= 5500 & 
    nCount_raw <= 5500 & nCount_raw >= 3500 &percent.mt <= 15
)
toc_qc <- GetAssayData(seu, "RNA", "counts")

common_genes <- intersect(rownames(tod), rownames(toc_qc))

tod <- tod[common_genes, , drop = FALSE]
toc_qc <- toc_qc[common_genes, , drop = FALSE]



sc <- SoupChannel(tod, toc_qc)

get_soup_groups <- function(sobj){
  sobj <- NormalizeData(sobj, verbose = FALSE)
  sobj <- FindVariableFeatures(object = sobj, nfeatures = 2000, verbose = FALSE, selection.method = 'vst')
  sobj <- ScaleData(sobj, verbose = FALSE)
  sobj <- RunPCA(sobj, npcs = 20, verbose = FALSE)
  sobj <- FindNeighbors(sobj, dims = 1:20, verbose = FALSE)
  sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE)
  
  return(sobj@meta.data[['seurat_clusters']])
  
}
add_soup_groups <- function(sobj){
  sobj$soup_group <- get_soup_groups(sobj)
  return(sobj)
}

seu <- add_soup_groups(seu)

clusters <- seu$soup_group
names(clusters) <- colnames(seu)
sc <- setClusters(sc, clusters)

sp <- data.frame(
  est    = Matrix::rowSums(tod) / sum(tod),
  counts = Matrix::rowSums(tod)
)
sp <- sp[is.finite(sp$est) & sp$counts > 0, , drop = FALSE]
rownames(sp) <- rownames(tod)
sc <- setSoupProfile(sc, sp)

# 4) оценка и коррекция
sc <- autoEstCont(sc, doPlot = FALSE)               # при необходимости: maxCont = 0.2, priorRho = 0.05
adj <- adjustCounts(sc, roundToInt = TRUE)

seu[["RNA"]] <- SetAssayData(seu[["RNA"]], slot = "counts", new.data = adj)
gcf <- seu
#---------------main_pipeline-----------------------------------------------------------
 
gcf <- NormalizeData(gcf)                             
gcf <- FindVariableFeatures(gcf, selection.method = "vst", nfeatures = 2000)   
#plan(sequential) 
gcf <- ScaleData(gcf, features = rownames(gcf)) 
gcf <- RunPCA(gcf, features = VariableFeatures(object = gcf))
ElbowPlot(gcf)
gcf <- FindNeighbors(gcf, dims = 1:20)
gcf <- FindClusters(gcf, resolution = 0.01)
gcf <- RunUMAP(gcf, dims = 1:20)

DimPlot(gcf, label=T, group.by = 'final_annot')
feature_plot(gcf, 'PECAM') 

gcf <- deduplex(gcf)

#------------------------DEG----------------------------------------
 
 

feature_plot(gcf, 'PECAM1')
DimPlot(gcf)
FAM <- FindAllMarkers(gcf,
                      logfc.threshold = 0.25,
                      min.pct = 0.1,
                      only.pos = TRUE,
                      test.use = 'wilcox',
                      slot = 'data')
library(writexl)
write_xlsx(FAM, "/home/amismailov/FAM_validation_3.xlsx")

cluster.names <- c('CAFs', 'Endothelial cells','Macrophages/Microglia',
               'Oligodendrocytes','MES1-like', 'MES2-like', 'NPC1-like', 'NPC2-like', 'AC-like',
               'OPC-like', 'Glioblastoma Stem Cells') 
top.genes = list()
top.genes <- character()           # итоговый вектор

for (cl in cluster.names) {
  cluster_genes <- FAM$gene[FAM$cluster == cl]   # гены текущего кластера
  added <- 0                                     # сколько уже добавили из этого кластера
  idx   <- 1                                     # индекс по списку cluster_genes
  
  # добавляем ровно 5 уникальных генов
  while (added < 5 && idx <= length(cluster_genes)) {
    g <- cluster_genes[idx]
    if (!(g %in% top.genes)) {                   # берём только новый ген
      top.genes <- c(top.genes, g)
      added <- added + 1
    }
    idx <- idx + 1                               # переходим к следующему гену
  }
}

new_cluster_names <- c(
  '0' = 'GBM',
  '1' = 'CAFs',
  '2' = 'GBM',
  '3' = 'GBM',
  '4' = 'Macrophages/Microglia',
  '5' = 'GBM',
  '6' = 'GBM',
  '7' = 'Endothelial cells',
  '8' = 'Oligodendrocytes',
  '9' = 'GBM',
  '10' = 'GBM',
  '11' = 'GBM',
  '12' = 'GBM',
  '13' = 'GBM')  
new_cluster_names <- c('CAF-like GBM' = 'GBM')
gcf <- RenameIdents(gcf, new_cluster_names)
DimPlot(gcf)


DimPlot(gcf, label=T)

feature_plot(gcf, 'BEX1')

#--------------------------CAFS----------------------------------------------------------------
cafs <- WhichCells(gcf, idents = 'Endothelial cells')
cafs <- subset(gcf, cells = cafs)

cafs <- RunPCA(cafs, features = VariableFeatures(object = cafs))
cafs <- FindNeighbors(cafs, dims = 1:5)
cafs <- FindClusters(cafs, resolution = 0.18)
cafs <- RunUMAP(object = cafs, dims = 1:5)


DimPlot(cafs, label=T)
FAM_cafs <- FindAllMarkers(cafs,
                      logfc.threshold = 0.25,
                      min.pct = 0.1,
                      only.pos = TRUE,
                      test.use = 'wilcox',
                      slot = 'data')
top.genes = list()
set.of.top.genes = vector()

feature_plot(cafs, 'PECAM1')
feature_plot(gcf, 'RGS5')
DimPlot(cafs, label=T)

for (cl in 0:4) {
  cluster_genes <- FAM_cafs$gene[FAM_cafs$cluster == cl][1:30]   # гены текущего кластера
  curr_name <- as.character(cl)#paste0('cluster', as.character(cl))
  top.genes[[curr_name]] = cluster_genes
  set.of.top.genes = c(set.of.top.genes, cluster_genes)
}

new_cluster_names <- c(
  '0' = 'CAFs',
  '1' = 'Endothelial cells',
  '2' = 'Endothelial cells',
  '3' = 'CAFs',
  '4' = 'CAFs')  

cafs <- RenameIdents(cafs, new_cluster_names)


cafs_clusters <- Idents(cafs)
cafs_annotations <- cafs_clusters
gcf$celltype_refined <- as.character(Idents(gcf))
gcf$celltype_refined[colnames(cafs)] <- as.character(cafs_annotations)
Idents(gcf) <- gcf$celltype_refined


#-------------------------NEW_ORDER--------------------
current_idents <- Idents(gcf)
new.order <- c('CAFs', 'Endothelial cells','Macrophages/Microglia',
               'Oligodendrocytes','GBM') 
new.order <- c('CAFs', 'Endothelial cells','Macrophages/Microglia',
               'Oligodendrocytes','MES1-like', 'MES2-like', 'NPC1-like', 'NPC2-like', 'AC-like',
               'OPC-like', 'Glioblastoma Stem Cells') 

gcf$final_annot <- factor(current_idents, levels = new.order)
Idents(gcf) <- 'final_annot'
DimPlot(gcf, label=F) 

#-------------------------SHORT_ANNOT--------------------
current_idents <- as.character(Idents(gcf))
new_idents <- current_idents

new_idents[new_idents == "Oligodendrocytes"] <- "OLs"
new_idents[new_idents == "Endothelial cells"] <- "ECs"
new_idents[new_idents == "Glioblastoma Stem Cells"] <- "GSCs"
new_idents[new_idents == "Macrophages/Microglia"] <- "M/M"

new.order <- c('CAFs', 'ECs', 'M/M', 'OLs',
               'MES1-like','MES2-like','NPC1-like','NPC2-like', 'AC-like','OPC-like','GSCs')
ordered_idents <- factor(new_idents, levels = new.order)
gcf$short_annotation <- ordered_idents
Idents(gcf) <- "short_annotation"
DimPlot(gcf)
#--------------------------PLOTS----------------------------------------------------------------
 
colors_initial <- c(
  "CAFs" = 'red',                
  "T lymphocytes" = 'forestgreen',
  "Macrophages/Microglia" = 'purple3',
  "Oligodendrocytes" = 'dodgerblue',
  "GBM" = 'gold1')
colors_final <- c(
  "CAFs" = 'red',                
  "Endothelial cells" = 'blue',
  "T lymphocytes" = 'forestgreen',
  "Macrophages/Microglia" = 'purple3',
  "Oligodendrocytes" = 'dodgerblue',
  "MES1-like" ='chocolate1',         
  "MES2-like" = 'chocolate4',          
  "NPC1-like" = 'cyan2',               
  "NPC2-like" = 'slateblue1',
  "AC-like" = 'black',
  "OPC-like" = 'orchid1',
  "Glioblastoma Stem Cells" = 'gold1' 
)

#Idents(gcf) <- 'final_annot'

svg('vln.svg', width = 8, height = 7)
DimPlot(gcf, reduction = "umap", group.by = 'final_annot', label = F, cols=colors_final)
dev.off()

Idents(gcf) <- 'short_annotation'
table(Idents(gcf))[['GSCs']]/sum(table(Idents(gcf)))
stromal_n <- 1138+812+4099+757
sum(table(Idents(gcf)))-stromal_n

cols <- c( 
  "CAFs" = 'red',       
  "Endothelial cells" = 'blue', 
  "Macrophages/Microglia" = 'purple3',
  "Oligodendrocytes" = 'deepskyblue',
  "GBM" = 'gold1' 
)
svg('vln.svg', width = 8, height = 7)
DimPlot(gcf, reduction = "umap", label = F, cols=cols)
dev.off()
#ACTA2, VIM, PDGFRB, COL1A1, COL3A1, COL4A2, 
#NR2F2, BGN, FAP, THY1, RGS5, VIM, FN1
feature <- 'FN1'
path <- paste0('CAF_markers/', feature, '.jpg')
jpeg(path, width = 1850, height = 1250, res=200)
VlnPlot(gcf, features=feature, cols=cols, log=TRUE, pt.size = 0,)
dev.off()



dp <- DotPlot(gcf, features = genes_validation_2, #group.by = 'init_annot',
              col.min = -2.5, col.max=2.5, scale.min = 0, scale.max = 100) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 1)) + # Поворот подписей на оси X
  ggtitle("") 


#svg("vln.svg", width = 14, height = 6) # for lot of  
svg("vln.svg", width = 8, height = 4) # cafs
dp + scale_color_gradientn(
  colours = c("black", "blue1", "deeppink1"),  # Четыре цвета
  values = c(0, 0.1, .55,  1),  # Точки привязки (нормализованные от 0 до 1)
  limits = c(-2.5, 2.5)  # Соответствует col.min и col.max в DotPlot
)
dev.off()

#--------------------CAF-MARKERS--------------------
caf.markers <- c('ACTA2', 'VIM', 'PDGFRB', 
                 'COL1A1', 'COL3A1', 'COL4A2', 
                 'NR2F2', 'BGN', 'FAP', 
                 'THY1', 'RGS5', 'FN1')
all <- c(caf.markers, c('VWF', 'PECAM1', 'ENG', 'NRP2', 'CD34',  'CDH5'))
dp <- DotPlot(gcf, features = all, #group.by = 'init_annot',
              col.min = -2.5, col.max=2.5, scale.min = 0, scale.max = 100) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 1)) + # Поворот подписей на оси X
  ggtitle("") 
svg("vln.svg", width = 8, height = 4) # cafs
dp + scale_color_gradientn(
  colours = c("black", "blue1", "deeppink1"),  # Четыре цвета
  values = c(0, 0.1, .55,  1),  # Точки привязки (нормализованные от 0 до 1)
  limits = c(-2.5, 2.5)  # Соответствует col.min и col.max в DotPlot
)
dev.off()

#--------------------------CNV----------------------------------------------------------------
 
counts <- GetAssayData(gcf, assay = "RNA", layer = "counts")
chromosome_means <- get_chromosome_means(counts)

write.csv(
  data.frame(chromosome_means),
  file = "chromosome_means_validation_v3.csv",
  row.names = TRUE,            # не пишем номер строки
  quote = FALSE                 # без кавычек, если не нужны
)



####### ML RES:
res <- read_csv('result_ML_val3.csv')
tail(res)
res_df <- as.data.frame(res)
head(res_df)

ann <- ifelse(res_df$pred == 1, "Glioblastoma cells", "Stromal cells")
names(ann) <- res_df$label
gcf$ML_annotation <- ann[Cells(gcf)]

svg('vln.svg', width = 8, height = 7)
DimPlot(gcf, group.by = 'ML_annotation', cols = c('#FF4500', 'deepskyblue'))
dev.off()

DimPlot(gcf, group.by = 'init_annot', cols = cols)

#-------CELL_CHAT--------------------------------------------------------------
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
cellchat <- cellchat_workflow(gcf)

print(cellchat@netP$pathways) # все пары взаимодействий

VlnPlot(gcf, 'SPP1', pt.size = 0, log=T)

cafs_pairs <- data.frame(interaction_name = c('PTN_PTPRZ1',
                                              'PPIA_BSG', 
                                              'PGF_VEGFR1',
                                              'MDK_PTPRZ1',
                                              'MDK_LRP1',
                                              'MDK_ITGA6_ITGB1', 
                                              'JAG1_NOTCH4.1', 
                                              'IGFBP7_CD93',
                                              'FGF1_FGFR1',
                                              'ANXA1_FPR1',
                                              'ANGPT2_ITGA5_ITGB1'
                                              ))
ec_pairs <- data.frame(interaction_name = c('PGF_VEGFR1',#
                                            'TGFB1_TGFBR2_ACVRL1',
                                            'PDGFB_PDGFRB', 
                                            'MDK_ITGA4_ITGB1',
                                            'MDK_ITGA6_ITGB1', 
                                            'JAG1_NOTCH4.1',
                                            'IGFBP7_CD93',
                                            'ANXA1_FPR1',
                                            'GRN_SORT1',
                                            'HBEGF_EGFR',
                                            'ANGPT2_ITGA5_ITGB1',
                                            'ADM_CALCRL',
                                            'MDK_PTPRZ1',
                                            'PTN_PTPRZ1'))

svg("vln.svg", width = 4, height = 2.6) 
netVisual_bubble(cellchat, 
                 sources.use = c('ECs'), 
                 pairLR.use = ec_pairs,
                 remove.isolate = T,
                 font.size=10,
                 n.colors = 11,
                 color.grid = "black",  
                 angle.x = 90) +  
  scale_color_gradientn(
    colors = c("black", 'cyan', "blue", "magenta"),
    limits = c(0, 1))
dev.off() 

feature_plot(val, 'PECAM1')
 

extractEnrichedLR(cellchat, signaling = 'PTN', geneLR.return = FALSE)



#All communications with related information:
comm_df <- subsetCommunication(cellchat)   
comm_df_pathway <- subsetCommunication(cellchat, slot.name = "netP")
write.xlsx(comm_df, file = "CellChat_Communications_validation_3.xlsx", sheetName = "Ligand_Receptor", rowNames = FALSE)
 