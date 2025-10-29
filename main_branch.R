library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(gridExtra)
library(writexl)
library(ggplot2) 
library(devtools)
library(CellChat)
library(DoubletFinder)
library(GSVA)
library(readxl) 
library(harmony)
library(celda)
library(SingleCellExperiment)
library(scran)
# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135045
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
# Violin plots of features, counts and mito RNA
count_plots <- function(obj, 
                        feature_lower=600,
                        feature_upper=6000,
                        count_lower=200,
                        count_upper=1e5,
                        mito=15,
                        pt.size=0){
  
  vln_plot1 <- VlnPlot(obj, features = "nFeature_RNA", layer = "counts", pt.size=pt.size) + 
    theme(panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),       # Убрать подпись оси x
          axis.text.x = element_blank(),        # Убрать метки оси x
          axis.ticks.x = element_blank()) + 
    geom_hline(yintercept = feature_upper, linetype = "dashed", color = "red", cex=1)+ 
    geom_hline(yintercept = feature_lower, linetype = "dashed", color = "red", cex=1)+ 
    NoLegend()
  
  vln_plot2 <- VlnPlot(obj, features = "nCount_RNA", layer = "counts", pt.size=pt.size) + 
    theme(panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),       # Убрать подпись оси x
          axis.text.x = element_blank(),        # Убрать метки оси x
          axis.ticks.x = element_blank()) + 
    geom_hline(yintercept = count_upper, linetype = "dashed", color = "red", cex=1)+ 
    geom_hline(yintercept = count_lower, linetype = "dashed", color = "red", cex=1)+ 
    NoLegend()
  
  vln_plot3 <- VlnPlot(obj, features = "mitoPercent", layer = "counts", pt.size=pt.size) + 
    theme(panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),       # Убрать подпись оси x
          axis.text.x = element_blank(),        # Убрать метки оси x
          axis.ticks.x = element_blank()) + 
    geom_hline(yintercept = mito, linetype = "dashed", color = "red", cex=1) +
    NoLegend()
  combined_plot <- vln_plot1 | vln_plot2 | vln_plot3
  return(combined_plot)
}
decontex_workflow <- function(seu){
  counts_mat <- GetAssayData(seu, layer = "counts", assay = "RNA") # raw counts matrix extraction
  cell_md <- seu@meta.data # metadata extraction

  sce <- SingleCellExperiment(
    assays = list(counts = counts_mat),
    colData = cell_md
  ) 
  set.seed(42) # fix seed to reproducing 
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
  seu <- add_groups(seu)
  sce <- decontX(sce, z = seu$soup_group) # run decontX
  # Extraction clear count and metadata
  clean_counts <- decontXcounts(sce)
  clean_md     <- as.data.frame(colData(sce))
  # create new clean Seurat-object
  seu_clean <- CreateSeuratObject(
    counts    = clean_counts,
    project   = "GBM_decontaminated",
    meta.data = clean_md
  )
  return(seu_clean)
} # DecontX pipeline 
stantard_workflow <- function(seu, n_components=15, res=0.1){
  seu <- JoinLayers(seu)
  seu <- decontex_workflow(seu)
  seu <- NormalizeData(seu)                             
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)   
  seu <- ScaleData(seu, features = rownames(seu)) 
  seu <- RunPCA(seu, features = VariableFeatures(object = seu))
  seu <- RunHarmony(
    object = seu,
    group.by.vars = "Sample"
  )
  seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:n_components)
  seu <- FindClusters(seu, resolution = res)
  seu <- RunUMAP(seu, reduction = "harmony", dims = 1:n_components)
  seu <- deduplex(seu)
  return(seu)
} # Norm, Scale, PCA, UMAP

# save and load preprocessed data
#saveRDS(merged_filtered, file = "/mnt/jack-5/amismailov/CAF_study/merged_filtered_v2.rds")
#merged_filtered <- readRDS(file = "/mnt/jack-5/amismailov/CAF_study/merged_filtered_v2.rds")

#--------------------------------------LOAD DATA-----------------------------------------
path = 'path_to_dir'
dirs <- c("sample1", "sample2","sample3","sample4","sample5","sample6","sample7")
for (x in dirs) {
  file_path <- paste0(path, '/', x, '/data.txt')
  expression_matrix <- read.delim(file_path, row.names = 1)
  assign(x, CreateSeuratObject(counts = expression_matrix))
}
merged <- merge(sample1, y = c(sample2, sample3, sample4, sample5, sample6, sample7),
                add.cell.id  = dirs,
                project = 'GBM')

merged$sample <- rownames(merged@meta.data)
#now we separate previous share colums into sample name and barcode 
merged@meta.data <- separate(merged@meta.data, col = 'sample', into = c('Sample', 'Barcode'),
                             sep = '_')
#check 
unique(merged$Sample)

#-----------------------------------------QC-and-filtering-----------------------------------------------
# adding quality metrics
merged$mitoPercent <- PercentageFeatureSet(merged, pattern = '^MT-')
# plot
count_plots(merged, feature_lower = 600, feature_upper = 4500, 
            count_lower = 200, count_upper = 10e3, mito=20) # plot nFeatures, nCounts, MT percent
# filter cells
merged_filtered <- subset(merged, subset = nFeature_RNA > 600 & nCount_RNA <= 10e3 &
                            nFeature_RNA < 4500 & mitoPercent < 20)

#print differences in sizes
merged_filtered
merged

#-----------------------Main preprocess pipeline-----------------------------------------------
# DeconteX, Normalization, Scaling, PCA, Harmony, UMAP, clustering, DoubledFinder filtering 
merged_filtered <- stantard_workflow(merged_filtered, n_components = 20, res=0.1) 
#---------------------------------Basic plots-----------------------------------------------
DimPlot(merged_filtered, reduction = 'umap', group.by = 'Sample')
DimPlot(merged_filtered, reduction = 'umap', label = TRUE)

feature_plot(merged_filtered, 'HBEGF')
#---------------------------------Differentially expressed genes-----------------------------------------------
FAM <- FindAllMarkers(merged_filtered,
                      logfc.threshold = 0.25, 
                      min.pct = 0.1,
                      only.pos = TRUE,
                      test.use = 'wilcox',
                      layer = 'data')
FAM <- read_xlsx("path_to_dir/DEG_main_final.xlsx")
write_xlsx(FAM, "path_to_dir/DEG_main_final.xlsx")

#---------------------------------Rename idents-----------------------------------------------
#merged_filtered <- FindClusters(merged_filtered, resolution = 0.1)

new_cluster_names <- c(
  '0' = 'Macrophages/Microglia',
  '1' = 'GBM',
  '2' = 'GBM',
  '3' = 'GBM',
  '4' = 'Oligodendrocytes',
  '5' = 'T lymphocytes',
  '6' = 'Macrophages/Microglia',
  '7' = 'CAFs',
  '8' = 'Endothelial cells',
  '9' = 'Undefined'
)
merged_filtered <- RenameIdents(merged_filtered, new_cluster_names)
merged_filtered$init_annot <- as.character(Idents(merged_filtered)) # initial annotation
merged_filtered <- subset(merged_filtered, idents = "Undefined", invert = TRUE) # DELETE INDEFINED CLUSTER

#---------------------SET YOUR ORDER-------------------------------
# initial annotation:
new.order <- c('CAFs', 'Endothelial cells', 'T lymphocytes', 'Macrophages/Microglia', 'Oligodendrocytes', 'GBM') 
current_idents <- as.character(Idents(merged_filtered)) # get current idents
ordered_idents <- factor(current_idents, levels = new.order) # set your order
merged_filtered$init_annot <- ordered_idents # add new annotation
Idents(merged_filtered) <- "init_annot" # make in main
# final annotation:
# To add glioblastoma subtypes use script GSVA.R -> final annotation
# Set order after GSVA annotation (GBM subclustering)
new.order <- c('CAFs', 'Endothelial cells', 'T lymphocytes', 'Macrophages/Microglia', 'Oligodendrocytes',
               'MES1-like','MES2-like','NPC1-like','NPC2-like', 'AC-like','OPC-like','Glioblastoma Stem Cells') 
current_idents <- as.character(Idents(merged_filtered)) # get current idents
ordered_idents <- factor(current_idents, levels = new.order)
merged_filtered$final_annot <- ordered_idents
Idents(merged_filtered) <- "final_annot"

#---------------------------------MAIN PLOTS-----------------------------------------------
colors_samples <- c('red', 'blue', 'gold', 'purple3', 'dodgerblue', 'forestgreen', 'orchid1')
colors_initial <- c('red','blue','forestgreen','purple3', 'dodgerblue','gold')
colors_final <- c(
  'red','blue','forestgreen','purple3','dodgerblue',
  'chocolate1','chocolate4',
  'cyan2','slateblue1','black','orchid1','gold' 
)

p <- DimPlot(merged_filtered, reduction = 'umap', group.by = 'Sample', cols = colors_samples)
p <- DimPlot(merged_filtered, reduction = 'umap', group.by = 'init_annot', cols = colors_initial)
# final annotation after GSVA
p <- DimPlot(merged_filtered, reduction = 'umap', group.by = 'final_annot', cols = colors_final) # After GSVA glioblastoma cells will be subclusteres

# save SVG
svg('vln.svg', width = 8, height = 7)
p
dev.off()  






# GET TOP 5 DEG for each cluster
top.genes <- character() # top DEG genes for each cluster (unique)
for (cl in new.order) {
  cluster_genes <- FAM$gene[FAM$cluster == cl]    
  added <- 0                                      
  idx   <- 1                                      
  
  # add top 5 unique genes
  while (added < 5 && idx <= length(cluster_genes)) {
    g <- cluster_genes[idx]
    if (!(g %in% top.genes)) {                    
      top.genes <- c(top.genes, g)
      added <- added + 1
    }
    idx <- idx + 1                             
  }
}
# VIOLIN PLOTS (CAF MARKERS)
caf.markers <- c('ACTA2', 'VIM', 'PDGFRB', 
                 'COL1A1', 'COL3A1', 'COL4A2', 
                 'NR2F2', 'BGN', 'FAP', 
                 'THY1', 'RGS5', 'FN1')
# dotplot
dp <- DotPlot(merged_filtered, features = top.genes, group.by='short_annotation',
              col.min = -2.5, col.max=2.5, scale.min = 0, scale.max = 100) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(angle = 0, hjust = 1)) + # Поворот подписей на оси X
  ggtitle("") 

svg("vln.svg", width = 11, height = 5) # for DEG 
dp + scale_color_gradientn(
  colours = c("black", "blue1", "deeppink1"),  # Четыре цвета
  values = c(0, 0.1, .55,  1),  # Точки привязки (нормализованные от 0 до 1)
  limits = c(-2.5, 2.5)  # Соответствует col.min и col.max в DotPlot
)
dev.off()

colors_final <- c(
  'red','blue','forestgreen','purple3','dodgerblue',
  'chocolate1','chocolate4',
  'cyan2','slateblue1','black','orchid1','gold' 
)
### PLOT certain clusters:
target_cluster <- "Glioblastoma Stem Cells" 
highlight_color <- "gold" 

svg('vln.svg', width = 8, height = 7)
DimPlot(merged_filtered, reduction = 'umap', group.by = 'final_annot', pt.size = 0.01,
        cells.highlight = WhichCells(merged_filtered, idents = target_cluster),
        cols.highlight = highlight_color,  
        sizes.highlight = 0.01,
        cols = "gray")  + 
  scale_color_manual(labels = c("Other Cells", target_cluster),
                     values = c("gray", highlight_color)) +
  labs(color = "Cell Type") +
  ggtitle(target_cluster)
dev.off() 

#----------plots
feature_plot <- function(data, feature){
  FeaturePlot(data, feature, order=F, slot = 'data')  +
    scale_color_gradientn(
      colours = c("gray", "blue1", "deeppink1"),  
      #values = c(0, 0.3, .65,  1),  
    )
}
feature_plot(merged_filtered, 'HBEGF')
feature_plot(merged_filtered_old, 'HBEGF')

FeaturePlot(merged_filtered, 'SOX2')
FeaturePlot(merged_filtered_old, 'SOX2')

feat <- 'PDGFRA'
VlnPlot(merged_filtered, features = feat, 
        cols  = colors_final, group.by = 'final_annot',
        log = T, pt.size = 0)

VlnPlot(merged_filtered_old, features = feat,
        cols  = colors_final, group.by = 'final_annot',
        log = T, pt.size = 0)



#----------------------ML ANNOTATION------------------------------
res <- read_csv('path_to_dir/result_ML_main_v2.csv')
tail(res)
res_df <- as.data.frame(res)
head(res_df)

ann <- ifelse(res_df$pred == 1, "Glioblastoma cells", "Stromal cells")
names(ann) <- res_df$label
merged_filtered$ML_annotation <- ann[Cells(merged_filtered)]

svg('vln.svg', width = 8, height = 7)
DimPlot(merged_filtered, group.by = 'ML_annotation', cols = c('#FF3E96', '#1E90FF'))
dev.off()
