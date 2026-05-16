source('functions.R')
library(harmony)
library(Seurat)
library(Matrix)
# LOAD [Neftel et al. 2019] GSM3828672
main_dataset_part_1 <- function(path){
  genes <- read.delim(paste0(path, 'Genes.txt'), header = FALSE)
  cells <- read.csv(paste0(path, 'Cells.csv'), header = TRUE)
  tpm <- readMM(paste0(path, 'Exp_data_TPM.mtx'))
  
  rownames(tpm) <- genes$V1
  colnames(tpm) <- paste0('GSM3828672_', cells$cell_name)
  
  obj <- CreateSeuratObject(counts = tpm, project = 'GSM3828672')
  obj@meta.data$orig.ident <- paste0('GSM3828672_', cells$sample)
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")   
  
  obj <- subset(obj, subset = nFeature_RNA > 300 & 
                        nFeature_RNA < 10e3)
  
  obj <- deduplex_v3(obj, join=T)
  obj <- decontex_workflow(obj)
  gc()
  return(obj)
}

# LOAD 10x
main_dataset_part_2 <- function(path){
  counts <- Read10X(data.dir = path)
  
  x10 <- CreateSeuratObject(
    counts = counts,
    project = "10x",
    min.cells = 3,
    min.features = 0
  )
  x10[["percent.mt"]] <- PercentageFeatureSet(x10, pattern = "^MT-")  
  x10 <- subset(x10, subset = nFeature_RNA > 300 & 
                  nCount_RNA < 80e3 &
                  percent.mt < 20)
  x10 <- deduplex_v3(x10)
  x10 <- decontex_workflow(x10)
  gc()
  return(x10)
}



# LOAD SMARTSEQ2 (GSE135045)
main_dataset_part_3 <- function(path, project_name='GSE135045'){
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
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")  
  obj <- subset(obj, subset = percent.mt < 20 & 
                  nFeature_RNA > 300 &  
                  nCount_RNA < 25e3)
  
  obj <- deduplex_v3(obj, join = T)
  obj <- decontex_workflow(obj)
  gc()
  
  return(obj)
}



# LOAD weizmann (GSE103224)
main_dataset_part_4 <- function(path){
  genes <- read.delim(paste0(path, 'Genes.txt'), header = FALSE)
  cells <- read.csv(paste0(path, 'Cells.csv'), header = TRUE)
  umi <- readMM(paste0(path, 'Exp_data_UMIcounts.mtx'))
  
  rownames(umi) <- genes$V1
  colnames(umi) <- paste0('GSE103224_', cells$cell_name)
  
  obj <- CreateSeuratObject(
    counts = umi,
    project = 'GSE103224',
    min.cells = 3,
    min.features = 0
  )
  
  obj@meta.data$orig.ident <- paste0('GSE103224', cells$sample)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = nFeature_RNA > 300 & 
                        nFeature_RNA < 6e3 &
                        nCount_RNA < 15e3
  )
  obj <- deduplex_v3(obj, join=T)
  obj <- decontex_workflow(obj)
  return(obj)
}

main_dataset_merge <- function(){
  main_1 <- main_dataset_part_1(path='/')
  main_2 <- main_dataset_part_2(path='/')
  main_3 <- main_dataset_part_3(path='/')
  main_4 <- main_dataset_part_4(path='/')
  
  combined_all <- merge(main_1, y = c(main_2, main_3, main_4), add.cell.ids = c("GSM3828672", "10X", "GSE135045", "GSE103224"))
  obj_list <- SplitObject(combined_all, split.by = "orig.ident") 
  
  obj_list <- lapply(X = obj_list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    return(x)
  })
  
  features <- SelectIntegrationFeatures(object.list = obj_list)
  
  combined_final <- merge(obj_list[[1]], y = obj_list[2:length(obj_list)])
  
  VariableFeatures(combined_final) <- features
  combined_final <- ScaleData(combined_final, vars.to.regress = "percent.mt")
  combined_final <- RunPCA(combined_final, npcs = 30)
  
  
  combined_final <- RunHarmony(combined_final, group.by.vars = "orig.ident")
  
  
  combined_final <- RunUMAP(combined_final, reduction = 'harmony', dims=1:10)
  combined_final <- FindNeighbors(combined_final, reduction = 'harmony', dims=1:10)
  combined_final <- FindClusters(combined_final, res=0.1, algorithm = 1)
  
  gc()
  return(combined_final)
}


# LOAD [GSE173278] 10x
validation_dataset_load <- function(){
  path <- '/mnt/jack-5/amismailov/CAF_study/GCF/filtered'
  
  counts <- Read10X(data.dir = path)
  seu <- CreateSeuratObject(counts = counts, project = "filtered", min.cells = 0, min.features = 0)
  rm(gcf)
  gc()
  
  cts_seu <- seu[["RNA"]]$counts
  seu$nFeature_raw <- Matrix::colSums(cts_seu > 0)   # Кол-во детектированных генов
  seu$nCount_raw   <- Matrix::colSums(cts_seu)       # Суммарные UMI
  seu$percent.mt   <- PercentageFeatureSet(seu, pattern = "^MT-")  # Митохондрии
  
  metadata <- read.csv(gzfile(paste0(path, "/metadata.csv.gz")),
                       row.names = 1)
  
  seu <- AddMetaData(seu, metadata = metadata)
  seu <- deduplex_v3(seu)
  seu <- decontex_workflow(seu)
  
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)   
  seu <- ScaleData(seu, features = rownames(obj_simple)) 
  seu <- RunPCA(seu, features = VariableFeatures(object = obj_simple))
  
  seu <- RunHarmony(seu, group.by.vars = c('patient'))
  
  

  seu <- FindNeighbors(seu, reduction='pca', dims = 1:5)
  seu <- FindClusters(seu, resolution = 0.1)
  seu <- RunUMAP(seu, reduction='pca', dims = 1:7)
  
  rm(cts_seu, metadata)
  gc()
  return(seu)
}












