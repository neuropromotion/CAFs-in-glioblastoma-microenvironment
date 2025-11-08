library(CellChat)
library(openxlsx) 
make_short_annotation <- function(obj, new.order){
  current_idents <- as.character(Idents(obj))
  new_idents <- current_idents
  
  new_idents[new_idents == "Oligodendrocytes"] <- "OLs"
  new_idents[new_idents == "Endothelial cells"] <- "ECs"
  new_idents[new_idents == "Glioblastoma Stem Cells"] <- "GSCs"
  new_idents[new_idents == "Macrophages/Microglia"] <- "M/M"
  new_idents[new_idents == "T lymphocytes"] <- "T cells"
   
  ordered_idents <- factor(new_idents, levels = new.order) 
  obj$short_annotation <- ordered_idents
  Idents(obj) <- "short_annotation"
  return(obj)
} # optinal (for ex: Endothelial cells = ECs)
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
  interaction_input <- read.csv(file = 'cellchat_configue/interaction_input_CellChatDB.csv')
  complex_input <- read.csv(file = 'cellchat_configue/complex_input_CellChatDB.csv', row.names = 1)
  cofactor_input <- read.csv(file = 'cellchat_configue/cofactor_input_CellChatDB.csv', row.names = 1)
  geneInfo <- read.csv(file = 'cellchat_configue/geneInfo_CellChatDB.csv', row.names = 1)
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

cellchat <- cellchat_workflow(merged_filtered)

print(cellchat@netP$pathways) # all pairs

# contribution pathways in whole network 
netAnalysis_contribution(cellchat, c(cellchat@netP$pathways)[1:20], 
                         title = 'Contribution of each LR pair')
# SCATTER PLOT
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")  
netAnalysis_signalingRole_scatter(cellchat, signaling = cellchat@netP$pathways)

# BUBBLE PLOT 

pairLR_df <- data.frame(interaction_name = c("PTN_NCL", 'MDK_NCL'))
pairLR_df <- data.frame(interaction_name = c("PTN_PTPRZ1", 'MDK_PTPRZ1'))

svg("vln.svg", width = 3.5, height = 3.5) 
netVisual_bubble(cellchat,
                 #signaling = 'GRN',
                 sources.use = c('CAFs', 'ECs', 'OLs', 'M/M'),
                 #targets.use = c('ECs'),
                 targets.use = c('T cells', 'M/M'),
                 #sources.use = c('ECs', 'M/M', 'OLs', 'T cells'),
                 sort.by.source= T,
                 #pairLR.use = pairLR_df, 
                 remove.isolate = T,
                 font.size=10,
                 #n.colors = 20,
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
write.xlsx(comm_df, file = "CellChat_Communications.xlsx", sheetName = "Ligand_Receptor", rowNames = FALSE)
#### #### #### #### #### #### #### #### #### 
