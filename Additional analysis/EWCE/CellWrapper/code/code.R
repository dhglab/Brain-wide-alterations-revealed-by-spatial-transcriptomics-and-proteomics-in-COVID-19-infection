options(stringsAsFactors = F)

library(EWCE)
library(data.table)
library(tidyverse)
library(Seurat)

generate.celltype.wrapper <- function(expr_matrix, cell_annotation , thresh=0.01, data_set_name = "", output_dir){
  previous_wd <- getwd()
  on.exit(setwd(previous_wd))
  level2_types = unique(cell_annotation$level2class)
  
  remove_gene= rep(NA,nrow(expr_matrix))
  names(remove_gene) = rownames(expr_matrix)
  ind=0
  remove_genes_mat <- matrix(NA, nrow = nrow(expr_matrix), ncol = length(level2_types) )
  rownames(remove_genes_mat) <- row.names(expr_matrix)
  colnames(remove_genes_mat) <- level2_types
  
  for(cell in level2_types){
    cell_expr <- expr_matrix[,which(cell_annotation$level2class==cell)]
    remove_genes_mat[,cell] <- rowMeans(as.matrix(cell_expr))
  }
  
  remove_gene <- apply(remove_genes_mat>thresh, 1, sum) == 0
  expr_matrix = expr_matrix[!remove_gene,]
  
  print(paste(nrow(expr_matrix),"genes remaining for EWCE specificity calculation",sep=" "))
  
  ### Equal Subtypes
  print("Begin Equal Subtypes Method: now getting means across genes")
  
  level2_types = unique(cell_annotation$level2class)
  expr_average = matrix(NA,nrow=nrow(expr_matrix),ncol=length(level2_types))
  rownames(expr_average) = rownames(expr_matrix)
  colnames(expr_average) = level2_types
  
  for(cell in level2_types){
    expr_average[,cell] = rowMeans(expr_matrix[,which(cell_annotation$level2class==cell)])
  }
  expr_annot_average = data.frame("cell_id"=level2_types,
                                  "level2class"=level2_types,
                                  "level1class"=cell_annotation$level1class[match(level2_types,cell_annotation$level2class)])
  expr_matrix1 = expr_average
  expr_annot1 = expr_annot_average
  
  annotLevels = list(level1class = expr_annot1$level1class, level2class = expr_annot1$level2class)
  print("generating cell type data")
  current_output <- file.path(output_dir, "Equal_Level2_CellType")
  dir.create(current_output,showWarnings = F)
  setwd(current_output)
  generate_celltype_data(exp = expr_matrix1, annotLevels = annotLevels,groupName = data_set_name,savePath = current_output)
  
  
  ### Size Biased Subtypes
  print("Begin Size Biased Subtypes Method")
  annotLevels = list(level1class = cell_annotation$level1class, level2class = cell_annotation$level2class)
  print("generating cell type data")
  current_output <- file.path(output_dir, "SizeBiased_Level2_CellType")
  dir.create(current_output,showWarnings = F)
  setwd(current_output)
  generate_celltype_data(exp = expr_matrix,annotLevels = annotLevels,groupName = data_set_name,savePath = current_output)
}

setwd("./Azimuth_Bakken")
#import Azimuth data (downloaded from Azimuth website, matrix_Bakken, metadata_Bakken)
Bakkenctrl <- vroom::vroom("/u/project/geschwind/tzhang/tzhang/Dementia/EWCE/AdditionalCelltypeEWCE_SC/Azimuth_Bakken/matrix_Bakken.csv")
save(Bakkenctrl, file="Bakkenctrl.rda")

head(Bakkenctrl)[1:10]

Bakkenctrl <- Bakkenctrl[order(Bakkenctrl$sample_name),]
metadata_Bakken <- metadata_Bakken[order(metadata_Bakken$sample_name),]

cell_annotation = data.frame("cell_id"= row.names(Bakkenctrl), 
                             "level1class"= metadata_Bakken$class_label,
                             "level2class" = metadata_Bakken$subclass_label)

Bakkenctrl <- as.data.frame(Bakkenctrl)
rownames(Bakkenctrl) <- Bakkenctrl$sample_name
head(Bakkenctrl)[1:10]

Bakkenctrl <- Bakkenctrl[, -1]
head(Bakkenctrl)[1:10]

Bakkenctrl_t <- as.data.frame(t(Bakkenctrl))

generate.celltype.wrapper(expr_matrix = Bakkenctrl_t,
                          cell_annotation = cell_annotation,
                          data_set_name = "Bakken_anno_ctrl",
                          output_dir = "/u/project/geschwind/tzhang/tzhang/Dementia/EWCE/AdditionalCelltypeEWCE_SC/Azimuth_Bakken/")
print("Done")


setwd("./Yang_COVID")
#import Yang_et_al seurat_obj
seurat_obj <- `COVID-19_brain_snRNA-seq_parenchyma_cortex_final_seurat_v3.2.3`
Idents(object = seurat_obj) <- "cellID2"
ctrl_cells <- seurat_obj[,seurat_obj$Biogroup =="Control"]
saveRDS(ctrl_cells, file="Yang_COVID_scdata_ctrl.rds")

ctrl_cells <- Yang_COVID_scdata_ctrl
DefaultAssay(object = ctrl_cells) <- "RNA"
sc_data <- GetAssayData(object =ctrl_cells , slot = "data")
colnames(sc_data) <- Idents(ctrl_cells)
dim(sc_data)
#[1] 30861 15946

sc_data_mat <- as.matrix(sc_data)

unique(colnames(sc_data_mat))

ctl_anno = data.frame("cell_id"=colnames(sc_data_mat), 
                      "level1class"=ctrl_cells$cellID,
                      "level2class" =ctrl_cells$cellID2)

generate.celltype.wrapper(expr_matrix = sc_data_mat,
                          cell_annotation = ctl_anno,
                          data_set_name = "Yang_anno_ctl",
                          output_dir = "/u/project/geschwind/tzhang/tzhang/Dementia/EWCE/AdditionalCelltypeEWCE_SC/Yang_COVID/")
print("Done")

setwd("./WamsleyASD") #contributed by Lucy Bricks
library(Seurat)
ASDobject <- readRDS("/u/project/geschwind/bwamsley/Seurat/ASDobject.rds")
ctrl_cells <- ASDobject[,ASDobject$Diagnosis == "CTL"]

sc_data <- GetAssayData(object = ctrl_cells, slot = "data")
colnames(sc_data) <- Idents(ctrl_cells)
dim(sc_data)

idx <- floor(runif(50000, min=1, max=284696))
idx <- unique(idx)
sc_data <- sc_data[,idx]

sc_data_mat <- as.matrix(sc_data)

unique(colnames(sc_data_mat))

cell_class <- data.frame(level2class = unique(colnames(sc_data_mat))) %>% 
  mutate(level1class = case_when(
    grepl("ODC_1|ODC_2|ODC_3", level2class) ~ "ODC",
    grepl("OPCS_1|OPCS_2|OPCS_3|OPCS_4", level2class) ~ "OPC",
    grepl("EXT_2_L56|EXT_4_L56|EXT_12_L56|EXT_13_L45|EXT_9_L6|EXT_5_L56|EXT_8_L6|EXT_14_L5B", level2class) ~ "ExNeu_deep",
    grepl("EXT_3_L23|EXT_6_L23|EXT_7_L23|EXT_10_L34|EXT_11_L34|EXT_1_L4", level2class) ~ "ExNeu_superficial",
    grepl("INT_1_KIT|INT_2_PV_SULF1|INT_3_RELN_NDNF|INT_4_VIP_CR|INT_5_SST|INT_6_PV", level2class) ~ "InNeu",
    grepl("MG1|MG2", level2class) ~ "MG",
    grepl("BBB_Pericytes|BBB_EndoMuralMicro|BBB_EndoMural|BBB_Endo_1", level2class) ~ "BBB",
    grepl("ASTRO_1|ASTRO_2|ASTRO_3", level2class) ~ "Astro",
  ) )

table(cell_class)
ASDobect_anno = data.frame("cell_id"=rownames(metadata_Bakken), 
                           "level2class" = colnames(sc_data_mat)) %>% 
  left_join(cell_class , by = "level2class")

generate.celltype.wrapper(expr_matrix = sc_data_mat,
                          cell_annotation = ASDobect_anno,
                          data_set_name = "ASDobect_anno_ctl",
                          output_dir = "/u/project/geschwind/tzhang/tzhang/Dementia/EWCE/AdditionalCelltypeEWCE_SC/WamsleyASD/")
print("Done")