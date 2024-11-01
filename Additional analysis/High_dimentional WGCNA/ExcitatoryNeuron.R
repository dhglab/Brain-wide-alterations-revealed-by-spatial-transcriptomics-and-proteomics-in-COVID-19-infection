
install.packages("BiocManager")
BiocManager::install()
install.packages(c("Seurat", "WGCNA", "igraph", "devtools"))

devtools::install_github('smorabit/hdWGCNA', ref='dev')
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggplot2)
library(sctransform)
library(WGCNA)
library(hdWGCNA)
library(dplyr)

set.seed(12345)
enableWGCNAThreads()

setwd("./Fig3SourceData/hdWGCNA")

#import Yang, et al seurat_obj
seurat_obj <- readRDS("/u/project/geschwind/tzhang/tzhang/Dementia/EWCE/AdditionalCelltypeEWCE_SC/Yang_COVID/COVID-19_brain_snRNA-seq_parenchyma_cortex_final_seurat_v3.2.3.rds")

#select non-viral control and COVID-19 patients only
Idents(object = seurat_obj) <- "Disease"
seurat_obj_COVID <- subset(x = seurat_obj, idents = c("Non-viral", "COVID-19"), invert = FALSE)
saveRDS(seurat_obj_COVID, "seurat_obj_COVID.rds")

seurat_obj <- seurat_obj_COVID

Idents(object = seurat_obj) <- "cellID"
DefaultAssay(seurat_obj) <- "RNA"

seurat_obj@meta.data$Age_s =scale(as.integer(seurat_obj@meta.data$Age))
seurat_obj@meta.data$Batch = as.factor(seurat_obj@meta.data$Batch)
seurat_obj@meta.data$Sex = as.factor(seurat_obj@meta.data$Sex)
seurat_obj@meta.data$percent.mt_s = scale(seurat_obj@meta.data$percent.mt)
seurat_obj@meta.data$nCount_RNA_s = scale(seurat_obj@meta.data$nCount_RNA)
seurat_obj@meta.data$nFeature_RNA_s = scale(seurat_obj@meta.data$nFeature_RNA)
seurat_obj@meta.data$nNuclei_s = scale(as.numeric(seurat_obj@meta.data$nNuclei))

seurat_obj <- ScaleData(
  seurat_obj,
  features = rownames(seurat_obj),
  vars.to.regress = c("Age_s","Sex","Batch","percent.mt_s","nCount_RNA_s","nFeature_RNA_s","nNuclei_s"),
  model.use = "linear",
  do.scale = TRUE,
  do.center = TRUE,
  block.size = 32000,
  verbose = TRUE
)

saveRDS(seurat_obj, "seurat_obj_scaleregress.rds")

#set assay with scaledata
seurat_obj <- seurat_obj_scaleregress
scaledata <- GetAssayData(object = seurat_obj, slot = "scale.data")
dim(scaledata)

scaledata <- as(scaledata, "dgCMatrix")
dim(scaledata)
seurat_obj <- SetAssayData(
  object = seurat_obj,
  assay = 'RNA',
  new.data = scaledata,
  slot = 'data'
)

seurat_obj <- SetupForWGCNA(seurat_obj, gene_select = "fraction",
                            fraction = 0.05,wgcna_name = "COVID_sc")

length(seurat_obj@misc$COVID_sc$wgcna_genes)

seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- RunUMAP(seurat_obj, dims = 1:8)

DimPlot(seurat_obj, reduction = "pca")
DimPlot(seurat_obj, reduction = "umap")

#Generate metacells by CellID and Sample
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cellID","Sample"),
  reduction = "pca",
  k = 25, 
  max_shared = 0, 
  ident.group = "cellID" 
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
  seurat_obj,
  assay = "RNA",
  slot = "data",
  group_name = "Excitatory neuron", 
  group.by='cellID'
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' 
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
pdf("thresh_hold_Excitatory.pdf")
wrap_plots(plot_list, ncol=2)
dev.off()

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=8,
  setDatExpr=FALSE,
  tom_name = 'Excitatory neuron_regress'
)

#Figure3A
pdf("Excitatory_Dendrogram_scaledata.pdf", height=4, width=6)
PlotDendrogram(seurat_obj, main='Excitatory neuron_Dendrogram')
dev.off()

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="Sample"
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cellID', group_name = 'Excitatory neuron',
  sparse=FALSE
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Excitatory-M"
)

# plot genes ranked by kME for each module
pdf("kME_Scale_Excitatory.pdf")
PlotKMEs(seurat_obj, ncol=5)
dev.off()

# get the module assignment table:
modules <- GetModules(seurat_obj)

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)

saveRDS(seurat_obj, file='hdWGCNA_Scaledata_Excitatory_neuron.rds')
############################################################################################
# network and visualization packages
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(igraph)
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# compute gene scoring for the top 25 hub genes by kME for each module
# with Seurat method
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)

# make a featureplot of hMEs for each module

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', 
  order=TRUE 
)

# stitch together with patchwork
#Figure A1
pdf("Celltype_dimention.pdf", height=4, width=5.5)
DimPlot(seurat_obj, reduction = "umap")
dev.off()

#Figure A1_1
Idents(object = seurat_obj) <- "cellID2"
pdf("Celltype_dimention_2.pdf", height=8, width=9)
DimPlot(seurat_obj, reduction = "umap")
dev.off()

Idents(object = seurat_obj) <- "cellID"

pdf("Excitatory_Modulefeature_scale.pdf")
wrap_plots(plot_list, ncol=4)
dev.off()

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', 
  order='shuffle', 
  ucell = FALSE 
)

# stitch together with patchwork
pdf("Excitatory_Modulefeature_hubscores.pdf")
wrap_plots(plot_list, ncol=2)
dev.off()

# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

#Plot module network
#M3 networkplot is used for Figure3G
ModuleNetworkPlot(seurat_obj)

# hubgene network
pdf("Hubgenenetworkplot.pdf", height=7, width=7)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, 
  n_neighbors=25, 
  min_dist=0.4 
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
pdf('ModulefeatureUMAP_unsupervised_scale.pdf')
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, 
    size=umap_df$kME*2 
  ) +
  umap_theme()
dev.off()

pdf("ModuleUMAPplot_Excitatory_with_hubgenes.pdf")
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.4, # proportion of edges to sample (20% here)
  label_hubs=1 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  )
dev.off()

#Figure3D
pdf("ModuleUMAPplot_Excitatory_with_no_hubgenes.pdf")
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.4, 
  label_hubs=0,
  keep_grey_edges=FALSE
)
dev.off()

# run supervised UMAP:
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10,
  n_neighbors=25,
  min_dist=0.4,
  supervised=TRUE,
  target_weight=0.5
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
pdf("SupervisedUMAPmoduleExcitatory_2.pdf")
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()
dev.off()

saveRDS(seurat_obj, file='hdWGCNA_Scaledata_Excitatory_neuron.rds')

#####################Differential module eigengene (DME) analysis#############
unique(seurat_obj@meta.data$Disease)
#[1] "Non-viral" "COVID-19" 

# get MEs from seurat object
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "ME"
)

MEs <- GetMEs(seurat_obj, harmonized=FALSE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

New <- seurat_obj@meta.data[1:46]
New <- cbind(New, MEs)
ME_names <- colnames(MEs)

#scale covariants to numeric form#
unique(New$Disease)

New$Disease <- fct_relevel(New$Disease, "Non-viral")

Expr <- as.matrix(MEs)

New$nNuclei <- scale(as.numeric(New$nNuclei))
New$nFeature_RNA <- scale(as.numeric(New$nFeature_RNA))
formula <- as.formula(paste0("Expr ~ Disease + Age_s + Sex + Batch + percent.mt_s + nCount_RNA_s + nNuclei + nFeature_RNA"))
lm <- lm(formula, data = New)
Coefficentlist <- coef(summary(lm))
Coefficentlist[[1]]

lm_results <- lapply(ME_names, 
                     function(ME){lm(as.formula(paste0(ME,"~ Disease")), data = New)
                     }
)

lapply(lm_results,summary)
Res_1<-as.data.frame(coefficients(summary(lm_results[[1]])))
Res_1<- Res_1[!(row.names(Res_1) %in% c("(Intercept)")),]
Res_2<-as.data.frame(coefficients(summary(lm_results[[2]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res_1, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results[[3]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results[[4]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results[[5]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results[[6]])))
Res_2<-Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results[[7]])))
Res_2<-Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results[[8]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results[[9]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results[[10]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results[[11]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res, Res_2)

colnames(Res)<- c("Beta","Std","tvalue","pvalue")
Res$ME<- ME_names

library(dplyr)
Res$fdr <- p.adjust(Res$pvalue, method = "fdr")
Res$Mlogfdr <- (-log10(Res$fdr))
write.csv(Res, file="Res_linearmodelExmodules.csv")

#ModuleKME
ModulekME <- seurat_obj@misc$COVID_sc$wgcna_modules
write.csv(ModulekME, file="hdWGCNA_Excitatory_modules_kME.csv")

#Plot M3 module celltype
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

p <- VlnPlot(
  seurat_obj,
  features = 'ExNeuN3',
  group.by = 'cellID',
  pt.size = 0 # don't show actual data points
)

# add box-and-whisker plots on top:
p <- p + geom_boxplot(width=.25, fill='white')

# change axis labels and remove legend:
p <- p + xlab('') + ylab('ME') + NoLegend()

# plot output
pdf("Scale_Excitatory_M3_celltype.pdf", height=4.5, width=5)
p
dev.off()

seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))
saveRDS(seurat_obj, file='hdWGCNA_Scaledata_Excitatory_neuron_module.rds')

#Celltype_dimention.pdf and Celltype_dimention_2.pdf are used for Figure 3A
#Excitatory_Dendrogram_scaledata.pdf is used for Figure 3B
#ModuleUMAPplot_Excitatory_with_no_hubgenes.pdf is used for Figure 3C
#Scale_Excitatory_M3_celltype.pdf is used for Figure 3E
#ME3 hubgenenetwork plot from ModuleNetwork plot folder is used for Figure 3G
