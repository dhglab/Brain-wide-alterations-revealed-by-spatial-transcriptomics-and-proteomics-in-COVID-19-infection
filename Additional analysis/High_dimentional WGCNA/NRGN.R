library(Seurat)
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

setwd("./NRGN")
#load seurat_obj_COVID_metacell_regress_newPCA.rds from ExNeuN processing

seurat_obj <- seurat_obj_COVID_metacell_regress_newPCA

seurat_obj <- SetDatExpr(
  seurat_obj,
  assay = "RNA",
  slot = "data",
  group_name = "NRGN neuron", 
  group.by='cellID'
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed'
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
pdf("thresh_hold_NRGN.pdf")
wrap_plots(plot_list, ncol=2)
dev.off()

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=8,
  setDatExpr=FALSE,
  tom_name = 'NRGN neuron_regress' # name of the topoligical overlap matrix written to disk
)

#Supplementary Fig2A
pdf("NRGN_Dendrogram_scaledata.pdf", height=4, width=6)
PlotDendrogram(seurat_obj, main='NRGN neuron_Dendrogram')
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
  group.by = 'cellID', group_name = 'NRGN neuron',
  sparse=FALSE
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "NRGN-M"
)

# plot genes ranked by kME for each module
pdf("kME_Scale_NRGN.pdf")
PlotKMEs(seurat_obj, ncol=5)
dev.off()

# get the module assignment table:
modules <- GetModules(seurat_obj)

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)

saveRDS(seurat_obj, file='hdWGCNA_Scaledata_NRGN_neuron.rds')
############################################################################################

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
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)


# stitch together with patchwork
pdf("Celltype_dimention.pdf", height=4, width=5.5)
DimPlot(seurat_obj, reduction = "umap")
dev.off()

pdf("NRGN_Modulefeature_scale.pdf")
wrap_plots(plot_list, ncol=4)
dev.off()

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = FALSE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
pdf("NRGN_Modulefeature_hubscores.pdf")
wrap_plots(plot_list, ncol=2)
dev.off()

# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
saveRDS(seurat_obj, file='hdWGCNA_object_with_NRGN_meta.rds')

# plot with Seurat's DotPlot function
# flip the x/y axes, rotate the axis labels, and change color scheme:
pdf("Dotplot_celltype_module_features_t.pdf", height = 4.5, width=6)
DotPlot(seurat_obj, features=mods, group.by = 'cellID')+ 
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')
dev.off()
# plot output

# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# network analysis & visualization package:
library(igraph)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)


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
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=25, # neighbors parameter for UMAP
  min_dist=0.4 # min distance between points in UMAP space
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
pdf('ModulefeatureUMAP_unsupervised_scale.pdf')
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()
dev.off()
###worked#######

pdf("ModuleUMAPplot_NRGN_with_hubgenes.pdf", height=12,width=12)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=1 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  )
dev.off()

pdf("ModuleUMAPplot_NRGN_with_hubgenes_0.pdf", height=12,width=12)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=0 ,# how many hub genes to plot per module?
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
pdf("SupervisedUMAPmoduleNRGN.pdf")
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()
dev.off()

#####################Differential module eigengene (DME) analysis#############
# get hMEs from seurat object
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "ME"
)

MEs <- GetMEs(seurat_obj, harmonized=FALSE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

New <- seurat_obj@meta.data[1:35]
New <- cbind(New, MEs)
ME_names <- colnames(MEs)

unique(New$Disease)

New$Disease <- fct_relevel(New$Disease, "Non-viral")

Expr <- as.matrix(MEs)

New$nNuclei <- scale(as.numeric(New$nNuclei))
New$nFeature_RNA <- scale(as.numeric(New$nFeature_RNA))
formula <- as.formula(paste0("Expr ~ Disease + Age_s + Sex + Batch + percent.mt_s + nCount_RNA_s + nNuclei + nFeature_RNA"))
lm <- lm(formula, data = New)
Coefficentlist <- coef(summary(lm))
Coefficentlist[[1]]
lm_results <- lapply(mods, 
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
Res_2<-as.data.frame(coefficients(summary(lm_results[[12]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results[[13]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results[[14]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results[[15]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results[[16]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results[[17]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results[[18]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results[[19]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)")),]
Res<-rbind(Res, Res_2)

colnames(Res)<- c("Beta","Std","tvalue","pvalue")
Res$ME<- ME_names

library(dplyr)
Res$fdr <- p.adjust(Res$pvalue, method = "fdr")
Res$Mlogfdr <- (-log10(Res$fdr))
write.csv(Res, file="Res_linearmodelNRGNmodules.csv")
#input file for Supplementary Fig2B

sc_NRGN_modules_kM <- seurat_obj@misc$COVID_sc$wgcna_modules
write.csv(sc_NRGN_modules_kM, file="sc_NRGN_modules_kM.csv")

library(enrichR)
library(GeneOverlap)
# using the cowplot theme for ggplot
theme_set(theme_gray())
set.seed(12345)

# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021','KEGG2021')

# perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 200 # number of genes per module to test
)

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "NRGN-M"
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)
save(enrich_df, file="enrich_df_NRGN_scaleregress.rda")

# enrichr dotplot Supplementary Fig2C
pdf("Enrich_NRGN_scaleregress.pdf", height=8, width=9.5)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", 
  database = "GO_Biological_Process_2021", 
  n_terms=1 
)
dev.off()
