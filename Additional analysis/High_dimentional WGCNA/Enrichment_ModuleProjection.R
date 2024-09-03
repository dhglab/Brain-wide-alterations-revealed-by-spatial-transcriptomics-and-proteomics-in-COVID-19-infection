library(enrichR)
library(GeneOverlap)
theme_set(theme_gray())

# set random seed for reproducibility
set.seed(12345)

#load seurat_obj
seurat_obj <- hdWGCNA_Scaledata_Excitatory_neuron_module.rds
# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, 
  max_genes = 200 
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)
save(enrich_df, file="enrich_df_Ex.rda")

library(Seurat)
# enrichr dotplot, Extended Figure 4E
pdf("Enrich_ExNeuN.pdf", height=6, width=7.5)
EnrichrDotPlot(
  seurat_obj,
  mods = "all", 
  database = "GO_Biological_Process_2021", 
  n_terms=1 
)
dev.off()

#####Projection GeoMX FrontalGM and Pons module to hdWGCNA#########
#load concensuskME_table from GeoMX
consensus_modules = data.frame("gene_name"=consensus_kM_s$Gene,
                               "module"=consensus_kM_s$module,
                               "color"=concensus_kM_s$module)

consensus_modules <- subset(consensus_modules, gene_name %in% rownames(seurat_obj))
rownames(seurat_obj)

consensus_modules <- consensus_modules[match(unique(consensus_modules$gene_name), consensus_modules$gene_name),]
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures=2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- ProjectModules(
  seurat_obj,
  modules = consensus_modules,
  group.by.vars = "Batch",
  seurat_ref = NULL,
  wgcna_name = "None",
  wgcna_name_proj = 'concensusFrontalGMPN'
)

saveRDS(seurat_obj, file = 'consensusFrontalGMPN_projected.rds')

library(ggrastr)
seurat_obj <- consensusFrontalGMPN_projected
seurat_obj <- SetActiveWGCNA(seurat_obj, 'concensusFrontalGMPN')

pdf("Figure3_A2_A3.pdf")
ModuleFeaturePlot(seurat_obj, module_names = c('consensusME_turquoise'))
ModuleFeaturePlot(seurat_obj, module_names = c('consensusME_greenyellow'))
dev.off()

