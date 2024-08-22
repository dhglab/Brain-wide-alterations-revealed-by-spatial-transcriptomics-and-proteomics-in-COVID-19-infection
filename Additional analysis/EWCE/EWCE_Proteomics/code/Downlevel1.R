install.packages('matrixStats')
getOption("repos")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

load("/CellTypeData_ASDobect_anno_ctl.rda")

library(EWCE); library(ggplot2);library(tidyverse)
install.packages('RNOmni', lib = "/u/home/t/tzhang/R/APPTAINER/h2-rstudio_4.1.0")
marker_vec=ctd

test_vec = All_proteomics_Down$Frontal_Down
level_vec=c(1,2)
markers= ctd

test=names(test_vec)
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test_vec
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "Frontal_Down_Enrichment-P_SizeBiased.RData")

Frontal_Down <-full_results$results
Frontal_Down$area = c("Frontal")
Frontal_Down_p <- filter(Frontal_Down, p<0.05)
save (Frontal_Down, Frontal_Down_p, file="Frontal_Down_celltype_proteomics.Rdata")
#####################################################################################
test_vec = All_proteomics_Down$Temporal_Down
level_vec=c(1,2)
markers= ctd

test=names(test_vec)
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test_vec
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "Temporal_Down_Enrichment-P_SizeBiased.RData")

Temporal_Down <-full_results$results
Temporal_Down$area = c("Temporal")
Temporal_Down_p <- filter(Temporal_Down, p<0.05)
save (Temporal_Down, Temporal_Down_p, file="Temporal_Down_celltype_proteomics.Rdata")
#########################################################################################
test_vec = All_proteomics_Down$Hippocampus_Down
level_vec=c(1,2)
markers= ctd

test=names(test_vec)
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test_vec
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "Hippocampus_Down_Enrichment-P_SizeBiased.RData")

Hippocampus_Down <-full_results$results
Hippocampus_Down$area = c("Hippocampus")
Hippocampus_Down_p <- filter(Hippocampus_Down, p<0.05)
save (Hippocampus_Down, Hippocampus_Down_p, file="Hippocampus_Down_celltype_proteomics.Rdata")
#####################################################################################
test_vec = All_proteomics_Down$Midbrain_Down
level_vec=c(1,2)
markers= ctd

test=names(test_vec)
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test_vec
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "Midbrain_Down_Enrichment-P_SizeBiased.RData")

Midbrain_Down <-full_results$results
Midbrain_Down$area = c("Midbrain")
Midbrain_Down_p <- filter(Midbrain_Down, p<0.05)
save (Midbrain_Down, Midbrain_Down_p, file="Midbrain_Down_celltype_proteomics.Rdata")
#########################################################################################
test_vec = All_proteomics_Down$`Basal Ganglia_Down`
level_vec=c(1,2)
markers= ctd

test=names(test_vec)
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test_vec
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "BG_Down_Enrichment-P_SizeBiased.RData")

BG_Down <-full_results$results
BG_Down$area = c("BG")
BG_Down_p <- filter(BG_Down, p<0.05)
save (BG_Down, BG_Down_p, file="BG_Down_celltype_proteomics.Rdata")
###############################################################################
test_vec = All_proteomics_Down$Pons_Down
level_vec=c(1,2)
markers= ctd

test=names(test_vec)
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test_vec
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "Pons_Down_Enrichment-P_SizeBiased.RData")

Pons_Down <-full_results$results
Pons_Down$area = c("Pons")
Pons_Down_p <- filter(Pons_Down, p<0.05)
save (Pons_Down, Pons_Down_p, file="Pons_Down_celltype_proteomics.Rdata")
###################################################################################
test_vec = All_proteomics_Down$Occipital_Down
level_vec=c(1,2)
markers= ctd

test=names(test_vec)
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test_vec
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "Occipital_Down_Enrichment-P_SizeBiased.RData")

Occipital_Down <-full_results$results
Occipital_Down$area = c("Occipital")
Occipital_Down_p <- filter(Occipital_Down, p<0.05)
save (Occipital_Down, Occipital_Down_p, file="Occipital_Down_celltype_proteomics.Rdata")
##########################################################################################
test_vec = All_proteomics_Down$Thalamus_Down
level_vec=c(1,2)
markers= ctd

test=names(test_vec)
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test_vec
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "Thalamus_Down_Enrichment-P_SizeBiased.RData")

Thalamus_Down <-full_results$results
Thalamus_Down$area = c("Thalamus")
Thalamus_Down_p <- filter(Thalamus_Down, p<0.05)
save (Thalamus_Down, Thalamus_Down_p, file="Thalamus_Down_celltype_proteomics.Rdata")
##########################################################################################
Proteomics_Down <- rbind(Frontal_Down_p, Temporal_Down_p, Hippocampus_Down_p, BG_Down_p, Midbrain_Down_p, 
                       Pons_Down_p, Occipital_Down_p, Thalamus_Down_p)
Proteomics_Down$FDR <- p.adjust(Proteomics_Down$p, method = "fdr")
write.csv(Proteomics_Down, file="Proteomics_Down_EWCE_level1.csv")
