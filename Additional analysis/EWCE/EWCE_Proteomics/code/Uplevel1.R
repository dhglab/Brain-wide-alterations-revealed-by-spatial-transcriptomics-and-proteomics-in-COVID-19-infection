load("/CellTypeData_ASDobect_anno_ctl.rda")

library(EWCE); library(ggplot2);library(tidyverse)
install.packages('RNOmni', lib = "/u/home/t/tzhang/R/APPTAINER/h2-rstudio_4.1.0")
marker_vec=ctd

test_vec = All_proteomics_up$Frontal_Up
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

save(EWCE_enrichment,EWCE_contingency,background,file= "Frontal_Up_Enrichment-P_SizeBiased.RData")

Frontal_Up <-full_results$results
Frontal_Up$area = c("Frontal")
Frontal_Up_p <- filter(Frontal_Up, p<0.05)
save (Frontal_Up, Frontal_Up_p, file="Frontal_Up_celltype_proteomics.Rdata")
#####################################################################################
test_vec = All_proteomics_up$Temporal_Up
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

save(EWCE_enrichment,EWCE_contingency,background,file= "Temporal_Up_Enrichment-P_SizeBiased.RData")

Temporal_Up <-full_results$results
Temporal_Up$area = c("Temporal")
Temporal_Up_p <- filter(Temporal_Up, p<0.05)
save (Temporal_Up, Temporal_Up_p, file="Temporal_Up_celltype_proteomics.Rdata")
#########################################################################################
test_vec = All_proteomics_up$Hippocampus_Up
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

save(EWCE_enrichment,EWCE_contingency,background,file= "Hippocampus_Up_Enrichment-P_SizeBiased.RData")

Hippocampus_Up <-full_results$results
Hippocampus_Up$area = c("Hippocampus")
Hippocampus_Up_p <- filter(Hippocampus_Up, p<0.05)
save (Hippocampus_Up, Hippocampus_Up_p, file="Hippocampus_Up_celltype_proteomics.Rdata")
#####################################################################################
test_vec = All_proteomics_up$Midbrain_Up
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

save(EWCE_enrichment,EWCE_contingency,background,file= "Midbrain_Up_Enrichment-P_SizeBiased.RData")

Midbrain_Up <-full_results$results
Midbrain_Up$area = c("Midbrain")
Midbrain_Up_p <- filter(Midbrain_Up, p<0.05)
save (Midbrain_Up, Midbrain_Up_p, file="Midbrain_Up_celltype_proteomics.Rdata")
#########################################################################################
test_vec = All_proteomics_up$`Basal Ganglia_Up`
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

save(EWCE_enrichment,EWCE_contingency,background,file= "BG_Up_Enrichment-P_SizeBiased.RData")

BG_Up <-full_results$results
BG_Up$area = c("BG")
BG_Up_p <- filter(BG_Up, p<0.05)
save (BG_Up, BG_Up_p, file="BG_Up_celltype_proteomics.Rdata")
###############################################################################
test_vec = All_proteomics_up$Pons_Up
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

save(EWCE_enrichment,EWCE_contingency,background,file= "Pons_Up_Enrichment-P_SizeBiased.RData")

Pons_Up <-full_results$results
Pons_Up$area = c("Pons")
Pons_Up_p <- filter(Pons_Up, p<0.05)
save (Pons_Up, Pons_Up_p, file="Pons_Up_celltype_proteomics.Rdata")
###################################################################################
test_vec = All_proteomics_up$Occipital_Up
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

save(EWCE_enrichment,EWCE_contingency,background,file= "Occipital_Up_Enrichment-P_SizeBiased.RData")

Occipital_Up <-full_results$results
Occipital_Up$area = c("Occipital")
Occipital_Up_p <- filter(Occipital_Up, p<0.05)
save (Occipital_Up, Occipital_Up_p, file="Occipital_Up_celltype_proteomics.Rdata")
##########################################################################################
test_vec = All_proteomics_up$Thalamus_Up
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

save(EWCE_enrichment,EWCE_contingency,background,file= "Thalamus_Up_Enrichment-P_SizeBiased.RData")

Thalamus_Up <-full_results$results
Thalamus_Up$area = c("Thalamus")
Thalamus_Up_p <- filter(Thalamus_Up, p<0.05)
save (Thalamus_Up, Thalamus_Up_p, file="Thalamus_Up_celltype_proteomics.Rdata")
##########################################################################################
Proteomics_up <- rbind(Frontal_Up_p, Temporal_Up_p, Hippocampus_Up_p, BG_Up_p, Midbrain_Up_p, 
                       Pons_Up_p, Occipital_Up_p, Thalamus_Up_p)
Proteomics_up$FDR <- p.adjust(Proteomics_up$p, method = "fdr")
write.csv(Proteomics_up, file="Proteomics_up_EWCE_level1.csv")
