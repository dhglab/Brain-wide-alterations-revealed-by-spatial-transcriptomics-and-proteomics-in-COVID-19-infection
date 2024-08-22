
setwd("./WamsleyASD")
load("./SizeBiased_Level2_CellType/CellTypeData_ASDobect_anno_ctl.rda")
load("./DE_EWCE_input.rda")

set.seed(12211991)
options(stringsAsFactors = F)

#install.packages('RNOmni')

library(EWCE); library(ggplot2); library(cowplot); library(limma); library(readxl); library(biomaRt) ; library(tidyverse)

avg_method="SizeBiased"

markers=ctd

########FrontalGM
test = FrontalGM56_down_R$`Target name`
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "FrontalGM56_downWamsleyEWCE-Enrichment-P_SizeBiased.RData")

FrontalGM56_down_cell <-full_results$results
FrontalGM56_down_cell$area <- c("FrontalGM56")
FrontalGM56_down_cell$trend <- c("down")
FrontalGM56_down_cell$FDR <- p.adjust(FrontalGM56_down_cell$p, method = "fdr")
save(FrontalGM56_down_cell, file="FrontalGM56_down_cell.Rdata")

########

test = FrontalGM56_up_R$`Target name`
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "FrontalGM56_up_WamsleyEWCE-Enrichment-P_SizeBiased.RData")

FrontalGM56_up_cell <-full_results$results
FrontalGM56_up_cell$area <- c("FrontalGM56")
FrontalGM56_up_cell$trend <- c("up")
FrontalGM56_up_cell$FDR <- p.adjust(FrontalGM56_up_cell$p, method = "fdr")
save(FrontalGM56_up_cell, file="FrontalGM56_up_cell.Rdata")

#######FrontalWM
test = FrontalWM56_down_R$`Target name`
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "FrontalWM56_downWamsleyEWCE-Enrichment-P_SizeBiased.RData")

FrontalWM56_down_cell <-full_results$results
FrontalWM56_down_cell$area <- c("FrontalWM56")
FrontalWM56_down_cell$trend <- c("down")
FrontalWM56_down_cell$FDR <- p.adjust(FrontalWM56_down_cell$p, method = "fdr")
save(FrontalWM56_down_cell, file="FrontalWM56_down_cell.Rdata")

test = FrontalWM56_up_R$`Target name`
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "FrontalWM56_up_WamsleyEWCE-Enrichment-P_SizeBiased.RData")

FrontalWM56_up_cell <-full_results$results
FrontalWM56_up_cell$area <- c("FrontalWM56")
FrontalWM56_up_cell$trend <- c("up")
FrontalWM56_up_cell$FDR <- p.adjust(FrontalWM56_up_cell$p, method = "fdr")
save(FrontalWM56_up_cell, file="FrontalWM56_up_cell.Rdata")

#####PN
test = PN46_down_R$`Target name`
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "PN46_downWamsleyEWCE-Enrichment-P_SizeBiased.RData")

PN46_down_cell <-full_results$results
PN46_down_cell$area <- c("PN46")
PN46_down_cell$trend <- c("down")
PN46_down_cell$FDR <- p.adjust(PN46_down_cell$p, method = "fdr")
save(PN46_down_cell, file="PN46_down_cell.Rdata")

test = PN46_up_R$`Target name`
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "PN46_up_WamsleyEWCE-Enrichment-P_SizeBiased.RData")

PN46_up_cell <-full_results$results
PN46_up_cell$area <- c("PN46")
PN46_up_cell$trend <- c("up")
PN46_up_cell$FDR <- p.adjust(PN46_up_cell$p, method = "fdr")
save(PN46_up_cell, file="PN46_up_cell.Rdata")

####PCT
test = PCT46_down_R$`Target name`
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "PCT46_downWamsleyEWCE-Enrichment-P_SizeBiased.RData")

PCT46_down_cell <-full_results$results
PCT46_down_cell$area <- c("PCT46")
PCT46_down_cell$trend <- c("down")
PCT46_down_cell$FDR <- p.adjust(PCT46_down_cell$p, method = "fdr")
save(PCT46_down_cell, file="PCT46_down_cell.Rdata")

test = PCT46_up_R$`Target name`
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "PCT46_up_WamsleyEWCE-Enrichment-P_SizeBiased.RData")

PCT46_up_cell <-full_results$results
PCT46_up_cell$area <- c("PCT46")
PCT46_up_cell$trend <- c("up")
PCT46_up_cell$FDR <- p.adjust(PCT46_up_cell$p, method = "fdr")
save(PCT46_up_cell, file="PCT46_up_cell.Rdata")

##########CST
test = CST46_down_R$`Target name`
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "CST46_downWamsleyEWCE-Enrichment-P_SizeBiased.RData")

CST46_down_cell <-full_results$results
CST46_down_cell$area <- c("CST46")
CST46_down_cell$trend <- c("down")
CST46_down_cell$FDR <- p.adjust(CST46_down_cell$p, method = "fdr")
save(CST46_down_cell, file="CST46_down_cell.Rdata")

test = CST46_up_R$`Target name`
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "CST46_up_WamsleyEWCE-Enrichment-P_SizeBiased.RData")

CST46_up_cell <-full_results$results
CST46_up_cell$area <- c("CST46")
CST46_up_cell$trend <- c("up")
CST46_up_cell$FDR <- p.adjust(CST46_up_cell$p, method = "fdr")
save(CST46_up_cell, file="CST46_up_cell.Rdata")

####Teg
test = Teg46_down_R$`Target name`
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "Teg46_downWamsleyEWCE-Enrichment-P_SizeBiased.RData")

Teg46_down_cell <-full_results$results
Teg46_down_cell$area <- c("Teg46")
Teg46_down_cell$trend <- c("down")
Teg46_down_cell$FDR <- p.adjust(Teg46_down_cell$p, method = "fdr")
save(Teg46_down_cell, file="Teg46_down_cell.Rdata")

test = Teg46_up_R$`Target name`
level=1 
reps=100000

EWCE_enrichment <- EWCE_contingency <- matrix(NA,nrow=dim(ctd[[level]]$specificity)[2],ncol=1)
rownames(EWCE_enrichment) <- rownames(EWCE_contingency) <- colnames(ctd[[level]]$specificity)
colnames(EWCE_enrichment) <- colnames(EWCE_contingency) <- "hits"

hits <- test
background = union(hits,rownames(ctd[[level]]$specificity))

full_results = bootstrap_enrichment_test(sct_data=ctd,hits=hits,bg=background,
                                         genelistSpecies="human",sctSpecies="human",
                                         reps=reps,annotLevel=level)
for(cell in rownames(EWCE_enrichment)){
  idx = which(full_results$results$CellType==cell)
  EWCE_enrichment[cell,1]=full_results$results$p[idx]
  EWCE_contingency[cell,1]=length(hits)}

save(EWCE_enrichment,EWCE_contingency,background,file= "Teg46_up_WamsleyEWCE-Enrichment-P_SizeBiased.RData")

Teg46_up_cell <-full_results$results
Teg46_up_cell$area <- c("Teg46")
Teg46_up_cell$trend <- c("up")
Teg46_up_cell$FDR <- p.adjust(Teg46_up_cell$p, method = "fdr")
save(Teg46_up_cell, file="Teg46_up_cell.Rdata")

COVID_down_Wamsley <- rbind (FrontalGM56_down_cell, FrontalWM56_down_cell, PN46_down_cell,PCT46_down_cell,CST46_down_cell,Teg46_down_cell)
save(COVID_down_Wamsley, file = "COVID_down_Wamsley_level1.Rdata")
COVID_down_Wamsley_fdr <- subset(COVID_down_Wamsley, FDR< 0.05)
write.csv(COVID_down_Wamsley_fdr, file = "COVID_down_Wamsley_level1.csv")

COVID_up_Wamsley <- rbind (FrontalGM56_up_cell, FrontalWM56_up_cell, PN46_up_cell,PCT46_up_cell,CST46_up_cell,Teg46_up_cell)
save(COVID_up_Wamsley, file = "COVID_up_Wamsley_level1.Rdata")
COVID_up_Wamsley_fdr <- subset(COVID_up_Wamsley, FDR< 0.05)
write.csv(COVID_up_Wamsley_fdr, file = "COVID_up_Wamsley_level1.csv")

#Results are input files for Figure 1F
