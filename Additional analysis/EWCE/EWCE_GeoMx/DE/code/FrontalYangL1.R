
setwd("./Yang_COVID")
load("./SizeBiased_Level2_CellType/CellTypeData_Yang_anno_ctrl.rda")
load("./FrontalGM56_rank_up_down.RData")
load("./FrontalWM56_rank_up_down.RData")

set.seed(12211991)
options(stringsAsFactors = F)

#install.packages('RNOmni')

library(EWCE); library(ggplot2); library(cowplot); library(limma); library(readxl); library(biomaRt) ; library(tidyverse)

avg_method="SizeBiased"

markers=ctd

######FrontalGM
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

save(EWCE_enrichment,EWCE_contingency,background,file= "FrontalGM56_downYangEWCE-Enrichment-P_SizeBiased.RData")

FrontalGM56_down_cell <-full_results$results
FrontalGM56_down_cell$area <- c("FrontalGM")
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

save(EWCE_enrichment,EWCE_contingency,background,file= "FrontalGM56_up_YangEWCE-Enrichment-P_SizeBiased.RData")

FrontalGM56_up_cell <-full_results$results
FrontalGM56_up_cell$area <- c("FrontalGM")
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

save(EWCE_enrichment,EWCE_contingency,background,file= "FrontalWM56_downYangEWCE-Enrichment-P_SizeBiased.RData")

FrontalWM56_down_cell <-full_results$results
FrontalWM56_down_cell$area <- c("FrontalWM")
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

save(EWCE_enrichment,EWCE_contingency,background,file= "FrontalWM56_up_YangEWCE-Enrichment-P_SizeBiased.RData")

FrontalWM56_up_cell <-full_results$results
FrontalWM56_up_cell$area <- c("FrontalWM")
FrontalWM56_up_cell$trend <- c("up")
FrontalWM56_up_cell$FDR <- p.adjust(FrontalWM56_up_cell$p, method = "fdr")
save(FrontalWM56_up_cell, file="FrontalWM56_up_cell.Rdata")

COVID_down_Yang <- rbind (FrontalGM56_down_cell, FrontalWM56_down_cell)
save(COVID_down_Yang, file = "COVID_down_Yang_level1.Rdata")
write.csv(COVID_down_Yang, file = "COVID_down_Yang_level1.csv")

COVID_up_Yang <- rbind (FrontalGM56_up_cell, FrontalWM56_up_cell)
save(COVID_up_Yang, file = "COVID_up_Yang_level1.Rdata")
write.csv(COVID_up_Yang, file = "COVID_up_Yang_level1.csv")

#The results are input files for Figure 1_G1