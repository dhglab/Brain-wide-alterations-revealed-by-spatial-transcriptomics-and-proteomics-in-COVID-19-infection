
setwd("./Yang_COVID")
load("./SizeBiased_Level2_CellType/CellTypeData_Yang_anno_ctrl.rda")
load("./FrontalGM56_rank_up_down.RData")

set.seed(12211991)
options(stringsAsFactors = F)

#install.packages('RNOmni')

library(EWCE); library(ggplot2); library(cowplot); library(limma); library(readxl); library(biomaRt) ; library(tidyverse)

avg_method="SizeBiased"

markers=ctd

test = FrontalGM56_down_R$`Target name`
level=2 
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
save(FrontalGM56_down_cell, file="FrontalGM56_down_cell_Yang_level2.Rdata")
FrontalGM56_down_cell_fdr <- subset(FrongalGM56_down_cell, FDR < 0.05)
write.csv(FrontalGM56_down_cell_fdr, file = "COVID_down_Yang_FrontalGM_level2.csv")

#The results are input files for Figure 1_G1