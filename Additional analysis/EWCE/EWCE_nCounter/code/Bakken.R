load("./SizeBiased_Level2_CellType/CellTypeData_Bakken_anno_ctrl.rda")

library(EWCE); library(ggplot2); library(cowplot); library(limma); 
library(readxl); library(biomaRt) ; library(tidyverse)
library(readr)
setwd("./nCounter/Bakken")
#Differential results of Frontallboe nCounter, res_FrontalP_vs_CMG218_for_EWCE 
nCounter_up_R <- subset(res_FrontalP_vs_CMG218_for_EWCE, rank > 2)
nCounter_down_R <- subset(res_FrontalP_vs_CMG218_for_EWCE, rank < -2)
save(nCounter_down_R, nCounter_up_R, file = "Frontal_nCounterDE_rank.rda")

avg_method="SizeBiased"

marker_vec=ctd

test_vec = nCounter_up_R$Gene
level_vec=c(1,2)
ops = expand.grid(marker_vec,test_vec,level_vec)
markers= ctd

test=names(test_vec)
level=2 
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

save(EWCE_enrichment,EWCE_contingency,background,file= "nCounter_up-Enrichment-P_SizeBiased.RData")

nCounter_up <-full_results$results
nCounter_up$Trend = c("up")
nCounter_up$FDR <- p.adjust(nCounter_up$p, method = "fdr")
nCounter_up_p <- filter(nCounter_up, FDR<0.05)
nCounter_up_p$abundance <- abs(nCounter_up_p$sd_from_mean)
save (nCounter_up, nCounter_up_p, file="nCounter_up_celltype_BakkenleveL2.Rdata")
#########################################################################
test_vec = nCounter_down_R$Gene
level_vec=c(1,2)
ops = expand.grid(marker_vec,test_vec,level_vec)
markers= ctd

test=names(test_vec)
level=2 
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

save(EWCE_enrichment,EWCE_contingency,background,file= "nCounter_down-Enrichment-P_SizeBiased.RData")

nCounter_down <-full_results$results
nCounter_down$Trend = c("down")
nCounter_down$FDR <- p.adjust(nCounter_down$p, method = "fdr")
nCounter_down_p <- filter(nCounter_down, FDR<0.05)
nCounter_down_p$abundance <- -abs(nCounter_down_p$sd_from_mean)
save (nCounter_down, nCounter_down_p, file="nCounter_down_celltype_BakkenleveL2.Rdata")

nCounter_Bakken_L2 <- rbind(nCounter_down_p, nCounter_up_p)

write.csv(nCounter_Bakken_L2, file="nCounter_Bakken_L2.csv")

#The results will be input file for Figure 5A
