
setwd("./NRGNEWCEwith Wamsley and Bakken")

load("./CellTypeData_Yang_anno_ctl.rda")
Yanglevel1specificity <- ctd[[1]]$specificity
Yanglevel2specificity <- ctd[[2]]$specificity

#FoundNRGN_celltype_markers
YangNRGNL1 <- subset(Yangsclevel1typespecificity, `NRGN neuron` > 0.25)
YangNRGNL2 <- subset(Yangsclevel2typespecificity, `NRGN neuron` > 0.25)

library(EWCE); library(ggplot2); library(cowplot); library(limma); library(readxl); library(biomaRt) ; library(tidyverse)

load("./SizeBiased_Level2_CellType/CellTypeData_Bakken_anno_ctrl.rda")

avg_method="SizeBiased"

markers=ctd

test = YangNRGNL2$...1
level=2 # <- ie. Use level 1 annotations (i.e. Interneurons)
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

save(EWCE_enrichment,EWCE_contingency,background,file= "NRGN_BakkenBakkenEWCE-Enrichment-P_SizeBiased.RData")

NRGN_Bakken_cell <-full_results$results
NRGN_Bakken_cell$area <- c("Bakken")
NRGN_Bakken_cell$FDR <- p.adjust(NRGN_Bakken_cell$p, method = "fdr")
save(NRGN_Bakken_cell, file="NRGN_Bakken_cell.Rdata")
##################################################################################################
load("./CellTypeData_ASDobect_anno_ctl.rda")
avg_method="SizeBiased"

markers=ctd

test = YangNRGNL2$...1
level=2 # <- ie. Use level 1 annotations (i.e. Interneurons)
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

save(EWCE_enrichment,EWCE_contingency,background,file= "NRGN_WamsleyL2BakkenEWCE-Enrichment-P_SizeBiased.RData")

NRGN_WamsleyL2_cell <-full_results$results
NRGN_WamsleyL2_cell$area <- c("WamsleyL2")
NRGN_WamsleyL2_cell$FDR <- p.adjust(NRGN_WamsleyL2_cell$p, method = "fdr")
save(NRGN_WamsleyL2_cell, file="NRGN_WamsleyL2_cell.Rdata")
###################################################################################################
test = YangNRGNL1$...1
level=1 # <- ie. Use level 1 annotations (i.e. Interneurons)
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

save(EWCE_enrichment,EWCE_contingency,background,file= "NRGN_WamsleyL1BakkenEWCE-Enrichment-P_SizeBiased.RData")

NRGN_WamsleyL1_cell <-full_results$results
NRGN_WamsleyL1_cell$area <- c("WamsleyL1")
NRGN_WamsleyL1_cell$FDR <- p.adjust(NRGN_WamsleyL1_cell$p, method = "fdr")
save(NRGN_WamsleyL1_cell, file="NRGN_WamsleyL1_cell.Rdata")

NRGN_celltype <- rbind(NRGN_WamsleyL1_cell, NRGN_WamsleyL2_cell, NRGN_Bakken_cell)
write.csv(NRGN_celltype, file="NRGN_celltype.csv")

#replace negative sd_from_mean with 0 for plotting 
NRGN_celltype_2 <- NRGN_celltype
NRGN_celltype_2$sdfrommean <- replace(NRGN_celltype$sd_from_mean, NRGN_celltype$sd_from_mean < 0, 0) 

NRGN_celltype_2_WamsleyL1 <- subset(NRGN_celltype_2, area %in% c("WamsleyL1"))
NRGN_celltype_2_WamsleyL2 <- subset(NRGN_celltype_2, area %in% c("WamsleyL2"))
NRGN_celltype_2_Bakken <- subset(NRGN_celltype_2, area %in% c("Bakken"))

#Extended Figure1_D1
pdf("NRGNcelltype_WamsleyL1.pdf", height = 4, width=3)
ggplot(NRGN_celltype_2_WamsleyL1, aes(x=CellType, y=sdfrommean))+
  geom_bar(stat='identity', fill="#56B4E9")+
  facet_wrap(~area,  ncol=1)+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_hline(yintercept=2.5, linetype="dashed", color = "red")+
  theme_cowplot()
dev.off()

#Extended Figure1_D2
pdf("NRGNcelltype_Bakken.pdf", height=4, width=7.5)
ggplot(NRGN_celltype_2_Bakken, aes(x=CellType, y=sdfrommean))+
  geom_bar(stat='identity', fill="#56B4E9")+
  facet_wrap(~area,  ncol=1)+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_hline(yintercept=2.5, linetype="dashed", color = "red")+
  theme_cowplot()
dev.off()

#Extended Figure1_D3
pdf("NRGNcelltype_WamsleyL2_2.pdf", height=4, width=12)
ggplot(NRGN_celltype_2_WamsleyL2, aes(x=CellType, y=sdfrommean))+
  geom_bar(stat='identity', fill="#56B4E9")+
  facet_wrap(~area,  ncol=1)+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_hline(yintercept=2.5, linetype="dashed", color = "red")+
  theme_cowplot()
dev.off()


