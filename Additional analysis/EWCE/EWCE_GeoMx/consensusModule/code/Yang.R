setwd("./Yang_COVID")

library(EWCE); library(ggplot2); library(cowplot); library(limma); 
library(readxl); library(biomaRt) ; library(tidyverse)
library(readr)
load("./Yang_COVID/SizeBiased_Level2_CellType/CellTypeData_Yang_anno_ctl.rda")

All_frontalGMPN_modules <- read_csv("./All_frontalGMPN_modules.csv")
avg_method="SizeBiased"

marker_vec=ctd

test_vec = All_frontalGMPN_modules$black
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

save(EWCE_enrichment,EWCE_contingency,background,file= "black-Enrichment-P_SizeBiased.RData")

black <-full_results$results
black$MEcolors = c("MEblack")
black$FDR <- p.adjust(black$p, method = "fdr")
black_p <- filter(black, FDR<0.05)
save (black, black_p, file="black_celltype.Rdata")
#########################################################################
test_vec = All_frontalGMPN_modules$blue
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

save(EWCE_enrichment,EWCE_contingency,background,file= "blue-Enrichment-P_SizeBiased.RData")

blue <-full_results$results
blue$MEcolors = c("MEblue")
blue$FDR <- p.adjust(blue$p, method = "fdr")
blue_p <- filter(blue, FDR<0.05)
save (blue, blue_p, file="blue_celltype.Rdata")
#######################################################
test_vec = All_frontalGMPN_modules$brown
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

save(EWCE_enrichment,EWCE_contingency,background,file= "brown-Enrichment-P_SizeBiased.RData")

brown <-full_results$results
brown$MEcolors = c("MEbrown")
brown$FDR <- p.adjust(brown$p, method = "fdr")
brown_p <- filter(brown, FDR<0.05)
save (brown, brown_p, file="brown_celltype.Rdata")
#############################################################
test_vec = All_frontalGMPN_modules$cyan
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

save(EWCE_enrichment,EWCE_contingency,background,file= "cyan-Enrichment-P_SizeBiased.RData")

cyan <-full_results$results
cyan$MEcolors = c("MEcyan")
cyan$FDR <- p.adjust(cyan$p, method = "fdr")
cyan_p <- filter(cyan,FDR<0.05)
save (cyan, cyan_p, file="cyan_celltype.Rdata")
##########################################################################
test_vec = All_frontalGMPN_modules$darkgreen
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

save(EWCE_enrichment,EWCE_contingency,background,file= "darkgreen-Enrichment-P_SizeBiased.RData")

darkgreen <-full_results$results
darkgreen$MEcolors = c("MEdarkgreen")
darkgreen$FDR <- p.adjust(darkgreen$p, method = "fdr")
darkgreen_p <- filter(darkgreen,FDR<0.05)
save (darkgreen, darkgreen_p, file="darkgreen_celltype.Rdata")
###########################################################################
test_vec = All_frontalGMPN_modules$darkred
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

save(EWCE_enrichment,EWCE_contingency,background,file= "darkred-Enrichment-P_SizeBiased.RData")

darkred <-full_results$results
darkred$MEcolors = c("MEdarkred")
darkred$FDR <- p.adjust(darkred$p, method = "fdr")
darkred_p <- filter(darkred, FDR<0.05)
save (darkred, darkred_p, file="darkred_celltype.Rdata")
########################################################################
test_vec = All_frontalGMPN_modules$green
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

save(EWCE_enrichment,EWCE_contingency,background,file= "green-Enrichment-P_SizeBiased.RData")

green <-full_results$results
green$MEcolors = c("MEgreen")
green$FDR <- p.adjust(green$p, method = "fdr")
green_p <- filter(green,FDR<0.05)
save (green, green_p, file="green_celltype.Rdata")
#########################################################################
test_vec = All_frontalGMPN_modules$greenyellow
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

save(EWCE_enrichment,EWCE_contingency,background,file= "greenyellow-Enrichment-P_SizeBiased.RData")

greenyellow <-full_results$results
greenyellow$MEcolors = c("MEgreenyellow")
greenyellow$FDR <- p.adjust(greenyellow$p, method = "fdr")
greenyellow_p <- filter(greenyellow, FDR<0.05)
save (greenyellow, greenyellow_p, file="greenyellow_celltype.Rdata")
#################################################################################
test_vec = All_frontalGMPN_modules$grey60
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

save(EWCE_enrichment,EWCE_contingency,background,file= "grey60-Enrichment-P_SizeBiased.RData")

grey60 <-full_results$results
grey60$MEcolors = c("MEgrey60")
grey60$FDR <- p.adjust(grey60$p, method = "fdr")
grey60_p <- filter(grey60,FDR<0.05)
save (grey60, grey60_p, file="grey60_celltype.Rdata")
##############################################################################
test_vec = All_frontalGMPN_modules$lightcyan
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

save(EWCE_enrichment,EWCE_contingency,background,file= "lightcyan-Enrichment-P_SizeBiased.RData")

lightcyan <-full_results$results
lightcyan$MEcolors = c("MElightcyan")
lightcyan$FDR <- p.adjust(lightcyan$p, method = "fdr")
lightcyan_p <- filter(lightcyan, FDR<0.05)
save (lightcyan, lightcyan_p, file="lightcyan_celltype.Rdata")
############################################################################
test_vec = All_frontalGMPN_modules$lightgreen
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

save(EWCE_enrichment,EWCE_contingency,background,file= "lightgreen-Enrichment-P_SizeBiased.RData")

lightgreen <-full_results$results
lightgreen$MEcolors = c("MElightgreen")
lightgreen$FDR <- p.adjust(lightgreen$p, method = "fdr")
lightgreen_p <- filter(lightgreen, FDR<0.05)
save (lightgreen, lightgreen_p, file="lightgreen_celltype.Rdata")
#######################################################################
test_vec = All_frontalGMPN_modules$lightyellow
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

save(EWCE_enrichment,EWCE_contingency,background,file= "lightyellow-Enrichment-P_SizeBiased.RData")

lightyellow <-full_results$results
lightyellow$MEcolors = c("MElightyellow")
lightyellow$FDR <- p.adjust(lightyellow$p, method = "fdr")
lightyellow_p <- filter(lightyellow, FDR<0.05)
save (lightyellow, lightyellow_p, file="lightyellow_celltype.Rdata")
###########################################################################
test_vec = All_frontalGMPN_modules$magenta
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

save(EWCE_enrichment,EWCE_contingency,background,file= "magenta-Enrichment-P_SizeBiased.RData")

magenta <-full_results$results
magenta$MEcolors = c("MEmagenta")
magenta$FDR <- p.adjust(magenta$p, method = "fdr")
magenta_p <- filter(magenta,FDR<0.05)
save (magenta, magenta_p, file="magenta_celltype.Rdata")
##############################################################################
test_vec = All_frontalGMPN_modules$midnightblue
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

save(EWCE_enrichment,EWCE_contingency,background,file= "midnightblue-Enrichment-P_SizeBiased.RData")

midnightblue <-full_results$results
midnightblue$MEcolors = c("MEmidnightblue")
midnightblue$FDR <- p.adjust(midnightblue$p, method = "fdr")
midnightblue_p <- filter(midnightblue,FDR<0.05)
save (midnightblue, midnightblue_p, file="midnightblue_celltype.Rdata")
###############################################################################
test_vec = All_frontalGMPN_modules$pink
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

save(EWCE_enrichment,EWCE_contingency,background,file= "pink-Enrichment-P_SizeBiased.RData")

pink <-full_results$results
pink$MEcolors = c("MEpink")
pink$FDR <- p.adjust(pink$p, method = "fdr")
pink_p <- filter(pink,FDR<0.05)
save (pink, pink_p, file="pink_celltype.Rdata")
#########################################################
test_vec = All_frontalGMPN_modules$purple
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

save(EWCE_enrichment,EWCE_contingency,background,file= "purple-Enrichment-P_SizeBiased.RData")

purple <-full_results$results
purple$MEcolors = c("MEpurple")
purple$FDR <- p.adjust(purple$p, method = "fdr")
purple_p <- filter(purple, FDR<0.05)
save (purple, purple_p, file="purple_celltype.Rdata")
#########################################################################
test_vec = All_frontalGMPN_modules$red
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

save(EWCE_enrichment,EWCE_contingency,background,file= "red-Enrichment-P_SizeBiased.RData")

red <-full_results$results
red$MEcolors = c("MEred")
red$FDR <- p.adjust(red$p, method = "fdr")
red_p <- filter(red, FDR<0.05)
save (red, red_p, file="red_celltype.Rdata")
###########################################################################
test_vec = All_frontalGMPN_modules$royalblue
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

save(EWCE_enrichment,EWCE_contingency,background,file= "royalblue-Enrichment-P_SizeBiased.RData")

royalblue <-full_results$results
royalblue$MEcolors = c("MEroyalblue")
royalblue$FDR <- p.adjust(royalblue$p, method = "fdr")
royalblue_p <- filter(royalblue, FDR<0.05)
save (royalblue, royalblue_p, file="royalblue_celltype.Rdata")
###############################################################################
test_vec = All_frontalGMPN_modules$salmon
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

save(EWCE_enrichment,EWCE_contingency,background,file= "salmon-Enrichment-P_SizeBiased.RData")

salmon <-full_results$results
salmon$MEcolors = c("MEsalmon")
salmon$FDR <- p.adjust(salmon$p, method = "fdr")
salmon_p <- filter(salmon,FDR<0.05)
save (salmon, salmon_p, file="salmon_celltype.Rdata")
##############################################################################
test_vec = All_frontalGMPN_modules$tan
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

save(EWCE_enrichment,EWCE_contingency,background,file= "tan-Enrichment-P_SizeBiased.RData")

tan <-full_results$results
tan$MEcolors = c("MEtan")
tan$FDR <- p.adjust(tan$p, method = "fdr")
tan_p <- filter(tan, FDR<0.05)
save (tan, tan_p, file="tan_celltype.Rdata")
#############################################################################
test_vec = All_frontalGMPN_modules$turquoise
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

save(EWCE_enrichment,EWCE_contingency,background,file= "turquoise-Enrichment-P_SizeBiased.RData")

turquoise <-full_results$results
turquoise$MEcolors = c("MEturquoise")
turquoise$FDR <- p.adjust(turquoise$p, method = "fdr")
turquoise_p <- filter(turquoise, FDR<0.05)
save (turquoise, turquoise_p, file="turquoise_celltype.Rdata")
write.csv(turquoise, file="turquoise_celltype_yang_level2.csv")
###############################################################################
test_vec = All_frontalGMPN_modules$yellow
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

save(EWCE_enrichment,EWCE_contingency,background,file= "yellow-Enrichment-P_SizeBiased.RData")

yellow <-full_results$results
yellow$MEcolors = c("MEyellow")
yellow$FDR <- p.adjust(yellow$p, method = "fdr")
yellow_p <- filter(yellow,FDR<0.05)
save (yellow, yellow_p, file="yellow_celltype.Rdata")

COVIDmodule_Yang <- rbind(black_p, blue_p,	brown_p,	cyan_p,	darkgreen_p,
                     darkred_p,	green_p,	greenyellow_p,	grey60_p,	lightcyan_p,	lightgreen_p,
                     lightyellow_p,	magenta_p,	midnightblue_p,	pink_p,	purple_p,	red_p,	
                     royalblue_p,	salmon_p,	tan_p,	turquoise_p,	yellow_p)

save(COVIDmodule_Yang, file= "COVIDmodule_Yang_level2.Rdata")
write.csv(COVIDmodule_Yang, file="COVIDModules_celltype_Yang_level2.csv")

#Results are input files for Extended Figure 4D