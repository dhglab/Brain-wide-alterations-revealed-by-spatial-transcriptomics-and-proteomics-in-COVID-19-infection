library(grid)
library(tidyverse)
library(shadowtext)
library(ggplot2)

setwd("./Fig2Sourcedata")
#load GO_inputfiles_barplots from Figure2_3F.rda

#barplots for Figure2
GO_inputfiles_barplots$MLogFDR <- -GO_inputfiles_barplots$LogFDR

A2 <- subset(GO_inputfiles_barplots, Module %in% "turquoise")
tiff("A2.tiff", units="in", width = 5, height = 2, res=300)
ggplot(A2, aes(x = reorder(Description, +MLogFDR), y = MLogFDR)) +
  geom_col(width = 0.5, color='royalblue',fill='royalblue')+
  geom_hline(yintercept = 1.3, color = "red", linetype="solid")+
  coord_flip() +
  labs(y = "-Log(FDR p value)", x = "Top GO")
dev.off()

C2 <- subset(GO_inputfiles_barplots, Module %in% "greenyellow")
tiff("C2.tiff", units="in", width = 5, height = 1.7, res=300)
ggplot(C2, aes(x = reorder(Description, +MLogFDR), y = MLogFDR)) +
  geom_col(width = 0.5, color='royalblue',fill='royalblue')+
  geom_hline(yintercept = 1.3, color = "red", linetype="solid")+
  coord_flip() +
  labs(y = "-Log(FDR p value)", x = "Top GO")
dev.off()

E2 <- subset(GO_inputfiles_barplots, Module %in% "green")
tiff("E2.tiff", units="in", width = 4, height = 1.5, res=300)
ggplot(E2, aes(x = reorder(Description, +MLogFDR), y = MLogFDR)) +
  geom_col(width = 0.5, color='royalblue',fill='royalblue')+
  geom_hline(yintercept = 1.3, color = "red", linetype="solid")+
  coord_flip() +
  labs(y = "-Log(FDR p value)", x = "Top GO")
dev.off()

G2 <- subset(GO_inputfiles_barplots, Module %in% "blue")
tiff("G2.tiff", units="in", width = 5, height = 1.7, res=300)
ggplot(G2, aes(x = reorder(Description, +MLogFDR), y = MLogFDR)) +
  geom_col(width = 0.5, color='royalblue',fill='royalblue')+
  geom_hline(yintercept = 1.3, color = "red", linetype="solid")+
  coord_flip() +
  labs(y = "-Log(FDR p value)", x = "Top GO")
dev.off() 

setwd("./Fig3Sourcedata")
#load RegressedGOExNeuNM3 from Figure2_3F.rda
#barplots for Figure3F
RegressedGOExNeuNM3$fdr <- -RegressedGOExNeuNM3$`Log(q-value)`
tiff("Fig3_F.tiff", units="in", width= 4.5, height=3.5, res=300)
ggplot(RegressedGOExNeuNM3, aes(x = reorder(Description, +fdr), y=fdr)) + 
  geom_bar(stat = "identity", fill = "royalblue",width=0.7) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "red") +
  coord_flip()
dev.off()

