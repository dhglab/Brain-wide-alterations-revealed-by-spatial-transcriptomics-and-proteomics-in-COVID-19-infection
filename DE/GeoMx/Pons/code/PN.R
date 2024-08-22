library(dplyr)
library(ggplot2)
library(EDASeq)
library(RColorBrewer)
library(RUVSeq)
library(Hmisc)
library(corrplot)
library(matrixStats)
library(clusterProfiler)
library(WGCNA)
library(Biobase)
library(PCAtools)
library('DESeq2')
library(devtools)

#Load raw(QC_remove non integer), targets, and differential expression from GeoDSP profiler, PNFDRQC_integers, PNQC, TargetPN, PN4_6DEG

#check rownames of target and colnames of matrix
all(rownames(TargetPN) == colnames(PNFDRQC)) 

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = PNFDRQC,
                              colData = TargetPN,
                              design = ~ Disease)
dds

vsd <- vst(dds, blind=TRUE)
boxplot(assay(vsd), las=2, main="vsd")
plotPCA(vsd, intgroup = c("Disease"))

mat <- assay(vsd)
mm <- model.matrix(~Disease, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$LOT, design=mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = c("Disease"))

pca <- pca(mat, metadata = TargetPN)

#Supplementary Figure 1A and 1B
tiff("PN5+6_PCA.tiff", units="in", width= 6.5, height=4, res=300)
biplot(pca, x = 'PC1', y = 'PC2',
       lab = NULL,
       colby = 'Disease', colkey = c('Control' = 'royalblue1', 'Patient' = 'coral'),
       colLegendTitle = 'Disease Status',
       encircle = TRUE,
       encircleFill = TRUE,
       hline = 0, vline = c(0, 0, 0),
       hlineWidth = 0.1, vlineWidth = 0.1,
       gridlines.major = FALSE, gridlines.minor = FALSE,
       legendPosition = 'right', legendLabSize = 16, legendIconSize = 8.0,
       title = "Pontine Nuclei")
dev.off()

tiff("PN4+6_correlates.tiff", units="in", width= 6, height=3.5, res=300)
eigencorplot(pca,
             components = getComponents(pca, 1:5),
             metavars = c('Age','PMI','Duration','LOT','Immunocompromised','Gender', 'Ventilation','Disease'),
             col = c('white', 'green3'),
             cexCorval = 1,
             fontCorval = 2,
             posLab = 'bottomleft', 
             rotLabX = 45,
             scale = TRUE,
             main = bquote(PN ~ PCA ~ Pearson ~ r^2 ~ clinical ~ correlates),
             cexMain = 1.5,
             plotRsquared = TRUE,
             corFUN = 'pearson',
             corUSE = 'pairwise.complete.obs',
             signifSymbols = c('***', '**', '*', ''),
             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 1))
dev.off()