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

#Load raw(QC_remove non integer) and Q3 normalized counts, targets, and differential expression from GeoDSP profiler, FrontalWMQC, WMDEG,TargetFrontalWM

#check rownames of target and colnames of matrix
all(rownames(TargetFrontalWM) == colnames(FrontalWMQC)) 

WMFDR01 <- subset(WMDEG, WMDEG$`Adjusted pvalue` < 0.1)
WMFDR01QC <- subset(FrontalWMQC, rownames(FrontalWMQC) %in% WMFDR01$`Target name`)

dds <- DESeqDataSetFromMatrix(countData = WMFDR01QC,
                              colData = TargetFrontalWM,
                              design = ~ Disease)
dds

vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
boxplot(assay(vsd), las=2, main="vsd")
plotPCA(vsd, intgroup = c("Disease"))

#remove batch/LOT effect
mat <- assay(vsd)
mm <- model.matrix(~Disease, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$LOT, design=mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = c("Disease"))

pca <- pca(mat, metadata = TargetFrontalWM)

#Supplementary Figure 1A and 1B
tiff("FrontalWM5+6_PCA.tiff", units="in", width= 6.5, height=4, res=300)
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
       title = "Frontal White Matter")
dev.off()

tiff("FrontalWM5+6_correlate.tiff", units="in", width= 6, height=3.5, res=300)
eigencorplot(pca,
             components = getComponents(pca, 1:5),
             metavars = c('Age','PMI','Duration','LOT','Immunocompromised','Gender', 'Ventilation','Disease'),
             col = c('white', 'green3'),
             cexCorval = 1,
             fontCorval = 2,
             posLab = 'bottomleft', 
             rotLabX = 45,
             scale = TRUE,
             main = bquote(FrontalWM ~ PCA ~ Pearson ~ r^2 ~ clinical ~ correlates),
             cexMain = 1.5,
             plotRsquared = TRUE,
             corFUN = 'pearson',
             corUSE = 'pairwise.complete.obs',
             signifSymbols = c('***', '**', '*', ''),
             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 1))
dev.off()
