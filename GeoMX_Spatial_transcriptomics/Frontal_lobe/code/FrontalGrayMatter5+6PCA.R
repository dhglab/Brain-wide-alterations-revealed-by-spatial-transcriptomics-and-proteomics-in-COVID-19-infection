library(dplyr)
library(ggplot2)
library(EDASeq)
library(RColorBrewer)
library(RUVSeq)
library(corrplot)
library(matrixStats)
library(WGCNA)
library(PCAtools)
library('DESeq2')
library(devtools)

#import GeoMX Frontal Gray Matter counts and DEG lists
FrontalGM <- as.data.frame(read_csv("C:/Users/tingz/OneDrive/USMLE/COVID Brain grant 2020/GeoMx DSP/Frontal GeoMX DSP/0_Preprocess_Frontal5+6/QC Grey Matter all WTA data_GroupName_removenonintegers.csv"))
rownames(FrontalGM) <- FrontalGM$TargetName
FrontalGM <- FrontalGM[, -c(1:1)]

DEG_GreyMatter_5_6 <- read_excel("C:/Users/tingz/Downloads/COVID_Revision/SupplementaryFig1Sourcedata/DEG_GreyMatter_5+6.xlsx")

#import metadata
TargetFrontalGM<- as.data.frame(read_csv("C:/Users/tingz/OneDrive/USMLE/COVID Brain grant 2020/GeoMx DSP/Frontal GeoMX DSP/0_Preprocess_Frontal5+6/5+6 Frontal Grey Target_GroupName_1.csv"))
rownames(TargetFrontalGM) <- TargetFrontalGM$Group
TargetFrontalGM <- TargetFrontalGM[, -c(1:1)]

#sorting Disease
TargetFrontalGM <- TargetFrontalGM[order(TargetFrontalGM$Disease),]

#Get FC > 0.3 and <- 0.3 and FDR < 0.1 for PCA plotting
FrontalGMFDR01 <- subset(DEG_GreyMatter_5_6, DEG_GreyMatter_5_6$`Adjusted pvalue` < 0.1)
FrontalGMFDR01_03 <- subset(FrontalGMFDR01, FrontalGMFDR01$Log2FC < -0.3)
FrontalGMFDR0103 <- subset(FrontalGMFDR01, FrontalGMFDR01$Log2FC > 0.3)
FrontalGMFDR0103 <- rbind(FrontalGMFDR0103,FrontalGMFDR01_03)
FrontalGMFDRQC <- subset(FrontalGMQC, rownames(FrontalGMQC) %in% FrontalGMFDR0103$`Target name`)

#check rownames of target and colnames of matrix
all(rownames(TargetFrontalGM) == colnames(FrontalGMFDRQC)) 

dds <- DESeqDataSetFromMatrix(countData = FrontalGMFDRQC,
                              colData = TargetFrontalGM,
                              design = ~ Disease)
dds

vsd <- vst(dds, blind=TRUE)
boxplot(assay(vsd), las=2, main="vsd")
plotPCA(vsd, intgroup = c("Disease"))

#Remove batch effect/LOT
mat <- assay(vsd)
mm <- model.matrix(~Disease, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$LOT, design=mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = c("Disease"))

pca <- pca(mat, metadata = TargetFrontalGM)

setwd("~/GitHub/Brain-wide-alterations-revealed-by-spatial-transcriptomics-and-proteomics-in-COVID-19-infection/GeoMX_Spatial_transcriptomics/Frontal_lobe/output")
tiff("FrontalGM5+6_PCA.tiff", units="in", width= 6.5, height=4, res=300)
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
       title = "Frontal Grey Matter")
dev.off()

tiff("FrontalGM5+6_correlates.tiff", units="in", width= 6, height=3.5, res=300)
eigencorplot(pca,
             components = getComponents(pca, 1:5),
             metavars = c('Age','PMI','Duration','LOT','Immunocompromised','Gender', 'Ventilation','Disease'),
             col = c('white', 'green3'),
             cexCorval = 1,
             fontCorval = 2,
             posLab = 'bottomleft', 
             rotLabX = 45,
             scale = TRUE,
             main = bquote(FrontalGM ~ PCA ~ Pearson ~ r^2 ~ clinical ~ correlates),
             cexMain = 1.5,
             plotRsquared = TRUE,
             corFUN = 'pearson',
             corUSE = 'pairwise.complete.obs',
             signifSymbols = c('***', '**', '*', ''),
             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 1))
dev.off()

setwd("C:/Users/tingz/OneDrive/Documents/GitHub/Brain-wide-alterations-revealed-by-spatial-transcriptomics-and-proteomics-in-COVID-19-infection/GeoMX_Spatial_transcriptomics/Frontal_lobe/data")
save(FrontalGMFDRQC, FrontalGMQC, DEG_GreyMatter_5_6,TargetFrontalGM, file = "FrontalGM5_6_PCA.rda")