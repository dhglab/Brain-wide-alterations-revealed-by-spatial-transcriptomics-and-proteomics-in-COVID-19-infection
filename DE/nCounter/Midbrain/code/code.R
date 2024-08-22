library(readxl)
library(dplyr)
library(ggplot2)
library(readxl)
library(EDASeq)
library(RColorBrewer)
library(RUVSeq)
library(Hmisc)
library(corrplot)
library(matrixStats)
library(openxlsx)
library(clusterProfiler)
library(WGCNA)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(GO.db)
library(org.Hs.eg.db)
library(sva)

set("./nCounter/Midbrain")
### INPUT: raw - p x n, and target
EXPR_Midbrain <- read.csv(file.choose())
target_Midbrain <- read.csv(file.choose())
colnames(target_Midbrain) <- c('Groups','Disease','Age','Gender',
                              'PMI','Stageposition','LaneNumber','Immunocompromised','BindingDensity')
rownames(target_Midbrain) <- target_Midbrain$Groups

# create raw counts dataframe Ex, gene annotation dataframe geneAnno
Ex_Midbrain <- EXPR_Midbrain[, -c(1:3)]; colnames(Ex_Midbrain) <- target_Midbrain$Groups
geneAnno_Midbrain <- EXPR_Midbrain[,c(1:3)]
rownames(Ex_Midbrain) <- rownames(geneAnno_Midbrain) <- geneAnno_Midbrain$Probe.Name

# re-order
target_Midbrain$Disease <- factor(target_Midbrain$Disease)
rownames(target_Midbrain) <- target_Midbrain$Groups

#sorting Disease
target_Midbrain <- target_Midbrain[order(target_Midbrain$Disease),]
Ex_Midbrain <- Ex_Midbrain[,rownames(target_Midbrain)]
all(rownames(target_Midbrain) == colnames(Ex_Midbrain)) #check rownames of target and colnames of matrix

## remove positive and negative control probels
idx <- grepl('ERCC', geneAnno_Midbrain$Accession)
Ex_Midbrain<- Ex_Midbrain[!idx,]
geneAnno_Midbrain <- geneAnno_Midbrain[!idx,]

#save rda
save(Ex_Midbrain, geneAnno_Midbrain, target_Midbrain, file= "Ex_Midbrain.rda")

#sample connectivity of RAW data
normadj <- (0.5+0.5*bicor(Ex_Midbrain))^2
netsummary= fundamentalNetworkConcepts(normadj)
connectivity=netsummary$Connectivity
connectivity.zscore = (connectivity-mean(connectivity))/sqrt(var(connectivity))
connectivity.plot <- data.frame(Sample.Name = names(connectivity.zscore), Z.score = connectivity.zscore, Sample.Num = 1:length(connectivity.zscore))

pdf("Midbrain_sample_connectivity.pdf")
p <- ggplot(connectivity.plot, aes(x=Sample.Num, y = Z.score, label = Sample.Name)) + geom_text(size = 4, colour = "red")
p <- p+ geom_hline(aes(yintercept = -2))
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle("Sample Connectivity")
print(p)
dev.off() 

#P4, P11 are outliers
#make new expression matrix set with pData target, fData geneAnno
set_Midbrain_P411<- newSeqExpressionSet(as.matrix(Ex_Midbrain_P411), phenoData = target_Midbrain_P411,  featureData = geneAnno_Midbrain)
set_Midbrain_P411

#Housekeeping normalization
housekeeping <- rownames(set_Midbrain_P411)[fData(set_Midbrain_P411)$Class.Name == 'Housekeeping']
set1_Midbrain_P411 <- RUVg(set_Midbrain_P411, housekeeping, k=1)

## save normalized counts(library size normalization+ RUVg normalization), fData and pData
normEx_Midbrain_P411 <- counts(set1_Midbrain_P411)
target_Midbrain_P411 <- pData(set1_Midbrain_P411)

save(normEx_Midbrain_P411, target_Midbrain_P411, geneAnno_Midbrain, set1_Midbrain_P411, file = "normEx_Midbrain_P411.rda")

#PCA
tiff("nCounterPCA_midbrain.tiff", units="in", width= 5, height=4, res=300)
plotPCA(set1_Midbrain_P411,col = c(rep("#005AB5", 12), rep("#DC3220",19)),cex=1.0,
        title = "Frontal lobe nCounter PCA")
dev.off()

#PCA analysis and correlation to Age, PMI, LaneNumber
Disease <- factor(set1_Midbrain_P411$Disease, levels = c("CONTROL", "PATIENT"))
pca = prcomp(t(normEx_Midbrain_P411[]))     
pcaData=as.data.frame(pca$x) 
prcomp_subset <- pcaData[,c(1:5)] 

#1.correlation to Age
Age <- target_Midbrain_P411$Age
COVB <- cbind(prcomp_subset, Age)
COVB<- COVB[sapply(COVB,is.numeric)]
res <- rcorr(as.matrix(COVB)) 

#2.correlation to PMI
PMI <- target_Midbrain_P411$PMI
COVB2 <- cbind(prcomp_subset, PMI)
COVB2<- COVB2[sapply(COVB2,is.numeric)]
res2 <- rcorr(as.matrix(COVB2)) 

#3. PCA analysis and correlation to LaneNumber
LaneNumber <- target_Midbrain_P411$LaneNumber
COVB3 <- cbind(prcomp_subset, LaneNumber)
COVB3<- COVB3[sapply(COVB3,is.numeric)]
res3 <- rcorr(as.matrix(COVB3)) 

corrplot(res$r, p.mat=res$P, sig.level = 0.05, insig="p-value", diag = FALSE, type = 'lower') 
corrplot(res2$r, p.mat=res2$P, sig.level = 0.05, insig="p-value", diag = FALSE, type = 'lower')   
corrplot(res3$r, p.mat=res3$P, sig.level = 0.05, insig="p-value", diag = FALSE, type = 'lower')

## DESeq2 differential analysis
dds_Midbrain_P411 <- DESeqDataSetFromMatrix(countData = normEx_Midbrain_P411,
                                            colData= target_Midbrain_P411,
                                            design = ~W_1+Disease)
dds_Midbrain_P411 <- DESeq(dds_Midbrain_P411)

rld <- rlog(dds_Midbrain_P411, blind = FALSE)
boxplot(assay(rld), las=2, main="rlog transformation", col=c(rep("Blue",12),rep("Red",19))) 

res_MidbrainP_vs_C <- results(dds_Midbrain_P411,contrast = c("Disease","PATIENT","CONTROL"))
res_Midbrain_P411 <- as.data.frame(res_MidbrainP_vs_C)
save(res_Midbrain_P411, file="res_MidbrainP_vs_C_P411.rda")
