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
library(clusterProfiler)
library(WGCNA)

set("./nCounter/Pons")
### INPUT: raw - p x n and targets
EXPR_Pons <- read.csv(file.choose())
target_Pons <- read.csv(file.choose())
colnames(target_Pons) <- c('Sample','Groups','Disease','Age','Gender','PMI','Stageposition','LaneNumber','Immunocomprimised')
rownames(target_Pons) <- target_Pons$Groups

# create raw counts dataframe Ex, gene annotation dataframe geneAnno
Ex_Pons <- EXPR_Pons[, -c(1:3)]; colnames(Ex_Pons) <- target_Pons$Groups
geneAnno_Pons <- EXPR_Pons[,c(1:3)]
rownames(Ex_Pons) <- rownames(geneAnno_Pons) <- geneAnno_Pons$Probe.Name

#sorting Disease
target_Pons <- target_Pons[order(target_Pons$Disease),]
Ex_Pons <- Ex_Pons[,rownames(target_Pons)]
all(rownames(target_Pons) == colnames(Ex_Pons)) #check rownames of target and colnames of matrix

## remove positive and negative control probels
idx <- grepl('ERCC', geneAnno_Pons$Accession)
Ex_Pons<- Ex_Pons[!idx,]
geneAnno_Pons <- geneAnno_Pons[!idx,]

save(Ex_Pons, geneAnno_Pons, target_Pons, file= "Ex_Pons.rda")

###sample connectivity
normadj <- (0.5+0.5*bicor(Ex_Pons))^2
netsummary= fundamentalNetworkConcepts(normadj)
connectivity=netsummary$Connectivity
connectivity.zscore = (connectivity-mean(connectivity))/sqrt(var(connectivity))
connectivity.plot <- data.frame(Sample.Name = names(connectivity.zscore), Z.score = connectivity.zscore, Sample.Num = 1:length(connectivity.zscore))

pdf(paste0(outdir, "/", "Pons_sample_connectivity.pdf"), height = 10, width = 10)
p <- ggplot(connectivity.plot, aes(x=Sample.Num, y = Z.score, label = Sample.Name)) + geom_text(size = 4, colour = "red")
p <- p+ geom_hline(aes(yintercept = -2))
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle("Sample Connectivity")
print(p)
dev.off()     
#P2 are outliers

#load Ex_Pons_2.rda, outlier removed
#make new expression matrix set with pData target, fData geneAnno
set_Pons_2<- newSeqExpressionSet(as.matrix(Ex_Pons_2), phenoData = target_Pons_2,  featureData = geneAnno_Pons_2)
set_Pons_2

#RUVseq housekeeping normalization
##getting housekeeping genes from set
housekeeping <- rownames(set_Pons_2)[fData(set_Pons_2)$Class.Name == 'Housekeeping']
set1_Pons_2 <- RUVg(set_Pons_2, housekeeping, k=1)

## save normalized counts
normEx_Pons_2 <- counts(set1_Pons_2)
target_Pons_2 <- pData(set1_Pons_2)

save(normEx_Pons_2, target_Pons_2, geneAnno_Pons_2, set1_Pons_2, file = "normEx_Pons_2.rda")

#PCA
tiff("nCounterPCA_Pons.tiff", units="in", width= 5, height=4, res=300)
plotPCA(set1_Pons_2,col = c(rep("#005AB5", 10), rep("#DC3220",16)),cex=1.0,
        title = "Frontal lobe nCounter PCA")
dev.off()

#PCA analysis and correlation to Age, PMI, batch
Disease <- factor(set1_Pons_2$Disease, levels = c("CONTROL", "PATIENT"))
pca = prcomp(t(normEx_Pons_2[]))     
pcaData=as.data.frame(pca$x) 
prcomp_subset <- pcaData[,c(1:5)] 

#1.correlation to Age
Age <- target_Pons_2$Age
COVB <- cbind(prcomp_subset, Age)
COVB<- COVB[sapply(COVB,is.numeric)]
res <- rcorr(as.matrix(COVB)) 

#2.correlation to PMI
PMI <- target_Pons_2$PMI
COVB2 <- cbind(prcomp_subset, PMI)
COVB2<- COVB2[sapply(COVB2,is.numeric)]
res2 <- rcorr(as.matrix(COVB2)) 

#3. PCA analysis and correlation to Batch
LaneNumber <- target_Pons_2$LaneNumber
COVB3 <- cbind(prcomp_subset, LaneNumber)
COVB3<- COVB3[sapply(COVB3,is.numeric)]
res3 <- rcorr(as.matrix(COVB3)) 

corrplot(res$r, p.mat=res$P, sig.level = 0.05, insig="p-value", diag = FALSE, type = 'lower') 
corrplot(res2$r, p.mat=res2$P, sig.level = 0.05, insig="p-value", diag = FALSE, type = 'lower')   
corrplot(res3$r, p.mat=res3$P, sig.level = 0.05, insig="p-value", diag = FALSE, type = 'lower')

#Correlate to Age

#adjust Age effect
dds_Pons_2 <- DESeqDataSetFromMatrix(countData = normEx_Pons_2,
                                 colData= target_Pons_2,
                                 design = ~W_1 + Disease)

dat  <- counts(dds_Pons_2)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ Disease, colData(dds_Pons_2))
mod0 <- model.matrix(~   1, colData(dds_Pons_2))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)
svseq$sv

for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ dds_Pons_2$Age, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}

ddssva_Pons_2 <- dds_Pons_2
ddssva_Pons_2$SV1 <- svseq$sv[,1]
ddssva_Pons_2$SV2 <- svseq$sv[,2]
design(ddssva_Pons_2) <- ~ SV1 + SV2 + Disease

## DESeq2 differential analysis
dds_Pons_2 <- DESeq(dds_Pons_2)

rld <- rlog(dds_Pons_2, blind = FALSE)
boxplot(assay(rld), las=2, main="rlog transformation", col=c(rep("Blue",10),rep("Red",16))) 

res_Pons_2 <- results(dds_Pons_2, contrast = c("Disease","PATIENT","CONTROL"))
res_Pons_2 <- as.data.frame(res_Pons_2)
save(res_Pons, file = "res_Pons_P_vs_C_2.rda")

