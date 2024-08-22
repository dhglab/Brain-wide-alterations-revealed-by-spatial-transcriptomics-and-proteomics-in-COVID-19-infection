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
library(sva)

set("./nCounter/BasalGanglia")
### INPUT: raw - p x n raw expressions with p genes and n samples
EXPR_BG1 <- read.csv(file.choose())
EXPR_BG2 <- read.csv(file.choose())

# get target.csv
#target <- read.csv(paste0(dir,'/target.csv'))
target_BG1 <- read.csv(file.choose())
colnames(target_BG1) <- c('Sample','File','Disease','Groups','Batch','Age','Gender','PMI','Immunocompromised','Stageposition','LaneNumber','BindingDensity')
rownames(target_BG1) <- target_BG1$Groups

target_BG2 <- read.csv(file.choose())
colnames(target_BG2) <- c('Sample','File','Disease','Groups','Batch','Age','Gender','PMI','Immunocompromised','Stageposition','LaneNumber','BindingDensity')
rownames(target_BG2) <- target_BG2$Groups

# create raw counts dataframe Ex, gene annotation dataframe geneAnno
Ex_BG1 <- EXPR_BG1[, -c(1:3)]; colnames(Ex_BG1) <- target_BG1$Groups
geneAnno_BG1<- EXPR_BG1[,c(1:3)]; colnames(ProbeNames_1)
rownames(Ex_BG1) <- rownames(geneAnno_BG1) <- geneAnno_BG1$Probe.Name

Ex_BG2 <- EXPR_BG2[, -c(1:3)]; colnames(Ex_BG2) <- target_BG2$Groups
geneAnno_BG2<- EXPR_BG2[,c(1:3)]; colnames(ProbeNames_2)
rownames(Ex_BG2) <- rownames(geneAnno_BG2) <- geneAnno_BG2$Probe.Name

# re-order
target_BG1$Disease <- factor(target_BG1$Disease)
target_BG2$Disease <- factor(target_BG2$Disease)

#sorting Disease
target_BG1 <- target_BG1[order(target_BG1$Disease),]
Ex_BG1 <- Ex_BG1[,rownames(target_BG1)]
all(rownames(target_BG1) == colnames(Ex_BG1))

target_BG2 <- target_BG2[order(target_BG2$Disease),]
Ex_BG2 <- Ex_BG2[,rownames(target_BG2)]
all(rownames(target_BG2) == colnames(Ex_BG2))

## remove positive and negative control probels
idx <- grepl('ERCC', geneAnno_BG1$Accession)
Ex_BG1<- Ex_BG1[!idx,]
geneAnno_BG1 <- geneAnno_BG1[!idx,]

idx <- grepl('ERCC', geneAnno_BG2$Accession)
Ex_BG2<- Ex_BG2[!idx,]
geneAnno_BG2 <- geneAnno_BG2[!idx,]

#save rda
save(Ex_BG1, geneAnno_BG1, target_BG1, ProbeNames_1, file= "Ex_BG1.rda")
save(Ex_BG2, geneAnno_BG2, target_BG2, ProbeNames_2, file= "Ex_BG2.rda")

#combine Ex_BG1, Ex_BG2
Sharedgenes <- intersect(rownames(Ex_BG1),rownames(Ex_BG2))
Ex_BG1s<- subset(Ex_BG1, subset = rownames(Ex_BG1) %in% Sharedgenes)
Ex_BG2s<- subset(Ex_BG2, subset = rownames(Ex_BG2) %in% Sharedgenes)
geneAnno_BG <- subset(geneAnno_BG2, subset = rownames(geneAnno_BG2) %in% Sharedgenes)

#input combined EXPR_BG and target_BG
EXPR_BG <- read.csv(file.choose())
target_BG <- read.csv(file.choose())
colnames(target_BG) <- c('Sample','File','Disease','Groups','Batch','Age','Gender','PMI','Immunocompromised','Stageposition','LaneNumber','BindingDensity')
rownames(target_BG) <- target_BG$Groups

# create raw counts dataframe Ex, gene annotation dataframe geneAnno
Ex_BG <- EXPR_BG[, -c(1:1)]
rownames(Ex_BG) <- geneAnno_BG$Probe.Name 

#check rownames of target and colnames of matrix
all(rownames(target_BG) == colnames(Ex_BG)) 

#save rda
save(Ex_BG, geneAnno_BG, target_BG, file="Ex_BG.rda")

normadj <- (0.5+0.5*bicor(Ex_BG))^2
netsummary= fundamentalNetworkConcepts(normadj)
connectivity=netsummary$Connectivity
connectivity.zscore = (connectivity-mean(connectivity))/sqrt(var(connectivity))
connectivity.plot <- data.frame(Sample.Name = names(connectivity.zscore), Z.score = connectivity.zscore, Sample.Num = 1:length(connectivity.zscore))

pdf("BG_combined_connectivity.pdf")
p <- ggplot(connectivity.plot, aes(x=Sample.Num, y = Z.score, label = Sample.Name)) + geom_text(size = 4, colour = "red")
p <- p+ geom_hline(aes(yintercept = -2))
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle("Sample Connectivity")
print(p)
dev.off() 
#P10, P12 are outliers

#load Ex_BG_P1012.rda, outliers removed
set_BG<- newSeqExpressionSet(as.matrix(Ex_BG_P1012), phenoData = target_BG_P1012,  featureData = geneAnno_BG)
set_BG

#Housekeeping normalization
housekeeping <- rownames(set_BG)[fData(set_BG)$Class.Name == 'Housekeeping']
set1_BG <- RUVg(set_BG, housekeeping, k=1)

normEx_BG <- counts(set1_BG)
target_BG <- pData(set1_BG)

save(normEx_BG, target_BG, geneAnno_BG, set1_BG, file = "normEx_BG_P1012.rda")

#PCA
tiff("nCounterPCA_BG.tiff", units="in", width= 5, height=4, res=300)
plotPCA(set1_BG,col = c(rep("#005AB5", 11), rep("#DC3220",17)),cex=1.0,
        title = "BG nCounter PCA")
dev.off()

#PCA analysis and correlation to Age, PMI, batch
Disease <- factor(set1_BG$Disease, levels = c("CONTROL", "PATIENT"))
pca = prcomp(t(normEx_BG[]))     
pcaData=as.data.frame(pca$x) 
prcomp_subset <- pcaData[,c(1:5)] 

#1.correlation to Age
Age <- target_BG$Age
COVB <- cbind(prcomp_subset, Age)
COVB<- COVB[sapply(COVB,is.numeric)]
res <- rcorr(as.matrix(COVB)) 

#2.correlation to PMI
PMI <- target_BG$PMI
COVB2 <- cbind(prcomp_subset, PMI)
COVB2<- COVB2[sapply(COVB2,is.numeric)]
res2 <- rcorr(as.matrix(COVB2)) 

#3. PCA analysis and correlation to Batch
Batch <- target_BG$Batch
COVB3 <- cbind(prcomp_subset, Batch)
COVB3<- COVB3[sapply(COVB3,is.numeric)]
res3 <- rcorr(as.matrix(COVB3)) 

corrplot(res$r, p.mat=res$P, sig.level = 0.05, insig="p-value", diag = FALSE, type = 'lower') 
corrplot(res2$r, p.mat=res2$P, sig.level = 0.05, insig="p-value", diag = FALSE, type = 'lower')   
corrplot(res3$r, p.mat=res3$P, sig.level = 0.05, insig="p-value", diag = FALSE, type = 'lower')
#correlate to Batch

## #DESeq2 SVA remove Batch effects
dds_BG <- DESeqDataSetFromMatrix(countData = normEx_BG,
                                 colData= target_BG,
                                 design = ~ W_1 + Disease)

dat  <- counts(dds_BG)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ Disease, colData(dds_BG))
mod0 <- model.matrix(~   1, colData(dds_BG))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)
svseq$sv

for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ dds_BG$Batch, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}

ddssva_BG <- dds_BG
ddssva_BG$SV1 <- svseq$sv[,1]
ddssva_BG$SV2 <- svseq$sv[,2]
design(ddssva_BG) <- ~ SV1 + SV2 + Disease

rld <- rlog(ddssva_BG, blind = FALSE)
boxplot(assay(rld), las=2, main="rlog transformation", col=c(rep("Blue",11),rep("Red",17)))

ddssva_BG <- DESeq(ddssva_BG)
res_sva_BG <- as.data.frame(results(ddssva_BG, contrast = c("Disease","PATIENT","CONTROL")))
save(res_sva_BG, file = "res_BG_P_vs_C_sva.rda")
