library(readxl)
library(dplyr)
library(ggplot2)
library(readxl)
library(EDASeq)
library(DESeq2)
library(RColorBrewer)
library(RUVSeq)
library(Hmisc)
library(corrplot)
library(matrixStats)
library(openxlsx)
library(clusterProfiler)
library(WGCNA)
library(magrittr)
library(FactoMineR)
library(factoextra)
library(ggplot2)

set("./nCounter/Frontal")

#import rawcounts, meta_target data 
EXPR_Frontal <- read.csv(file.choose())
target_Frontal <- read.csv(file.choose())
colnames(target_Frontal) <- c('Groups','Disease','Age','Gender',
                              'PMI','Immunocompromised','Stageposition', 'LaneNumber','BindingDensity')
rownames(target_Frontal) <- target_Frontal$Groups

# create raw counts dataframe Ex, gene annotation dataframe geneAnno
geneAnno_Frontal <- as.data.frame(EXPR_Frontal[,c(1:3)])
rownames(geneAnno_Frontal) <- geneAnno_Frontal$`Probe Name`
Ex_FRONTAL <- as.data.frame(EXPR_Frontal[, -c(1:3)])
rownames(Ex_FRONTAL) <- EXPR_Frontal$`Probe Name`

# re-order
target_Frontal$Disease <- factor(target_Frontal$Disease)

#sorting Disease
target_Frontal <- target_Frontal[order(target_Frontal$Disease),]
Ex_FRONTAL <- Ex_FRONTAL[,rownames(target_Frontal)]
all(rownames(target_Frontal) == colnames(Ex_FRONTAL)) #check rownames of target and colnames of matrix

## remove positive and negative control probels
idx <- grepl('ERCC', geneAnno_Frontal$`Accession #`)
Ex_FRONTAL<- Ex_FRONTAL[!idx,]
geneAnno_Frontal <- geneAnno_Frontal[!idx,]
save(Ex_FRONTAL, geneAnno_Frontal, target_Frontal,file ="Ex_Frontal.rda")

###sample connectivity
normadj <- (0.5+0.5*bicor(Ex_FRONTAL))^2
netsummary= fundamentalNetworkConcepts(normadj)
connectivity=netsummary$Connectivity
connectivity.zscore = (connectivity-mean(connectivity))/sqrt(var(connectivity))
connectivity.plot <- data.frame(Sample.Name = names(connectivity.zscore), Z.score = connectivity.zscore, Sample.Num = 1:length(connectivity.zscore))

pdf("Frontal_sample_connectivity.pdf")
p <- ggplot(connectivity.plot, aes(x=Sample.Num, y = Z.score, label = Sample.Name)) + geom_text(size = 4, colour = "red")
p <- p+ geom_hline(aes(yintercept = -2))
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle("Sample Connectivity")
print(p)
dev.off()
#Outlier P2, P18 

#load Ex_FrontalMG218.rda (Ex_FRONTAL, target_Frontal, geneAnno_Frontal), outlier removed
#make new expression matrix set 
set_Frontal<- newSeqExpressionSet(as.matrix(Ex_FRONTAL), phenoData = target_Frontal,  featureData = geneAnno_Frontal)
set_Frontal

#RUVseq housekeeping normalization
housekeeping <- rownames(set_Frontal)[fData(set_Frontal)$Class.Name == 'Housekeeping']
set1_Frontal <- RUVg(set_Frontal, housekeeping, k=1)

## save normalized counts(library size normalization+ RUVg normalization), fData and pData
normEx_Frontal <- counts(set1_Frontal)
target_Frontal <- pData(set1_Frontal)

save(normEx_Frontal, target_Frontal, geneAnno_Frontal, set1_Frontal, file = "normEx_Frontal_MG218.rda")

#PCA_Supplementary Figure 3A
tiff("nCounterPCA_Frontal.tiff", units="in", width= 5, height=4, res=300)
plotPCA(set1_Frontal,col = c(rep("#005AB5", 12), rep("#DC3220",19)),cex=1.0,
        title = "Frontal lobe nCounter PCA")
dev.off()

#PCA analysis and correlation to Age, PMI, LaneNumber
Disease <- factor(set1_Frontal$Disease, levels = c("CONTROL", "PATIENT"))
pca = prcomp(t(normEx_Frontal[]))     
pcaData=as.data.frame(pca$x) 
prcomp_subset <- pcaData[,c(1:5)] 

#1.correlation to Age
Age <- target_Frontal$Age
COVB <- cbind(prcomp_subset, Age)
COVB<- COVB[sapply(COVB,is.numeric)]
res <- rcorr(as.matrix(COVB)) 

#2.correlation to PMI
PMI <- target_Frontal$PMI
COVB2 <- cbind(prcomp_subset, PMI)
COVB2<- COVB2[sapply(COVB2,is.numeric)]
res2 <- rcorr(as.matrix(COVB2)) 

#3. PCA analysis and correlation to LaneNumber
LaneNumber <- target_Frontal$LaneNumber
COVB3 <- cbind(prcomp_subset, LaneNumber)
COVB3<- COVB3[sapply(COVB3,is.numeric)]
res3 <- rcorr(as.matrix(COVB3)) 

corrplot(res$r, p.mat=res$P, sig.level = 0.05, insig="p-value", diag = FALSE, type = 'lower') 
corrplot(res2$r, p.mat=res2$P, sig.level = 0.05, insig="p-value", diag = FALSE, type = 'lower')   
corrplot(res3$r, p.mat=res3$P, sig.level = 0.05, insig="p-value", diag = FALSE, type = 'lower')

#do not correlate with Age, PMI or LaneNumber

## DESeq2 differential analysis
library(DESeq2)
dds_Frontal <- DESeqDataSetFromMatrix(countData = normEx_Frontal,
                              colData= target_Frontal,
                              design = ~ W_1 + Disease)
dds_Frontal <- DESeq(dds_Frontal)

rld <- rlog(dds_Frontal, blind = FALSE)
boxplot(assay(rld), las=2, main="rlog transformation")

res_Frontal <- results(dds_Frontal, contrast = c("Disease","PATIENT","CONTROL"))
res_Frontal<- as.data.frame(res_Frontal)
save(res_Frontal, file="res_Frontal_P_vs_C_MG218.rda")

#heapmap cluster of FDR genes
rld <- rlog(dds_Frontal, blind = FALSE)
rnormEx_Frontal<- assay(rld)
rnormEx_Frontal_scaled = t(apply(rnormEx_Frontal, 1, scale))
colnames(rnormEx_Frontal_scaled) = colnames(normEx_Frontal)

FDR01 <- subset(res_Frontal, res_Frontal$padj <0.1)
rnormEx_Frontal_scaledFDR <- subset(rnormEx_Frontal_scaled, subset = rownames(rnormEx_Frontal_scaled) %in% rownames(FDR01))

library(ComplexHeatmap)
library(circlize)
target_Frontal_h <- data.frame(Age = target_Frontal$Age,
                               Gender = target_Frontal$Gender,
                               PMI = target_Frontal$PMI,
                               Immunocompromised = target_Frontal$Immunocompromised,
                               LaneNumber = target_Frontal$LaneNumber) 
colAnn <- HeatmapAnnotation(df = target_Frontal_h, which = 'col')

tiff("Supplementary_FigS3C.tiff", width = 9.5, height = 6.5, res = 300, units = "in")
ht1 <- Heatmap(rnormEx_Frontal_scaledFDR, show_row_names = FALSE, show_column_names = TRUE, top_annotation = colAnn, col = colorRamp2(c(-1.5, 0, 1.5), c("#0072B2", "white", "#CC79A7")), column_names_gp = gpar(col = c(rep("black", 12), rep("#D55E00", 17))),cluster_columns = TRUE, raster_by_magick=FALSE, name = "Gene Expr")
draw (ht1)
dev.off()

#Supplementary Figure 3B
tiff("Frontal_Age_nCounter.tiff", height = 4, width = 4, res = 300, units = "in")
p <- ggplot(target_Frontal, aes(x = Disease, y = Age, fill= Disease,)) + geom_boxplot()+ geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1)) + theme_minimal()+ stat_compare_means( aes(label = ..p.signif..), label.x = 1.5, label.y = 90)
p <- p + labs(y = "Age (years)")
p <- p + theme(text = element_text(size = 14))
p
dev.off()

tiff("Frontal_PMI_nCounter.tiff", height = 4, width = 4, res = 300, units = "in")
p2 <- ggplot(target_Frontal, aes(x = Disease, y = PMI, fill= Disease,)) + geom_boxplot()+ geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1)) + theme_minimal() + stat_compare_means( aes(label = ..p.signif..), label.x = 1.5, label.y = 8) + scale_fill_manual(values= c("royalblue1","pink"))
p2 <- p2 + labs(y = "PMI (days)")
p2 <- p2 + theme(text = element_text(size = 14))   
p2
dev.off()

tiff("Frontal_Gender_nCounter.tiff", height = 4, width = 3, res = 300, units = "in")
pl <- ggplot(data = target_Frontal,aes(x= Disease, fill = Gender)) +
  geom_bar(stat="count", position ="fill") +
  labs(y = "Percentage") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Pastel1") +
  geom_bar(colour = "black", position = "fill") +
  theme_minimal() +
  theme(legend.position="top") +
  theme(text = element_text(size = 14))
pl 
dev.off()

tiff("Frontal_Immunocomprimised_nCounter.tiff", height = 4, width = 3, res = 300, units = "in")
colnames(target_Frontal)[6] <- "ImmuCom"
pl2 <- ggplot(data = target_Frontal,aes(x= Disease, fill = ImmuCom)) +
  geom_bar(stat="count", position ="fill") +
  labs(y = "Percentage") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "jco") +
  geom_bar(colour = "black", position = "fill") +
  theme_minimal() +
  theme(legend.position="top") +
  theme(text = element_text(size = 14))
pl2
dev.off()

#Reviewer 1 question: in control groups, non_lung pathology vs lung pathology, load target_Frontal_control with lung pathology, and counts Ex_Frontal_MG218.rda
target_Frontal_control <- as.data.frame(target_control_trim)
rownames(target_Frontal_control) <- target_Frontal_control$Groups

#subset control out
Ex_FRONTAL_control <- Ex_FRONTAL[,1:12]

#make new expression matrix set with pData target, fData geneAnno
set_Frontal<- newSeqExpressionSet(as.matrix(Ex_FRONTAL_control), phenoData = target_Frontal_control,  featureData = geneAnno_Frontal)
set_Frontal

housekeeping <- rownames(set_Frontal)[fData(set_Frontal)$Class.Name == 'Housekeeping']
set1_Frontal <- RUVg(set_Frontal, housekeeping, k=1)

normEx_Frontal_control <- counts(set1_Frontal)
target_Frontal_control <- pData(set1_Frontal)

save(normEx_Frontal_control, target_Frontal_control, geneAnno_Frontal, set1_Frontal, geneAnno_Frontal, file ='normEx_Frontal_control.rda')

dds_Frontal <- DESeqDataSetFromMatrix(countData = normEx_Frontal_control,
                                      colData= target_Frontal_control,
                                      design = ~ W_1 + LungPath)
dds_Frontal <- DESeq(dds_Frontal)

rld <- rlog(dds_Frontal, blind = FALSE)
pdf("Rlog_Frontal_control.pdf", width = 6, height = 6)
boxplot(assay(rld), las=2, main="rlog transformation", col=c(rep("#005AB5",12)))
dev.off ()

res_Y_vs_N <- results(dds_Frontal,contrast = c("LungPath","Y","N"))
res_Y_vs_N <- as.data.frame(res_Y_vs_N)
save(res_Y_vs_N, file = "res_Y_vs_N_lungpath.rda")
#No DEG Genes, the results are Supplementary data 4. 

