library(dplyr)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(Hmisc)
library(corrplot)
library(matrixStats)
library(WGCNA)
library(Biobase)
library(PCAtools)
library('DESeq2')
library(devtools)
library(parallel)
library(abind)
library(progress)
library(lme4)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

setwd("./Fig2Sourcedata/WGCNA")

#import GeoMX Frontal Gray Matter normalized counts, Metadata and sequencing data from FrontalGMQ3_Meta_before.rda
#load FrontalGMQ3, TargetFrontalGM, FrontalGMMetaSeq 

#scale Age, PMI, Duration to numeric form#
Age_scaled <- scale(as.numeric(TargetFrontalGM$Age))
TargetFrontalGM$Age = as.numeric(paste(Age_scaled))

PMI_scaled <- scale(as.numeric(TargetFrontalGM$PMI))
TargetFrontalGM$PMI = as.numeric(paste(PMI_scaled))

Duration_scaled <- scale(as.numeric(TargetFrontalGM$Duration))
TargetFrontalGM$Duration = as.numeric(paste(Duration_scaled))

colnames(TargetFrontalGM)[17] = c("LOT")
colnames(TargetFrontalGM)[18] = c("ScanID")

TargetFrontalGM$Gender[TargetFrontalGM$Gender == 'Male'] <- 'M'
TargetFrontalGM$Gender[TargetFrontalGM$Gender == 'Female'] <- 'F'

TargetFrontalGM$Ventilation[TargetFrontalGM$Ventilation == 'Yes'] <- 'Y'
TargetFrontalGM$Ventilation[TargetFrontalGM$Ventilation == 'No'] <- 'N'

TargetFrontalGM$Immunocompromised[TargetFrontalGM$Immunocompromised == 'Yes'] <- 'Y'
TargetFrontalGM$Immunocompromised[TargetFrontalGM$Immunocompromised == 'No'] <- 'N'

##Import Metaseqdata
TargetFrontalGM_MetaSeq<- as.data.frame(FrontalGMMetaSeq)
rownames(TargetFrontalGM_MetaSeq) <-TargetFrontalGM_MetaSeq$Group
TargetFrontalGM_Seq <- TargetFrontalGM_MetaSeq[, -c(1:9)]
TargetFrontalGM_Seq <-scale(TargetFrontalGM_Seq)

###sample connectivity
normadj <- (0.5+0.5*bicor(FrontalGMQ3))^2
netsummary= fundamentalNetworkConcepts(normadj)
connectivity=netsummary$Connectivity
connectivity.zscore = (connectivity-mean(connectivity))/sqrt(var(connectivity))
connectivity.plot <- data.frame(Sample.Name = names(connectivity.zscore), Z.score = connectivity.zscore, Sample.Num = 1:length(connectivity.zscore))

pdf("GMQ35+6Connectivity.pdf", height = 10, width = 10)
p <- ggplot(connectivity.plot, aes(x=Sample.Num, y = Z.score, label = Sample.Name)) + geom_text(size = 4, colour = "red")
p <- p+ geom_hline(aes(yintercept = -2))
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle("Sample Connectivity")
print(p)
dev.off()

outliers <- (connectivity.zscore > mean(connectivity.zscore)+2*sd(connectivity.zscore))|(connectivity.zscore < mean(connectivity.zscore)-2*sd(connectivity.zscore))
print(colnames(FrontalGMQ3)[outliers])
#outliers P5_2, P5_4

FrontalGMQ3_after <- FrontalGMQ3[,!names(FrontalGMQ3) %in% c("P5_2","P5_4")]
TargetFrontalGM_after <- TargetFrontalGM[!(row.names(TargetFrontalGM) %in% c("P5_2","P5_4")),]

#linear regression#
#use get PC function for sequencing PC identification#
getTopPCs <- function(x, n, col_name="") {
  scaled_dat <- x
  prcomp_dat <- prcomp(scaled_dat,center=F);
  topPCs <- prcomp_dat$rotation[,1:n];
  varexp <- (prcomp_dat$sdev)^2 / sum(prcomp_dat$sdev^2)
  topvar <- varexp[1:n]
  colnames(topPCs) <- paste0(col_name, "\n", colnames(topPCs)," (",signif(100*topvar[1:n],2),"%)")
  return(topPCs)
}

Q3Expr_t <- t(FrontalGMQ3_after)
c <- rownames(TargetFrontalGM_after)
Target_Seq <- subset(TargetFrontalGM_Seq, subset = rownames(TargetFrontalGM_Seq) %in% c)

seqPCs <- data.frame(getTopPCs(t(scale(Target_Seq)), 5))
colnames(seqPCs) <- paste0("seqpc", 1:5)
seqnames <- colnames(seqPCs)

target <- cbind(TargetFrontalGM_after, seqPCs)
target_Frontal <- cbind(target, Target_Seq)

target_Frontal$Gender <- as.factor(target_Frontal$Gender) 
target_Frontal$Gender <- as.numeric(target_Frontal$Gender)

target_Frontal$Immunocompromised <- as.factor(target_Frontal$Immunocompromised)
target_Frontal$Immunocompromised <- as.numeric(target_Frontal$Immunocompromised) 

target_Frontal$Ventilation <- as.factor(target_Frontal$Ventilation) 
target_Frontal$Ventilation <- as.numeric(target_Frontal$Ventilation) 

Q3Expr_t <- t(FrontalGMQ3_after)
lmFrontal <- cbind(target_Frontal, Q3Expr_t)

formula <- as.formula(paste0("Q3Expr_t ~ Disease + ScanID + LOT + Age + PMI + Ventilation + Gender + Duration + Immunocompromised + ", paste0("seqpc", 1:5, collapse = " + ")))
lmu <- lm(formula, data = lmFrontal)
Coefficentlist <- coef(summary(lmu))
Coefficentlist[[1]]

# Correct for sequencing effects. -------- 
# Stop + remake cluster; regression is single-core, so we want all cores to be in use.
cl <- makeCluster(detectCores())
# Convert Batch to a binary matrix and export to cluster.
Gender_binary <- model.matrix(~0 + target_Frontal$Gender)
Immunocompromised_binary <- model.matrix(~0 + target_Frontal$Immunocompromised)
Ventilation_binary <- model.matrix(~0 + target_Frontal$Ventilation)
ScanID_binary <- model.matrix(~0 + target_Frontal$ScanID)
clusterExport(cl=cl, list('Q3Expr_t', 'target_Frontal', 'Coefficentlist', 'ScanID_binary','seqnames','seqPCs'))
clusterEvalQ(cl=cl, gc())
crange <- 1:nrow(FrontalGMQ3_after)

uncor <- Q3Expr_t
age_effect <- parSapply(cl, crange, function(i){age_effect <- as.matrix(target_Frontal[, "Age"]) %*% as.matrix(Coefficentlist[[i]]["Age", "Estimate"])})

scan_effect <- parSapply(cl, crange, function(i){scan_effect <- as.matrix(ScanID_binary[, -1]) %*% as.matrix(Coefficentlist[[i]][c("ScanID2_Frontal", "ScanID3_Frontal", "ScanID4_Frontal", "ScanIDFrontal","ScanIDPons"), "Estimate"])})

Gendereffect <- parSapply(cl, crange, function(i){gender_effect <-  as.numeric(target_Frontal[, "Gender"]) %*% as.matrix(Coefficentlist[[i]]["Gender", "Estimate"])})

PMI_effect <- parSapply(cl, crange, function(i){pmi_effect <- as.matrix(target_Frontal[, "PMI"]) %*% as.matrix(Coefficentlist[[i]]["PMI", "Estimate"])})

Ventilation_effect <- parSapply(cl, crange, function(i){Vent_effect <- as.numeric(target_Frontal[, "Ventilation"]) %*% as.matrix(Coefficentlist[[i]]["Ventilation", "Estimate"])})

seqeffect <- parSapply(cl, crange, function(i){seq_effect <- as.matrix(target_Frontal[,seqnames]) %*% as.matrix(Coefficentlist[[i]][seqnames, "Estimate"])})

coef_corr_mat <- (uncor - age_effect - scan_effect - Gendereffect - PMI_effect - Ventilation_effect - seqeffect)
FrontalGMQ3_reg <- t(coef_corr_mat)

stopCluster(cl)

save(FrontalGMQ3_reg, FrontalGMQ3_after, TargetFrontalGM_after, coef_corr_mat, lmFrontal, lmu, target_Frontal, TargetFrontalGM_Seq, Target_Seq, Coefficentlist, file = "FrontalGM_linearreg_for_WGCNA.RData")

#Continue with WGCNA###
#2.c.1 Choosing the soft-thresholding power: analysis of network topology
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(t(FrontalGMQ3_reg), powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
pdf("Power_connectivity_FrontalGM.pdf", height = 8, width = 8)
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#Stepwise constructions#
softPower = 12;

#2.b.3 Topological Overlap Matrix (TOM)
TOM=TOMsimilarityFromExpr(t(FrontalGMQ3_reg), power = 12, 
                          networkType = "signed", 
                          TOMType = "unsigned",
                          TOMDenom = "min",
                          verbose = 10)
dissTOM = 1-TOM

#2.b.4 Clustering using TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(15,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

save(TOM, file = "FrontalGM_TOM.RDATA")
