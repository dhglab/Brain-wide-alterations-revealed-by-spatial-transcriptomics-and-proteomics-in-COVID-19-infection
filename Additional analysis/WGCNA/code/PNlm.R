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

#import GeoMX normalized counts, Meta, sequencing data from PN_Meta_seq_before.rda
#load PNQ3_Groupname, TargetPN, TargetPN_MetaSeq

#scale the Age, PMI, Duration to numeric form#
Age_scaled <- scale(as.numeric(TargetPN$Age))
TargetPN$Age = as.numeric(paste(Age_scaled))

PMI_scaled <- scale(as.numeric(TargetPN$PMI))
TargetPN$PMI = as.numeric(paste(PMI_scaled))

Duration_scaled <- scale(as.numeric(TargetPN$Duration))
TargetPN$Duration = as.numeric(paste(Duration_scaled))

TargetPN$Gender[TargetPN$Gender == 'Male'] <- 'M'
TargetPN$Gender[TargetPN$Gender == 'Female'] <- 'F'

TargetPN$Ventilation[TargetPN$Ventilation == 'Yes'] <- 'Y'
TargetPN$Ventilation[TargetPN$Ventilation == 'No'] <- 'N'

TargetPN$Immunocompromised[TargetPN$Immunocompromised == 'Yes'] <- 'Y'
TargetPN$Immunocompromised[TargetPN$Immunocompromised == 'No'] <- 'N'

#scale Metaseqdata
TargetPN_MetaSeq  <-scale(TargetPN_MetaSeq)

###sample connectivity
normadj <- (0.5+0.5*bicor(PNQ3_Groupname))^2
netsummary= fundamentalNetworkConcepts(normadj)
connectivity=netsummary$Connectivity
connectivity.zscore = (connectivity-mean(connectivity))/sqrt(var(connectivity))
connectivity.plot <- data.frame(Sample.Name = names(connectivity.zscore), Z.score = connectivity.zscore, Sample.Num = 1:length(connectivity.zscore))

pdf("PNQ34+6Connectivity.pdf", height = 10, width = 10)
p <- ggplot(connectivity.plot, aes(x=Sample.Num, y = Z.score, label = Sample.Name)) + geom_text(size = 4, colour = "red")
p <- p+ geom_hline(aes(yintercept = -2))
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle("Sample Connectivity")
print(p)
dev.off()

outliers <- (connectivity.zscore > mean(connectivity.zscore)+2*sd(connectivity.zscore))|(connectivity.zscore < mean(connectivity.zscore)-2*sd(connectivity.zscore))
print(colnames(PNQ3_Groupname)[outliers])
#outliers "P6_4" "P6_5"

PNQ3_after <- PNQ3_Groupname[,!names(PNQ3_Groupname) %in% c("P6_4", "P6_5")]
TargetPN_after <- TargetPN[!(row.names(TargetPN) %in% c("P6_4", "P6_5")),]

#linear regression#
getTopPCs <- function(x, n, col_name="") {
  scaled_dat <- x
  prcomp_dat <- prcomp(scaled_dat,center=F);
  topPCs <- prcomp_dat$rotation[,1:n];
  varexp <- (prcomp_dat$sdev)^2 / sum(prcomp_dat$sdev^2)
  topvar <- varexp[1:n]
  colnames(topPCs) <- paste0(col_name, "\n", colnames(topPCs)," (",signif(100*topvar[1:n],2),"%)")
  return(topPCs)
}

c <- rownames(TargetPN_after)
Target_Seq <- subset(TargetPN_MetaSeq, subset = rownames(TargetPN_MetaSeq) %in% c)

seqPCs <- data.frame(getTopPCs(t(scale(Target_Seq)), 5))
colnames(seqPCs) <- paste0("seqpc", 1:5)
seqnames <- colnames(seqPCs)

target <- cbind(TargetPN_after, seqPCs)
target_PN <- cbind(target, Target_Seq)

target_PN$Gender <- as.factor(target_PN$Gender) 
target_PN$Gender <- as.numeric(target_PN$Gender)

target_PN$Immunocompromised <- as.factor(target_PN$Immunocompromised)
target_PN$Immunocompromised <- as.numeric(target_PN$Immunocompromised) 

target_PN$Ventilation <- as.factor(target_PN$Ventilation) 
target_PN$Ventilation <- as.numeric(target_PN$Ventilation) 


Q3Expr_t<- t(PNQ3_after)
lmPN <- cbind(target_PN, Q3Expr_t)

formula <- as.formula(paste0("Q3Expr_t ~ Disease + ScanID + LOT + Age + PMI + 
                            Ventilation + Gender + Duration + Immunocompromised + ", 
                             paste0("seqpc", 1:5, collapse = " + ")))
lm <- lm(formula, data = lmPN)
Coefficentlist <- coef(summary(lm))
Coefficentlist[[1]]

# Correct for sequencing effects. -------- 
# Stop + remake cluster; regression is single-core, so we want all cores to be in use.
cl <- makeCluster(detectCores())
# Convert Batch to a binary matrix and export to cluster.
ScanID_binary <- model.matrix(~0 + target_PN$ScanID)
clusterExport(cl=cl, list('Q3Expr_t', 'target_PN', 'Coefficentlist', 'ScanID_binary','seqnames','seqPCs'))
clusterEvalQ(cl=cl, gc())
crange <- 1:nrow(PNQ3_after)

uncor <- Q3Expr_t
age_effect <- parSapply(cl, crange, function(i){age_effect <- as.matrix(target_PN[, "Age"]) %*% as.matrix(Coefficentlist[[i]]["Age", "Estimate"])})

scan_effect <- parSapply(cl, crange, function(i){scan_effect <- as.matrix(ScanID_binary[, -1]) %*% as.matrix(Coefficentlist[[i]][c("ScanID6_Pons", "ScanID7_Pons", "ScanID7_Pons-Control", "ScanID9_MidBrain_Pons-COVID","ScanIDPons"), "Estimate"])})

PMI_effect <- parSapply(cl, crange, function(i){pmi_effect <- as.matrix(target_PN[, "PMI"]) %*% as.matrix(Coefficentlist[[i]]["PMI", "Estimate"])})

Ventilation_effect <- parSapply(cl, crange, function(i){Vent_effect <- as.numeric(target_PN[, "Ventilation"]) %*% as.matrix(Coefficentlist[[i]]["Ventilation", "Estimate"])})

seqeffect <- parSapply(cl, crange, function(i){seq_effect <- as.matrix(target_PN[,seqnames]) %*% as.matrix(Coefficentlist[[i]][seqnames, "Estimate"])})

coef_corr_mat <- (uncor - age_effect - scan_effect - PMI_effect - Ventilation_effect - seqeffect)
PNQ3_reg <- coef_corr_mat

save(PNQ3_reg, seqPCs, seqnames, PNQ3_after, TargetPN_after, Q3Expr_t, coef_corr_mat, lm, target_PN, Target_Seq, Coefficentlist, file = "PN_linearreg_for_WGCNA.RData")

stopCluster(cl)

#Continue with WGCNA#
#2.c.1 Choosing the soft-thresholding power: analysis of network topology
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(PNQ3_reg, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
pdf("Power connectivity PNQ35+6.pdf", height = 8, width = 8)
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
softPower = 24;
TOM_PN=TOMsimilarityFromExpr(PNQ3_reg , power = 24, 
                             networkType = "signed", 
                             TOMType = "unsigned",
                             TOMDenom = "min",
                             verbose = 10)
dissTOM_PN = 1-TOM_PN

#2.b.4 Clustering using TOM
# Call the hierarchical clustering function
geneTree_PN = hclust(as.dist(1-TOM_PN), method = "average");

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(15,9)
plot(geneTree_PN, xlab="", sub="", main = "Gene clustering on power 24",
     labels = FALSE, hang = 0.04)

save(TOM_PN, file = "PN_TOM.RData")