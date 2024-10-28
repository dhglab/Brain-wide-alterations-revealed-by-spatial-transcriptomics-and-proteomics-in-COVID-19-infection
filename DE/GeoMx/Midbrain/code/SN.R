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
library(parallel)
library(abind)
library(progress)
library(lme4)
library(Matrix)
library(PCAtools)
library(ggfortify)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

setwd("./midbrain/SN")
load("./midbrain/Midbrain4+4matrix+Meta.rda")

##SubsetSN data##
target_SN <- subset(TargetMid,TargetMid$Region %in% c('SN'))
SN <- unique(rownames(target_SN))
MidQ3t <- t(MidQ3)
SNQ3t <- subset(MidQ3t, rownames(MidQ3t) %in% SN)
SNQ3 <- t(SNQ3t)
dim(SNQ3)
pca_res <- prcomp(t(SNQ3), scale. = TRUE)
autoplot(pca_res, data = target_SN, colour = 'LOT')

###sample connectivity
normadj <- (0.5+0.5*bicor(SNQ3))^2
netsummary= fundamentalNetworkConcepts(normadj)
connectivity=netsummary$Connectivity
connectivity.zscore = (connectivity-mean(connectivity))/sqrt(var(connectivity))
connectivity.plot <- data.frame(Sample.Name = names(connectivity.zscore), Z.score = connectivity.zscore, Sample.Num = 1:length(connectivity.zscore))

pdf("SNConnectivity.pdf", height = 10, width = 10)
p <- ggplot(connectivity.plot, aes(x=Sample.Num, y = Z.score, label = Sample.Name)) + geom_text(size = 4, colour = "red")
p <- p+ geom_hline(aes(yintercept = -2))
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle("Sample Connectivity")
print(p)
dev.off()

outliers <- (connectivity.zscore > mean(connectivity.zscore)+2*sd(connectivity.zscore))|(connectivity.zscore < mean(connectivity.zscore)-2*sd(connectivity.zscore))
print(colnames(SNQ3)[outliers])
####outliers 0

#####linear regression#

getTopPCs <- function(x, n, col_name="") {
  scaled_dat <- x
  prcomp_dat <- prcomp(scaled_dat,center=F);
  topPCs <- prcomp_dat$rotation[,1:n];
  varexp <- (prcomp_dat$sdev)^2 / sum(prcomp_dat$sdev^2)
  topvar <- varexp[1:n]
  colnames(topPCs) <- paste0(col_name, "\n", colnames(topPCs)," (",signif(100*topvar[1:n],2),"%)")
  return(topPCs)
}

Target_SNSeq <- subset(TargetMid_Seq, subset = rownames(TargetMid_Seq) %in% SN)

seqPCs <- data.frame(getTopPCs(t(scale(Target_SNSeq)), 5))
colnames(seqPCs) <- paste0("seqpc", 1:5)
seqnames <- colnames(seqPCs)

target_SN <- cbind(target_SN, seqPCs)

target_SN$Gender <- as.factor(target_SN$Gender) 
target_SN$Gender <- as.numeric(target_SN$Gender)

target_SN$Immunocompromised <- as.factor(target_SN$Immunocompromised)
target_SN$Immunocompromised <- as.numeric(target_SN$Immunocompromised) 

target_SN$Ventilation <- as.factor(target_SN$Ventilation) 
target_SN$Ventilation <- as.numeric(target_SN$Ventilation) 

SNQ3Exprlog_t <- log2(SNQ3t)
dim(SNQ3Exprlog_t)
lmSNLog <- cbind(target_SN, SNQ3Exprlog_t)

formulaLog <- as.formula(paste0("SNQ3Exprlog_t ~ Disease + LOT + Age + PMI + Ventilation"))
lmL <- lm(formulaLog, data = lmSNLog)
Coefficentlist <- coef(summary(lmL))
Coefficentlist[[1]]

cl <- makeCluster(detectCores())
# Convert Batch to a binary matrix and export to cluster.
LOT_binary <- model.matrix(~0 + target_SN$LOT)
clusterExport(cl=cl, list('SNQ3Exprlog_t', 'target_SN', 'Coefficentlist', 'LOT_binary','seqnames','seqPCs'))
clusterEvalQ(cl=cl, gc())
crange <- 1:nrow(SNQ3)

Age_effect <- parSapply(cl, crange, function(i){age_effect <- as.matrix(target_SN[, "Age"]) %*% as.matrix(Coefficentlist[[i]]["Age", "Estimate"])})

PMI_effect <- parSapply(cl, crange, function(i){pmi_effect <- as.matrix(target_SN[, "PMI"]) %*% as.matrix(Coefficentlist[[i]]["PMI", "Estimate"])})

Ventilation_effect <- parSapply(cl, crange, function(i){Vent_effect <- as.numeric(target_SN[, "Ventilation"]) %*% as.matrix(Coefficentlist[[i]]["Ventilation", "Estimate"])})

Lot_effect <- parSapply(cl, crange, function(i){LOT_effect <- as.matrix(LOT_binary[, -1]) %*% as.matrix(Coefficentlist[[i]][c("LOTHWTA230003","LOTHWTA23003"), "Estimate"])})

SNQ3_reg_lot <- (SNQ3Exprlog_t - Lot_effect - Age_effect - PMI_effect - Ventilation_effect)

stopCluster(cl)

#PCA
library("factoextra")
library("FactoMineR")
library(forcats)

target_SN$Disease <- fct_relevel(target_SN$Disease, "Ctrl")#correct level#

pca_res <- PCA(SNQ3_reg_lot, graph = FALSE)
eig.val <- get_eigenvalue(pca_res)
fviz_eig(pca_res, addlabels = TRUE, ylim = c(0, 50))
pdf("SNQ3_PCA.pdf", width=5, height=4)
fviz_pca_ind(pca_res ,
             geom.ind = "point", # show points only
             col.ind = target_SN$Disease, # color by groups
             palette = c("royalblue1", "coral"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",
             title = "Substantia Nigra",
)
dev.off()

##############Proceed to DREAM###################
library(variancePartition)
library(ggfortify)
library(dplyr)

target_SN$Disease <- fct_relevel(target_SN$Disease, "Ctrl")

form2 <- ~ Disease + (1|ScanID)
fit2 = dream(t(SNQ3_reg_lot), form2,target_SN)
fit2 = eBayes(fit2)
limma_SN_reglog_lot2<- topTable(fit2, coef = NULL, number = 18677, genelist = fit2$genes,  adjust.method = "BH", sort.by = "p",resort.by = NULL,  p.value = 1, lfc = 0, confint = FALSE)

limma_SN_reglog_lot2$MlogFC <- -limma_SN_reglog_lot2$logFC
limma_SN_reglog_lot2$Mlog10adjpvalue <- -log10(limma_SN_reglog_lot2$adj.P.Val)
save(limma_SN_reglog_lot2, SNQ3_reg_lot, target_SN, file = "limma_SN_reglog_lot2.RData")
