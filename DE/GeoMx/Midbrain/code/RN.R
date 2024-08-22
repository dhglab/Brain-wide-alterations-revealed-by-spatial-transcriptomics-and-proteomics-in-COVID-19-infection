library(dplyr)
library(ggplot2)
library(EDASeq)
library(RColorBrewer)
library(Hmisc)
library(corrplot)
library(matrixStats)
library(clusterProfiler)
library(WGCNA)
library(Biobase)
library(PCAtools)
library(devtools)
library(parallel)
library(abind)
library(progress)
library(lme4)
library(factoextra)
library(FactoMineR)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

setwd("./midbrain/RN")

load("./midbrain/Midbrain4+4matrix+Meta.rda")

#Subset RN data#
target_RN <- subset(TargetMid,TargetMid$Region %in% c('RN'))
RN <- unique(rownames(target_RN))
MidQ3t <- t(MidQ3)
RNQ3t <- subset(MidQ3t, rownames(MidQ3t) %in% RN)
RNQ3 <- t(RNQ3t)
dim(RNQ3)
pca_res <- prcomp(t(RNQ3), scale. = TRUE)
autoplot(pca_res, data = target_RN, colour = 'LOT')

###sample connectivity
normadj <- (0.5+0.5*bicor(RNQ3))^2
netsummary= fundamentalNetworkConcepts(normadj)
connectivity=netsummary$Connectivity
connectivity.zscore = (connectivity-mean(connectivity))/sqrt(var(connectivity))
connectivity.plot <- data.frame(Sample.Name = names(connectivity.zscore), Z.score = connectivity.zscore, Sample.Num = 1:length(connectivity.zscore))

pdf("RNConnectivity.pdf", height = 10, width = 10)
p <- ggplot(connectivity.plot, aes(x=Sample.Num, y = Z.score, label = Sample.Name)) + geom_text(size = 4, colour = "red")
p <- p+ geom_hline(aes(yintercept = -2))
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle("Sample Connectivity")
print(p)
dev.off()

outliers <- (connectivity.zscore > mean(connectivity.zscore)+2*sd(connectivity.zscore))|(connectivity.zscore < mean(connectivity.zscore)-2*sd(connectivity.zscore))
print(colnames(RNQ3)[outliers])
#outliers 11

target_RN_after <- target_RN[!(row.names(target_RN) %in% c("11")),]
RNQ3t_after <- RNQ3t[!(row.names(RNQ3t) %in% c("11")),]
RNQ3_after <- t(RNQ3t_after)

#####linear regression####
getTopPCs <- function(x, n, col_name="") {
  scaled_dat <- x
  prcomp_dat <- prcomp(scaled_dat,center=F);
  topPCs <- prcomp_dat$rotation[,1:n];
  varexp <- (prcomp_dat$sdev)^2 / sum(prcomp_dat$sdev^2)
  topvar <- varexp[1:n]
  colnames(topPCs) <- paste0(col_name, "\n", colnames(topPCs)," (",signif(100*topvar[1:n],2),"%)")
  return(topPCs)
}

Target_RNSeq <- subset(TargetMid_Seq, subset = rownames(TargetMid_Seq) %in% rownames(target_RN_after))

seqPCs <- data.frame(getTopPCs(t(scale(Target_RNSeq)), 5))
colnames(seqPCs) <- paste0("seqpc", 1:5)
seqnames <- colnames(seqPCs)

target_RN <- cbind(target_RN_after, seqPCs)

target_RN$Gender <- as.factor(target_RN$Gender) 
target_RN$Gender <- as.numeric(target_RN$Gender)

target_RN$Immunocompromised <- as.factor(target_RN$Immunocompromised)
target_RN$Immunocompromised <- as.numeric(target_RN$Immunocompromised) 

target_RN$Ventilation <- as.factor(target_RN$Ventilation) 
target_RN$Ventilation <- as.numeric(target_RN$Ventilation) 

RNQ3Exprlog_t <- log2(RNQ3t_after)
dim(RNQ3Exprlog_t)
lmRNLog <- cbind(target_RN, RNQ3Exprlog_t)

formulaLog <- as.formula(paste0("RNQ3Exprlog_t ~ Disease + LOT + Age + PMI + Ventilation"))
lmL <- lm(formulaLog, data = lmRNLog)
Coefficentlist <- coef(summary(lmL))
Coefficentlist[[1]]

cl <- makeCluster(detectCores())
# Convert Batch to a binary matrix and export to cluster.
LOT_binary <- model.matrix(~0 + target_RN$LOT)
clusterExport(cl=cl, list('RNQ3Exprlog_t', 'target_RN','Coefficentlist', 'LOT_binary', 'seqnames','seqPCs'))
clusterEvalQ(cl=cl, gc())
crange <- 1:nrow(RNQ3)

age_effect <- parSapply(cl, crange, function(i){age_effect <- as.matrix(target_RN[, "Age"]) %*% as.matrix(Coefficentlist[[i]]["Age", "Estimate"])})

PMI_effect <- parSapply(cl, crange, function(i){pmi_effect <- as.matrix(target_RN[, "PMI"]) %*% as.matrix(Coefficentlist[[i]]["PMI", "Estimate"])})

Ventilation_effect <- parSapply(cl, crange, function(i){Vent_effect <- as.numeric(target_RN[, "Ventilation"]) %*% as.matrix(Coefficentlist[[i]]["Ventilation", "Estimate"])})

Loteffect <- parSapply(cl, crange, function(i){LOT_effect <- as.matrix(LOT_binary[, -1]) %*% as.matrix(Coefficentlist[[i]][c("LOTHWTA230003","LOTHWTA23003"), "Estimate"])})

RNQ3_regLog_lot  <- (RNQ3Exprlog_t - age_effect - Loteffect - PMI_effect - Ventilation_effect)

save(RNQ3_regLog_lot, seqPCs, seqnames, RNQ3_after, target_RN, Target_RNSeq, Coefficentlist, file = "RNlinearregression_lot.RData")

stopCluster(cl)

#PCA
pca_res <- PCA(RNQ3_regLog_lot, graph = FALSE)
eig.val <- get_eigenvalue(pca_res)
fviz_eig(pca_res, addlabels = TRUE, ylim = c(0, 50))
pdf("RNQ3_PCA.pdf", height=4, width=5)
fviz_pca_ind(pca_res ,
             geom.ind = "point", # show points only
             col.ind = target_RN$Disease, # color by groups
             palette = c("royalblue1", "coral"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
dev.off()

####Proceed to DREAM###
library(BiocParallel)
library(variancePartition)
library(ggfortify)
library(dplyr)

target_RN$Disease <- fct_relevel(target_RN$Disease, "Ctrl")#correct level#

form2 <- ~ Disease + (1|ScanID)
fit2 = dream(t(RNQ3_regLog_lot), form2,target_RN)
fit2 = eBayes(fit2)
limma_RN_reglog_lot2<- topTable(fit2, coef = NULL, number = 18677, genelist = fit2$genes,  adjust.method = "BH", sort.by = "p",resort.by = NULL,  p.value = 1, lfc = 0, confint = FALSE)

limma_RN_reglog_lot2$Mlog10adjpvalue <- -log10(limma_RN_reglog_lot2$adj.P.Val)

save(limma_RN_reglog_lot2, RNQ3_regLog_lot, target_RN, file = "limma_RN_reglog_lot2.RData")