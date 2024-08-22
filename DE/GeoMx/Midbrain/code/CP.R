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

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

setwd("./midbrain/CP")

load("./midbrain/Midbrain4+4matrix+Meta.rda")

#Subset CP data##
target_CP <- subset(TargetMid,TargetMid$Region %in% c('CP'))
CP <- unique(rownames(target_CP))
MidQ3t <- t(MidQ3)
CPQ3t <- subset(MidQ3t, rownames(MidQ3t) %in% CP)
CPQ3 <- t(CPQ3t)
dim(CPQ3)
pca_res <- prcomp(t(CPQ3), scale. = TRUE)
autoplot(pca_res, data = target_CP, colour = 'LOT')

###sample connectivity
normadj <- (0.5+0.5*bicor(CPQ3))^2
netsummary= fundamentalNetworkConcepts(normadj)
connectivity=netsummary$Connectivity
connectivity.zscore = (connectivity-mean(connectivity))/sqrt(var(connectivity))
connectivity.plot <- data.frame(Sample.Name = names(connectivity.zscore), Z.score = connectivity.zscore, Sample.Num = 1:length(connectivity.zscore))

pdf("CPConnectivity.pdf", height = 10, width = 10)
p <- ggplot(connectivity.plot, aes(x=Sample.Num, y = Z.score, label = Sample.Name)) + geom_text(size = 4, colour = "red")
p <- p+ geom_hline(aes(yintercept = -2))
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle("Sample Connectivity")
print(p)
dev.off()

outliers <- (connectivity.zscore > mean(connectivity.zscore)+2*sd(connectivity.zscore))|(connectivity.zscore < mean(connectivity.zscore)-2*sd(connectivity.zscore))
print(colnames(CPQ3)[outliers])
####outliers 0

#####linear regression#######
target_CP$Gender <- as.factor(target_CP$Gender) 
target_CP$Gender <- as.numeric(target_CP$Gender)

target_CP$Immunocompromised <- as.factor(target_CP$Immunocompromised)
target_CP$Immunocompromised <- as.numeric(target_CP$Immunocompromised) 

target_CP$Ventilation <- as.factor(target_CP$Ventilation) 
target_CP$Ventilation <- as.numeric(target_CP$Ventilation) 

CPQ3Exprlog_t <- log2(CPQ3t)
dim(CPQ3Exprlog_t)
lmCPLog <- cbind(target_CP, CPQ3Exprlog_t)

formulaLog <- as.formula(paste0("CPQ3Exprlog_t ~ Disease + LOT + Age + PMI + Ventilation + Gender + Duration + Immunocompromised"))
lmL <- lm(formulaLog, data = lmCPLog)
Coefficentlist <- coef(summary(lmL))
Coefficentlist[[1]]
#Only LOT and Age

cl <- makeCluster(detectCores())
LOT_binary <- model.matrix(~0 + target_CP$LOT)

clusterExport(cl=cl, list('CPQ3Exprlog_t', 'target_CP','Coefficentlist', 'LOT_binary'))
clusterEvalQ(cl=cl, gc())
crange <- 1:nrow(CPQ3)

age_effect <- parSapply(cl, crange, function(i){age_effect <- as.matrix(target_CP[, "Age"]) %*% as.matrix(Coefficentlist[[i]]["Age", "Estimate"])})

Loteffect <- parSapply(cl, crange, function(i){LOT_effect <- as.matrix(LOT_binary[, -1]) %*% as.matrix(Coefficentlist[[i]][c("LOTHWTA230003","LOTHWTA23003"), "Estimate"])})

CPQ3_regLog  <- (CPQ3Exprlog_t - age_effect - Loteffect)

save(CPQ3_regLog, CPQ3, target_CP, Target_CPSeq, Coefficentlist, file = "CPlinearregression.RData")

stopCluster(cl)

#PCA
pca_res <- PCA(CPQ3_regLog, graph = FALSE)
eig.val <- get_eigenvalue(pca_res)
fviz_eig(pca_res, addlabels = TRUE, ylim = c(0, 50))
pdf("RNQ3_PCA.pdf", height=4, width=5)
fviz_pca_ind(pca_res ,
             geom.ind = "point", # show points only
             col.ind = target_CP$Disease, # color by groups
             palette = c("royalblue1", "coral"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
dev.off()

###################################Proceed to DREAM#####################################
library(BiocParallel)
library(variancePartition)
library(ggfortify)
library(dplyr)

target_CP$Disease <- fct_relevel(target_CP$Disease, "Ctrl")#correct level#

form2 <- ~ Disease + (1|ScanID)
fit2 = dream(t(CPQ3_regLog), form2,target_CP)
fit2 = eBayes(fit2)
limma_CP_reglog2<- topTable(fit2, coef = NULL, number = 18677, genelist = fit2$genes,  adjust.method = "BH", sort.by = "p",resort.by = NULL,  p.value = 1, lfc = 0, confint = FALSE)

limma_CP_reglog2$Mlog10adjpvalue <- -log10(limma_CP_reglog2$adj.P.Val)

save(limma_CP_reglog2, CPQ3_regLog, target_CP, file = "limma_CP_reglog2.RData")
