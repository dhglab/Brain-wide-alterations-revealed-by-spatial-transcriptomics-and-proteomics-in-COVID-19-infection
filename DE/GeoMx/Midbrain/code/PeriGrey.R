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

setwd("./midbrain/Periaqueduct")
load("./midbrain/Midbrain4+4matrix+Meta.rda")

#Subset Peri data#
target_Peri <- subset(TargetMid,TargetMid$Region %in% c('Periaqueduct'))
Peri <- unique(rownames(target_Peri))
MidQ3t <- t(MidQ3)
PeriQ3t <- subset(MidQ3t, rownames(MidQ3t) %in% Peri)
PeriQ3 <- t(PeriQ3t)
dim(PeriQ3)
pca_res <- prcomp(t(PeriQ3), scale. = TRUE)
autoplot(pca_res, data = target_Peri, colour = 'LOT')

###sample connectivity
normadj <- (0.5+0.5*bicor(PeriQ3))^2
netsummary= fundamentalNetworkConcepts(normadj)
connectivity=netsummary$Connectivity
connectivity.zscore = (connectivity-mean(connectivity))/sqrt(var(connectivity))
connectivity.plot <- data.frame(Sample.Name = names(connectivity.zscore), Z.score = connectivity.zscore, Sample.Num = 1:length(connectivity.zscore))

pdf("PeriConnectivity.pdf", height = 10, width = 10)
p <- ggplot(connectivity.plot, aes(x=Sample.Num, y = Z.score, label = Sample.Name)) + geom_text(size = 4, colour = "red")
p <- p+ geom_hline(aes(yintercept = -2))
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle("Sample Connectivity")
print(p)
dev.off()

outliers <- (connectivity.zscore > mean(connectivity.zscore)+2*sd(connectivity.zscore))|(connectivity.zscore < mean(connectivity.zscore)-2*sd(connectivity.zscore))
print(colnames(PeriQ3)[outliers])
####outliers 0

#####linear regression#######
getTopPCs <- function(x, n, col_name="") {
  scaled_dat <- x
  prcomp_dat <- prcomp(scaled_dat,center=F);
  topPCs <- prcomp_dat$rotation[,1:n];
  varexp <- (prcomp_dat$sdev)^2 / sum(prcomp_dat$sdev^2)
  topvar <- varexp[1:n]
  colnames(topPCs) <- paste0(col_name, "\n", colnames(topPCs)," (",signif(100*topvar[1:n],2),"%)")
  return(topPCs)
}

Target_PeriSeq <- subset(TargetMid_Seq, subset = rownames(TargetMid_Seq) %in% rownames(target_Peri))

seqPCs <- data.frame(getTopPCs(t(scale(Target_PeriSeq)), 5))
colnames(seqPCs) <- paste0("seqpc", 1:5)
seqnames <- colnames(seqPCs)
target_Peri <- cbind(target_Peri, seqPCs)

target_Peri$Gender <- as.factor(target_Peri$Gender) 
target_Peri$Gender <- as.numeric(target_Peri$Gender)

target_Peri$Immunocompromised <- as.factor(target_Peri$Immunocompromised)
target_Peri$Immunocompromised <- as.numeric(target_Peri$Immunocompromised) 

target_Peri$Ventilation <- as.factor(target_Peri$Ventilation) 
target_Peri$Ventilation <- as.numeric(target_Peri$Ventilation) 

PeriQ3Exprlog_t <- log2(PeriQ3t)
dim(PeriQ3Exprlog_t)
lmPeriLog <- cbind(target_Peri, PeriQ3Exprlog_t)

formulaLog <- as.formula(paste0("PeriQ3Exprlog_t ~ Disease + LOT + Age + PMI + Ventilation + Gender + Duration + Immunocompromised"))
lmL <- lm(formulaLog, data = lmPeriLog)
Coefficentlist <- coef(summary(lmL))
Coefficentlist[[1]]
#Only LOT and Age

cl <- makeCluster(detectCores())
LOT_binary <- model.matrix(~0 + target_Peri$LOT)
clusterExport(cl=cl, list('PeriQ3Exprlog_t', 'target_Peri','Coefficentlist', 'LOT_binary', 'seqnames','seqPCs'))
clusterEvalQ(cl=cl, gc())
crange <- 1:nrow(PeriQ3)

age_effect <- parSapply(cl, crange, function(i){age_effect <- as.matrix(target_Peri[, "Age"]) %*% as.matrix(Coefficentlist[[i]]["Age", "Estimate"])})

Loteffect <- parSapply(cl, crange, function(i){LOT_effect <- as.matrix(LOT_binary[, -1]) %*% as.matrix(Coefficentlist[[i]][c("LOTHWTA230003","LOTHWTA23003"), "Estimate"])})

PeriQ3_regLog  <- (PeriQ3Exprlog_t - age_effect - Loteffect)

save(PeriQ3_regLog, PeriQ3, target_Peri, Target_PeriSeq, Coefficentlist, file = "Perilinearregression.RData")

stopCluster(cl)

#PCA
library("factoextra")
library("FactoMineR")
library(forcats)

target_Peri$Disease <- fct_relevel(target_Peri$Disease, "Ctrl")

pca_res <- PCA(PeriQ3_regLog, graph = FALSE)
eig.val <- get_eigenvalue(pca_res)
fviz_eig(pca_res, addlabels = TRUE, ylim = c(0, 50))
pdf("PeriQ3_PCA.pdf", height=4, width=5)
fviz_pca_ind(pca_res ,
             geom.ind = "point", # show points only
             col.ind = target_Peri$Disease, # color by groups
             palette = c("royalblue1", "coral"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",
             title = "PeriGM",
)
dev.off()

##############Proceed to DREAM###################
library(BioPeriarallel)
library(variancePartition)
library(ggfortify)
library(dplyr)

target_Peri$Disease <- fct_relevel(target_Peri$Disease, "Ctrl")

form2 <- ~ Disease + (1|ScanID)
fit2 = dream(t(PeriQ3_regLog), form2,target_Peri)
fit2 = eBayes(fit2)
limma_Peri_reglog2<- topTable(fit2, coef = NULL, number = 18677, genelist = fit2$genes,  adjust.method = "BH", sort.by = "p",resort.by = NULL,  p.value = 1, lfc = 0, confint = FALSE)

limma_Peri_reglog2$Mlog10adjpvalue <- -log10(limma_Peri_reglog2$adj.P.Val)

save(limma_Peri_reglog2, PeriQ3_regLog, target_Peri, file = "limma_Peri_reglog2.RData")
