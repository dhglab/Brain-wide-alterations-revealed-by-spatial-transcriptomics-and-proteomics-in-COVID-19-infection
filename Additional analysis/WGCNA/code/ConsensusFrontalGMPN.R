library(WGCNA)
library(igraph)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

setwd("./Fig2Sourcedata/WGCNA/consensusModule")

#Load TOM from FrontalGM and PontineNuclei_linear_regression
consensusTOM = pmin(TOM,TOM_PN)
save (consensusTOM, file = "consensusTOM for FrontalGM_PN.RData")

dissTOM = 1-consensusTOM

#2.b.4 Clustering using TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
sizeGrWindow(15,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on Frontal12PN24",
     labels = FALSE, hang = 0.04)

########Module selection##########
mColorh=NULL
for (ds in 0:4){
  tree = cutreeHybrid(dendro = geneTree, pamStage=FALSE,
                      minClusterSize = 50, cutHeight = 0.9999,
                      deepSplit = ds, distM = dissTOM)
  mColorh=cbind(mColorh,labels2colors(tree$labels));
}
pdf("FrontalGM_PN_Module_choices_12_24_50.pdf", height=10,width=25);
plotDendroAndColors(geneTree, mColorh, paste("dpSplt =",0:4), main = "FrontalGM12_PN24",dendroLabels=FALSE)
dev.off()

modules = mColorh[,5] #(Chosen based on plot below)
pdf("FrontalGM_PN_Module.pdf", height=6, width = 18);
plotDendroAndColors(geneTree, modules, main ="FrontalGM_PN", dendroLabels = FALSE)
dev.off()

#Extended Fig.3A
tiff("FrontalGM_PN_Module.tiff", units="in", width=7, height=5, res=300)
plotDendroAndColors(geneTree, modules, main ="FrontalGM_PN", dendroLabels = FALSE)
dev.off()

save(geneTree, mColorh, modules, tree, file = "FrontalGM and PN consensus_module_choice.RData")

########consensusKME############
#Load FrontalGMQ3_reg, PNQ3_reg
multiExpr = list(list(data=t(FrontalGMQ3_reg)), list(data=PNQ3_reg))
names(multiExpr) <- c("FrontalGM", "PNQ3")
moduleLabels <- modules
consensus_kM <- consensusKME(
  multiExpr,
  moduleLabels, 
  multiEigengenes = NULL, 
  consensusQuantile = 0, 
  signed = TRUE,
  useModules = NULL,
  metaAnalysisWeights = NULL,
  corAndPvalueFnc = corAndPvalue, corOptions = list(), corComponent = "cor",
  useRankPvalue = TRUE,
  rankPvalueOptions = list(pValueMethod = "scale"),
  setNames = NULL, 
  excludeGrey = TRUE, 
  greyLabel = if (is.numeric(moduleLabels)) 0 else "grey")
consensus_kM<- as.data.frame(consensus_kM)
consensus_kM_s <- consensus_kM[, grep("consensus", colnames(consensus_kM))]
consensus_kM_s$module <- modules
consensus_kM_s$Gene <- rownames(consensus_kM_s)
write.csv(consensus_kM_s, file="consensuskM_table.csv")

individual_kM <- consensus_kM[, grep("Z.", colnames(consensus_kM))]
individual_FrontalkM <-consensus_kM[, grep("Frontal", colnames(consensus_kM))]
individual_PNKM <-consensus_kM[, grep("PNQ3", colnames(consensus_kM))]
write.csv(individual_FrontalkM, file="kM_FrontalGM.csv")
write.csv(individual_PNKM, file="kM_PN.csv")
save(consensus_kM,consensus_kM_s, individual_FrontalkM, individual_PNKM, file="consensus_kM.RData")


#############Differential module expression analysis################
#Use the consensus modules to calculate eigengenes#
MEList_Fc = moduleEigengenes(t(FrontalGMQ3_reg), colors = modules)
MEs_Fc = MEList_Fc$eigengenes

MEList_Pc = moduleEigengenes(PNQ3_reg, colors = modules)
MEs_Pc = MEList_Pc$eigengenes

save(geneTree, mColorh, modules, tree,FrontalGMQ3_reg, target_Frontal, PNQ3_reg, target_PN, MEList_Fc, MEs_Fc, MEList_Pc, MEs_Pc, file = "FrontalGM and PN consensus_module_choice.RData")

#########FrontalGM##########
New_Fc<- cbind(target_Frontal, MEs_Fc)
ME_names_Fc <- colnames(MEs_Fc)

lm_results_Fc <- lapply(ME_names_Fc, 
                     function(ME){lm(as.formula(paste0(ME,"~ Disease + Age + PMI + Gender + Ventilation + Duration + Immunocompromised")), data = New_Fc)
                     }
)

lapply(lm_results_Fc,summary)

Res_1<-as.data.frame(coefficients(summary(lm_results_Fc[[1]])))
Res_1<- Res_1[!(row.names(Res_1) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[2]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_1, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[3]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[4]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[5]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[6]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[7]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[8]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[9]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[10]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[11]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[12]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[13]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[14]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[15]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[16]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[17]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[18]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[19]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[20]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[21]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[22]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Fc[[23]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Fc<-rbind(Res_Fc, Res_2)

colnames(Res_Fc)<- c("Beta","Std","tvalue","pvalue")
rownames(Res_Fc)<- ME_names_Fc

library(dplyr)
Res_filtered_Fc <- filter(Res_Fc, pvalue < 0.05)
Res_filtered_Fc$MEColors = rownames(Res_filtered_Fc)
Res_filtered_grey_Fc <- Res_filtered_Fc[!Res_filtered_Fc$MEColors %in% c("MEgrey"),]
Res_filtered_grey_Fc$Area <- c("FrontalGM")

save(Res_filtered_Fc,Res_filtered_grey_Fc, lm_results_Fc, ME_names_Fc, New_Fc, file = "DifferentialFrontalGMmodules_consensus.RData")

#########Pons#############
New_Pc<- cbind(target_PN, MEs_Pc)
ME_names_Pc <- colnames(MEs_Pc)

lm_results_pc <- lapply(ME_names_Pc, 
                        function(ME){lm(as.formula(paste0(ME,"~ Disease + Age + PMI + 
                                                          Gender + Ventilation + Duration 
                                                          + Immunocompromised")), data = New_Pc)
                        }
)

lapply(lm_results_Pc,summary)

####Pontine nuclei
New_Pc<- cbind(target_PN, MEs_Pc)
ME_names_Pc <- colnames(MEs_Pc)

lm_results_Pc <- lapply(ME_names_Pc, 
                        function(ME){lm(as.formula(paste0(ME,"~ Disease + Age + PMI + Gender + Ventilation + Duration + Immunocompromised")), data = New_Fc)
                        }
)

lapply(lm_results_Pc,summary)

Res_1<-as.data.frame(coefficients(summary(lm_results_Pc[[1]])))
Res_1<- Res_1[!(row.names(Res_1) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[2]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_1, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[3]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[4]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[5]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[6]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[7]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[8]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[9]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[10]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[11]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[12]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[13]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[14]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[15]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[16]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[17]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[18]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[19]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[20]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[21]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[22]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)
Res_2<-as.data.frame(coefficients(summary(lm_results_Pc[[23]])))
Res_2<- Res_2[!(row.names(Res_2) %in% c("(Intercept)","Age","PMI","Gender","Ventilation","Duration","Immunocompromised")),]
Res_Pc<-rbind(Res_Pc, Res_2)

colnames(Res_Pc)<- c("Beta","Std","tvalue","pvalue")
rownames(Res_Pc)<- ME_names_Pc

library(dplyr)
Res_filtered_Pc <- filter(Res_Pc, pvalue < 0.05)
Res_filtered_Pc$MEColors = rownames(Res_filtered_Pc)
Res_filtered_grey_Pc <- Res_filtered_Pc[!Res_filtered_Pc$MEColors %in% c("MEgrey"),]
Res_filtered_grey_Pc$Area <- c("PontineNuclei")

save(Res_filtered_Pc,Res_filtered_grey_Pc, lm_results_Pc, ME_names_Pc, New_Pc, file = "DifferentialPNmodules_consensus.RData")

write.csv(Res_Fc,Res_Pc, file = "FrontalGMPN_consensusmodules_DE_lm.csv")


