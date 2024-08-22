library(dplyr)
library(GeneOverlap)
library(ggplot2)

setwd("./Fig1Sourcedata")
#load data FDRGM0_1, FDRPN0_1
Frontal_merge1 <- full_join(FDRGM0_1, FDRPN0_1, by = c("Target name" = "Target name"))

#Filter out the NA
Frontal_merge1_clean <- Frontal_merge1 %>%
  filter(!is.na(Mlog2FC)) %>%
  filter(!is.na(`MLog2(P-C)`))

#Pearson's product-moment correlation
GMPNPearson<-cor.test(Frontal_merge1_clean$Mlog2FC,Frontal_merge1_clean$`MLog2(P-C)`, method = "pearson")

#GeneOverlap
go.obj_GMPN <- newGeneOverlap(FDRGM0_1$`Target name`,FDRPN0_1$`Target name`,genome.size = 18702)
go.obj_GMPN<-testGeneOverlap(go.obj_GMPN)
print(go.obj_GMPN)

write.csv(Frontal_merge1_clean, "Figure1_D_inputfile.csv")

tiff("Fig1D.tiff", units="in", width= 4, height=4, res=300)
ggplot(Frontal_merge1_clean, aes(x = Mlog2FC, y = `MLog2(P-C)`)) + 
  geom_point() +
  geom_smooth(method = "lm")
dev.off()

save(FDRGM0_1, FDRPN0_1, Frontal_merge1_clean, GMPNPearson, go.obj_GMPN, file="Fig1D.rda")

setwd("./ExtendedFig1Sourcedata")
#load data FDRGM0_1, FDRnCounter0_1
Frontal_merge2 <- full_join(FDRGM0_1, FDRnCounter0_1, by = c("Target name" = "Gene"))

#Filter out the NA
Frontal_merge2_clean <- Frontal_merge2 %>%
  filter(!is.na(Mlog2FC)) %>%
  filter(!is.na(log2FoldChange))

#Pearson's product-moment correlation
GMnCounterPearson<-cor.test(Frontal_merge2_clean$Mlog2FC,Frontal_merge2_clean$log2FoldChange, method = "pearson") 

#GeneOverlap
go.obj_GMnCounter <- newGeneOverlap(FDRGM0_1$`Target name`,FDRnCounter0_1$Gene,genome.size = 18702)
go.obj_GMnCounter <-testGeneOverlap(go.obj_GMnCounter )
print(go.obj_GMnCounter)

write.csv(Frontal_merge2_clean, file = "ExtendedFig1_B1.csv")

tiff("EDFigB_1.tiff", units="in", width= 4, height=4, res=300)
ggplot(Frontal_merge2_clean, aes(x = Mlog2FC, y = log2FoldChange)) + 
  geom_point() +
  geom_smooth(method = "lm")
dev.off()

#load data FDRWM0_1, FDRnCounter0_1
Frontal_merge3 <- full_join(FDRWM0_1, FDRnCounter0_1, by = c("Target name" = "Gene"))

#Filter out the NA
Frontal_merge3_clean <- Frontal_merge3 %>%
  filter(!is.na(`MLog2(P-C)`)) %>%
  filter(!is.na(log2FoldChange))

#Pearson's product-moment correlation
WMnCounterPearson<-cor.test(Frontal_merge3_clean$`MLog2(P-C)`,Frontal_merge3_clean$log2FoldChange, method = "pearson")

#GeneOverlap
go.obj_WMnCounter <- newGeneOverlap(FDRWM0_1$`Target name`,FDRnCounter0_1$Gene,genome.size = 18702)
go.obj_WMnCounter<-testGeneOverlap(go.obj_WMnCounter)
print(go.obj_WMnCounter)

write.csv(Frontal_merge3_clean, file = "ExtendedFigure1_B2.csv")

tiff("EDFigB_2.tiff", units="in", width= 4, height=4, res=300)
ggplot(Frontal_merge3_clean, aes(x = `MLog2(P-C)`, y = log2FoldChange)) + 
  geom_point() +
  geom_smooth(method = "lm")
dev.off()

save(FDRGM0_1, FDRWM0_1, FDRnCounter0_1,Frontal_merge2_clean, Frontal_merge3_clean, go.obj_GMnCounter, go.obj_WMnCounter, file = "ExtendedFig1B.rda")

setwd("./Fig5Sourcedata")
#load data Module, FrontalProDEG, PonProDEG
module_turquiose <- subset(Module, Module == "turquoise")
FrontalProDE<- subset(FrontalProDEG, pvalue < 0.05 )
PonProDE<- subset(PonProDEG, pvalue < 0.05)

go.obj_FPTur<- newGeneOverlap(FrontalProDE$`Gene Name`,module_turquiose$Gene,genome.size = 18702)
go.obj_FPTur<-testGeneOverlap(go.obj_FPTur)
print(go.obj_FPTur)

go.obj_PPTur <- newGeneOverlap(PonProDE$`Gene Name` ,module_turquiose$Gene,genome.size = 18702)
go.obj_PPTur<-testGeneOverlap(go.obj_PPTur)
print(go.obj_PPTur)

module_greenyellow <- subset(Module, Module == "greenyellow")
go.obj_FPGY <- newGeneOverlap(FrontalProDE$`Gene Name`,module_greenyellow$Gene,genome.size = 18702)
go.obj_FPGY<-testGeneOverlap(go.obj_FPGY)
print(go.obj_FPGY)

go.obj_PPGY <- newGeneOverlap(PonProDE$`Gene Name` ,module_greenyellow$Gene,genome.size = 18702)
go.obj_PPGY<-testGeneOverlap(go.obj_PPGY)
print(go.obj_PPGY)

module_green<- subset(Module, Module == "green")
go.obj_FPGreen <- newGeneOverlap(FrontalProDE$`Gene Name`,module_green$Gene,genome.size = 18702)
go.obj_FPGreen<-testGeneOverlap(go.obj_FPGreen)
print(go.obj_FPGreen)

go.obj_PPGreen <- newGeneOverlap(PonProDE$`Gene Name` ,module_green$Gene,genome.size = 18702)
go.obj_PPGreen <-testGeneOverlap(go.obj_PPGreen)
print(go.obj_PPGreen)

module_blue<- subset(Module, Module == "blue")
go.obj_FPBlue<- newGeneOverlap(FrontalProDE$`Gene Name`,module_blue$Gene,genome.size = 18702)
go.obj_FPBlue<-testGeneOverlap(go.obj_FPBlue)
print(go.obj_FPBlue)

go.obj_PPBlue<- newGeneOverlap(PonProDE$`Gene Name` ,module_blue$Gene,genome.size = 18702)
go.obj_PPBlue<-testGeneOverlap(go.obj_PPBlue)
print(go.obj_PPBlue)

save(Module, FrontalProDEG, PonProDEG, FrontalProDE, PonProDE,module_blue, module_green, module_greenyellow, module_turquiose, go.obj_PPTur, go.obj_FPTur, go.obj_FPGY, go.obj_FPGreen, go.obj_FPBlue, go.obj_PPGY, go.obj_PPGreen, go.obj_PPBlue, file = "Fig5C.rda")

setwd("./ExtendedFig6Sourcedata")
#load FDRGM05,Henghypoxia_2d_7d_500
FDRGM05 <- FDRGM05[order(FDRGM05$Mlog2FC, decreasing = TRUE),]
FDRGM05_up <- FDRGM05[1:500,]
go.obj_GMHypo2d <- newGeneOverlap(FDRGM05_up$`Target name`,
                            Heng_hypoxia_2d_7d_500$Hypoxia_2d,
                            genome.size = 18702)
go.obj_GMHypo2d <-testGeneOverlap(go.obj_GMHypo2d)
print(go.obj_GMHypo2d)

go.obj_GMHypo7d <- newGeneOverlap(FDRGM05_up$`Target name`,
                                  Heng_hypoxia_2d_7d_500$Hypoxia_7d,
                                  genome.size = 18702)
go.obj_GMHypo7d <-testGeneOverlap(go.obj_GMHypo7d)
print(go.obj_GMHypo7d)

#load data FDRPN05, FDRSN05, ZHang_hypoxiabrainDEG   
Zhang_up <- subset(ZHang_hypoxiabrainDEG, REGULATION == "UP")

go.obj_GMZhangup <- newGeneOverlap(FDRGM05_up$`Target name`,
                            Zhang_up$GENE_ID,
                            genome.size = 18702)
go.obj_GMZhangup <-testGeneOverlap(go.obj_GMZhangup)
print(go.obj_GMZhangup)

FDRPN05_up <- FDRPN05[1404:1903,]
go.obj_PNZhangup <- newGeneOverlap(FDRPN05_up$`Target name`,
                            Zhang_up$GENE_ID,
                            genome.size = 18702)
go.obj_PNZhangup <-testGeneOverlap(go.obj_PNZhangup)
print(go.obj_PNZhangup)

FDRSN05_up <- FDRSN05[1:500,]
go.obj_SNZhangup <- newGeneOverlap(FDRSN05_up$`Target name`,
                            Zhang_up$GENE_ID,
                            genome.size = 18702)
go.obj_SNZhangup <-testGeneOverlap(go.obj_SNZhangup)
print(go.obj_SNZhangup)

Zhang_down <- subset(ZHang_hypoxiabrainDEG, REGULATION == "DOWN")
FDRGM05_down <- FDRGM05[3932:4431,]
go.obj_GMZhangdown <- newGeneOverlap(FDRGM05_down$`Target name`,
                            Zhang_down$GENE_ID,
                            genome.size = 18702)
go.obj_GMZhangdown <-testGeneOverlap(go.obj_GMZhangdown)
print(go.obj_GMZhangdown)

FDRPN05_down <- FDRPN05[1:500,]
go.obj_PNZhangdown <- newGeneOverlap(FDRPN05_down$`Target name`,
                            Zhang_down$GENE_ID,
                            genome.size = 18702)
go.obj_PNZhangdown <-testGeneOverlap(go.obj_PNZhangdown)
print(go.obj_PNZhangdown)

FDRSN05_down <- FDRSN05[12006:12505,]
go.obj_SNZhangdown <- newGeneOverlap(FDRSN05_down$`Target name`,
                            Zhang_down$GENE_ID,
                            genome.size = 18702)
go.obj_SNZhangdown <-testGeneOverlap(go.obj_SNZhangdown)
print(go.obj_SNZhangdown)

save(Heng_hypoxia_2d_7d_500, ZHang_hypoxiabrainDEG, Zhang_down, Zhang_up, FDRGM05, FDRPN05, FDRSN05, go.obj_GMHypo2d,go.obj_GMHypo7d, FDRGM05_up, FDRGM05_down, FDRPN05_down, FDRPN05_up, FDRSN05_down, FDRSN05_up, go.obj_GMZhangdown, go.obj_GMZhangup, go.obj_PNZhangdown, go.obj_PNZhangup, go.obj_SNZhangdown, go.obj_SNZhangup, file = "ExtendedFig6.rda")
