library(dplyr)

options(stringsAsFactors = FALSE);
enableWGCNAThreads()

setwd("./midbrain/RN")

#import GeoMX normalized counts
MidQ3_Groupname<- as.data.frame(Midbrain_Q3_midbrain_NewID)
rownames(MidQ3_Groupname) <- MidQ3_Groupname$Gene
MidQ3 <- MidQ3_Groupname[, -c(1:1)]

##Import metadata
TargetMid<- as.data.frame(Target_midbrain_noSequence_clean)
rownames(TargetMid) <- TargetMid$Sample
TargetMid <- TargetMid[, -c(1:1)]

#scale covariant
Age_scaled <- scale(as.numeric(TargetMid$Age))
TargetMid$Age = as.numeric(paste(Age_scaled))

PMI_scaled <- scale(as.numeric(TargetMid$PMI))
TargetMid$PMI = as.numeric(paste(PMI_scaled))

Duration_scaled <- scale(as.numeric(TargetMid$Duration))
TargetMid$Duration = as.numeric(paste(Duration_scaled))

TargetMid$Gender[TargetMid$Gender == 'Male'] <- 'M'
TargetMid$Gender[TargetMid$Gender == 'Female'] <- 'F'

TargetMid$Ventilation[TargetMid$Ventilation == 'Yes'] <- 'Y'
TargetMid$Ventilation[TargetMid$Ventilation == 'No'] <- 'N'

TargetMid$Immunocompromised[TargetMid$Immunocompromised == 'Yes'] <- 'Y'
TargetMid$Immunocompromised[TargetMid$Immunocompromised == 'No'] <- 'N'

##Import Metaseqdata
TargetMid_MetaSeq<- as.data.frame(Target_midbrain_Seq_use)
rownames(TargetMid_MetaSeq) <-TargetMid_MetaSeq$Group
TargetMid_Seq <- TargetMid_MetaSeq[, -c(1:2)]
TargetMid_Seq <-as.data.frame(scale(TargetMid_Seq))

save(MidQ3, TargetMid, TargetMid_Seq, file="Midbrain4+4matrix+Meta.rda")
