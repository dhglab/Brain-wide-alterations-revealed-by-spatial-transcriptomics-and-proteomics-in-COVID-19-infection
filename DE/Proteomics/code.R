
# 2021.06.03_Frontal

library(MSstatsTMT)

evidence <- read.table("evidence.txt",sep = '\t',header = TRUE)

unique(evidence$Raw.file)

proteinGroups<-read.table("proteinGroups_Frontal-6cif.txt",sep = '\t',header = TRUE)

annotation.mq<-read.csv("annotation.csv")

?MaxQtoMSstatsTMTFormat

input.mq<-MaxQtoMSstatsTMTFormat(evidence, proteinGroups,annotation = annotation.mq)

?proteinSummarization
#msstats means medianpolish+impute, no internal normalize channel, 
#so reference norm=false. use this 06.03.2021
protein.summa<- proteinSummarization(data = input.mq, method = "msstats", global_norm = TRUE,
                                     reference_norm = FALSE)

#same with above but more information writen here, didn't use at 03.08
protein.summa<- proteinSummarization(data = input.mq, method = "msstats", global_norm = TRUE,
                                     reference_norm = FALSE, remove_norm_channel = TRUE,
                                     remove_empty_channel = TRUE, MBimpute = TRUE,
                                     maxQuantileforCensored = 0.999)

head(protein.summa)

write.csv(protein.summa, file = "protein.summa.csv")

?dataProcessPlotsTMT

####if do not want to draw the plots, do not need run this;
#profileplot shows peptides and summarization information together.
dataProcessPlotsTMT(data.peptide =input.mq,
                    data.summarization = protein.summa,
                    type ="profileplot" )

#didn't draw this time, this should be box plot. 06.18.2020
dataProcessPlotsTMT(data.peptide=input.mq,
                    data.summarization=protein.summa,
                    type='QCPlot')

unique(input.mq$Condition)

#contrast set test, but didn't use eventually
?groupComparisonTMT
comparison<-matrix(c(-1,1),nrow = 1)
final.compar <- groupComparisonTMT(data = protein.summa,
                                   contrast.matrix = comparison,
                                   moderated = TRUE)

final.compar <- groupComparisonTMT(data = protein.summa,
                                   contrast.matrix = "pairwise",
                                   moderated = TRUE)

View(final.compar)
#didn't write this
sig.final <- final.compar[final.compar$adj.pvalue<=0.05,]

save(final.compar, file = "final.compar.csv")

write.csv(final.compar, file = "final.compar.csv")




#2021.06.03 hippocampus

library(MSstatsTMT)

evidence <- read.table("evidence.txt",sep = '\t',header = TRUE)

unique(evidence$Raw.file)

proteinGroups<-read.table("proteinGroups.txt",sep = '\t',header = TRUE)

annotation.mq<-read.csv("annotation.csv")

input.mq<-MaxQtoMSstatsTMTFormat(evidence, proteinGroups,annotation = annotation.mq)

#msstats means medianpolish+impute, no internal normalize channel, 
#so reference norm=false. use this 06.03.2021
protein.summa<- proteinSummarization(data = input.mq, method = "msstats", global_norm = TRUE,
                                     reference_norm = FALSE)

head(protein.summa)

write.csv(protein.summa, file = "protein.summa.csv")

unique(input.mq$Condition)

final.compar <- groupComparisonTMT(data = protein.summa,
                                   contrast.matrix = "pairwise",
                                   moderated = TRUE)

View(final.compar)

save(final.compar, file = "final.compar.csv")

write.csv(final.compar, file = "final.compar.csv")


#2021.06.03 basal ganglia

library(MSstatsTMT)

evidence <- read.table("evidence.txt",sep = '\t',header = TRUE)

unique(evidence$Raw.file)

proteinGroups<-read.table("proteinGroups.txt",sep = '\t',header = TRUE)

annotation.mq<-read.csv("annotation .csv")

input.mq<-MaxQtoMSstatsTMTFormat(evidence, proteinGroups,annotation = annotation.mq)

#msstats means medianpolish+impute, no internal normalize channel, 
#so reference norm=false. use this 06.03.2021
protein.summa<- proteinSummarization(data = input.mq, method = "msstats", global_norm = TRUE,
                                     reference_norm = FALSE)

head(protein.summa)

write.csv(protein.summa, file = "protein.summa.csv")

unique(input.mq$Condition)

final.compar <- groupComparisonTMT(data = protein.summa,
                                   contrast.matrix = "pairwise",
                                   moderated = TRUE)

View(final.compar)

save(final.compar, file = "final.compar.csv")

write.csv(final.compar, file = "final.compar.csv")



#2021.06.21 Midbrains

library(MSstatsTMT)

evidence <- read.table("evidence.txt",sep = '\t',header = TRUE)

unique(evidence$Raw.file)

proteinGroups<-read.table("proteinGroups.txt",sep = '\t',header = TRUE)

annotation.mq<-read.csv("annotation.csv")

input.mq<-MaxQtoMSstatsTMTFormat(evidence, proteinGroups,annotation = annotation.mq)

#msstats means medianpolish+impute, no internal normalize channel, 
#so reference norm=false. use this 06.21.2021
protein.summa<- proteinSummarization(data = input.mq, method = "msstats", global_norm = TRUE,
                                     reference_norm = FALSE)

head(protein.summa)

write.csv(protein.summa, file = "protein.summa.csv")

unique(input.mq$Condition)

final.compar <- groupComparisonTMT(data = protein.summa,
                                   contrast.matrix = "pairwise",
                                   moderated = TRUE)

View(final.compar)

save(final.compar, file = "final.compar.csv")

write.csv(final.compar, file = "final.compar.csv")

#2021.06.21 Pons

library(MSstatsTMT)

evidence <- read.table("evidence.txt",sep = '\t',header = TRUE)

unique(evidence$Raw.file)

proteinGroups<-read.table("proteinGroups.txt",sep = '\t',header = TRUE)

annotation.mq<-read.csv("annotation.csv")

input.mq<-MaxQtoMSstatsTMTFormat(evidence, proteinGroups,annotation = annotation.mq)

#msstats means medianpolish+impute, no internal normalize channel, 
#so reference norm=false. use this 06.21.2021
protein.summa<- proteinSummarization(data = input.mq, method = "msstats", global_norm = TRUE,
                                     reference_norm = FALSE)

head(protein.summa)

write.csv(protein.summa, file = "protein.summa.csv")

unique(input.mq$Condition)

final.compar <- groupComparisonTMT(data = protein.summa,
                                   contrast.matrix = "pairwise",
                                   moderated = TRUE)

View(final.compar)

save(final.compar, file = "final.compar.csv")

write.csv(final.compar, file = "final.compar.csv")

#2021.07.09 occipital

library(MSstatsTMT)

evidence <- read.table("evidence_occipital.txt",sep = '\t',header = TRUE)

unique(evidence$Raw.file)

proteinGroups<-read.table("proteinGroups_occipital.txt",sep = '\t',header = TRUE)

annotation.mq<-read.csv("annotation.csv")

input.mq<-MaxQtoMSstatsTMTFormat(evidence, proteinGroups,annotation = annotation.mq)

#msstats means medianpolish+impute, no internal normalize channel, 
#so reference norm=false. use this 07.09.2021
protein.summa<- proteinSummarization(data = input.mq, method = "msstats", global_norm = TRUE,
                                     reference_norm = FALSE)

head(protein.summa)

write.csv(protein.summa, file = "protein.summa.csv")

unique(input.mq$Condition)

final.compar <- groupComparisonTMT(data = protein.summa,
                                   contrast.matrix = "pairwise",
                                   moderated = TRUE)

View(final.compar)

save(final.compar, file = "final.compar.csv")

write.csv(final.compar, file = "final.compar.csv")



#2021.07.09 thalamus

library(MSstatsTMT)

evidence <- read.table("evidence_thalamus.txt",sep = '\t',header = TRUE)

unique(evidence$Raw.file)

proteinGroups<-read.table("proteinGroups_thalamus.txt",sep = '\t',header = TRUE)

annotation.mq<-read.csv("annotation.csv")

input.mq<-MaxQtoMSstatsTMTFormat(evidence, proteinGroups,annotation = annotation.mq)

#msstats means medianpolish+impute, no internal normalize channel, 
#so reference norm=false. use this 07.09.2021
protein.summa<- proteinSummarization(data = input.mq, method = "msstats", global_norm = TRUE,
                                     reference_norm = FALSE)

head(protein.summa)

write.csv(protein.summa, file = "protein.summa.csv")

unique(input.mq$Condition)

final.compar <- groupComparisonTMT(data = protein.summa,
                                   contrast.matrix = "pairwise",
                                   moderated = TRUE)

View(final.compar)

save(final.compar, file = "final.compar.csv")

write.csv(final.compar, file = "final.compar.csv")



#2021.07.09 temporal

library(MSstatsTMT)

evidence <- read.table("evidence_temporal.txt",sep = '\t',header = TRUE)

unique(evidence$Raw.file)

proteinGroups<-read.table("proteinGroups_temporal.txt",sep = '\t',header = TRUE)

annotation.mq<-read.csv("annotation.csv")

input.mq<-MaxQtoMSstatsTMTFormat(evidence, proteinGroups,annotation = annotation.mq)

#msstats means medianpolish+impute, no internal normalize channel, 
#so reference norm=false. use this 08.10.2021
protein.summa<- proteinSummarization(data = input.mq, method = "msstats", global_norm = TRUE,
                                     reference_norm = FALSE)

head(protein.summa)

write.csv(protein.summa, file = "protein.summa.csv")

unique(input.mq$Condition)

final.compar <- groupComparisonTMT(data = protein.summa,
                                   contrast.matrix = "pairwise",
                                   moderated = TRUE)

View(final.compar)

save(final.compar, file = "final.compar.csv")

write.csv(final.compar, file = "final.compar.csv")
