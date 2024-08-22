library(ggplot2)
library(ggpubr)
library(readxl)
setwd("C:/Users/tingz/Downloads/COVID_Revision")
Boxplot_input <- read_excel("C:/Users/tingz/Downloads/COVID_Revision/For violinplot area measurement.xlsx")

Boxplot_input <- For_violinplot_area_measurement

Boxplot_input$Group <- as.factor(Boxplot_input$Group)

ACE<- subset(Boxplot_input, Antibody == "ACE")
tiff("ACE_2.tiff", units="in", width= 2.5, height=3, res=300)
g <- ggplot(ACE, aes(x=Group, y=Area, fill = Group)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(legend.position = "none")
g
dev.off()

Abeta <- subset(Boxplot_input, Antibody == "Abeta")
tiff("Abeta.tiff", units="in", width= 2.5, height=3, res=300)
g <- ggplot(Abeta, aes(x=Group, y=Area, fill = Group)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(legend.position = "none")
g
dev.off()

NRP <- subset(Boxplot_input, Antibody == "NRP1")
tiff("NRP.tiff", units="in", width= 2.5, height=3, res=300)
g <- ggplot(NRP, aes(x=Group, y=Area, fill = Group)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(legend.position = "none")
g
dev.off()

Furin<- subset(Boxplot_input, Antibody == "Furin")
tiff("Furin.tiff", units="in", width= 2.5, height=3, res=300)
g <- ggplot(Furin, aes(x=Group, y=Area, fill = Group)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(legend.position = "none")
g
dev.off()

Tau <- subset(Boxplot_input, Antibody == "Tau")
tiff("Tau.tiff", units="in", width= 2.5, height=3, res=300)
g <- ggplot(Tau, aes(x=Group, y=Area, fill = Group)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(legend.position = "none")
g
dev.off()

pdf("Tau.pdf", height = 2.5, width = 3)
a <- ggplot(Tau, aes(x=Group, y=Area, fill = Group)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1))
a
dev.off()

Fibrin <- subset(Boxplot_input, Antibody == "Fibrin")
tiff("Fibrin.tiff", units="in", width= 2.5, height=3, res=300)
g <- ggplot(Fibrin, aes(x=Group, y=Area, fill = Group)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme(legend.position = "none")
g+stat_compare_means(method = "t.test", label = "p", method.args = list(alternative = "less"), label.y = 0.1)
g
dev.off()

nCounter_viral_probes <- Specific_3_probes_

tiff("Extended Fig 1C.tiff", units="in", width= 6, height= 4, res=300)
g <- ggplot(nCounter_viral_probes, aes(x=Targets, y=Counts, fill = Disease)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1)) +
  scale_fill_manual(values=c("#999999","#E69F00"))+
  theme(legend.position = "top")
g+stat_compare_means(method = "t.test", label = "p", label.y = 62)
dev.off()

save(nCounter_viral_probes, Boxplot_input, file="ExtendedFigs.rda")
