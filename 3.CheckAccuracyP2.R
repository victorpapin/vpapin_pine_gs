library(ggplot2)
library(cowplot)

RESULTATS.real  = read.table("RESULTATS.REAL.txt", header=TRUE, sep="")
RESULTATS.simul = read.table("RESULTATS.SIMUL.txt", header=TRUE, sep="")

RESULTATS.simul2          = RESULTATS.simul
RESULTATS.simul2$DATA     = "simulGV"
RESULTATS.simul2$ACCURACY = RESULTATS.simul2$COR2

RESULTATS.real$COR2      = NA
RESULTATS.real$ACCURACY2 = NA
RESULTATS.real      = RESULTATS.real[,c("DATA","TYPE","POURCENT","IT","COR","COR2","h2","ACCURACY","ACCURACY2")]
RESULTATS           = rbind(RESULTATS.real, RESULTATS.simul, RESULTATS.simul2)

RESULTATS[which(RESULTATS$POURCENT==0.2),"POURCENT"] ="Tset_80%"
RESULTATS[which(RESULTATS$POURCENT==0.5),"POURCENT"] ="Tset_50%"

plot1 = ggplot(RESULTATS, aes(x=DATA, y=ACCURACY, group=DATA))+facet_grid(POURCENT ~ TYPE)+
          geom_violin(aes(fill=DATA))+
          geom_boxplot(width=0.1, fill="grey")+ stat_summary(fun=mean, fill="darkred", geom="point", shape=21, size=3, show.legend=FALSE)+
          theme_bw()+theme(legend.position="none")+
          xlab("")+ylab("GS accuracy")+ylim(0,1)

plot2 = ggplot(RESULTATS, aes(x=DATA, y=ACCURACY, fill=DATA))+facet_grid(POURCENT~TYPE)+
  geom_boxplot(width=0.2)+ stat_summary(fun=mean, fill="darkred", geom="point", shape=21, size=3, show.legend=FALSE)+
  theme_bw()+theme(legend.position="none")+
  xlab("")+ylab("GS accuracy")+ylim(0,1)


plot3 = ggplot(RESULTATS, aes(x=TYPE, y=ACCURACY, fill=DATA))+facet_grid(POURCENT~.)+
  geom_boxplot(width=0.2)+
  theme_bw()+
  xlab("")+ylab("GS accuracy")+ylim(0,1)


plot1
plot2
plot3

ggsave("CheckAccuracyP2_1.tiff", plot1, width= 20, height= 20, units= "cm", scale = 1, dpi= 600)
ggsave("CheckAccuracyP2_2.tiff", plot2, width= 20, height= 20, units= "cm", scale = 1, dpi= 600)
ggsave("CheckAccuracyP2_3.tiff", plot3, width= 20, height= 20, units= "cm", scale = 1, dpi= 600)


####################################################################

RESULTATS.FAM.real  = read.table("RESULTATS.REAL.FAM.txt", header=TRUE, sep="")
RESULTATS.FAM.simul = read.table("RESULTATS.SIMUL.FAM.txt", header=TRUE, sep="")

RESULTATS.FAM.simul2      = RESULTATS.FAM.simul
RESULTATS.FAM.simul2$DATA = "simulGV"
RESULTATS.FAM.simul2$COR  = RESULTATS.FAM.simul2$COR2

RESULTATS.FAM.real$COR2      = NA
RESULTATS.FAM.real      = RESULTATS.FAM.real[,c("DATA","TYPE","POURCENT","FAM","IT","COR","COR2","NB")]
RESULTATS.FAM           = rbind(RESULTATS.FAM.real, RESULTATS.FAM.simul, RESULTATS.FAM.simul2)

RESULTATS.FAM = RESULTATS.FAM[RESULTATS.FAM$POURCENT==0.5,]
RESULTATS.FAM[which(RESULTATS.FAM$POURCENT==0.5),"POURCENT"] ="Tset_50%"

RESULTATS.FAM$FAM = as.factor(RESULTATS.FAM$FAM)
RESULTATS.FAM     = RESULTATS.FAM[RESULTATS.FAM$TYPE=="genomic",]
gd.fam = as.character(unique(RESULTATS.FAM[RESULTATS.FAM$NB>10,"FAM"]))

RESULTATS.FAM.tmp = RESULTATS.FAM[RESULTATS.FAM$FAM %in% gd.fam,]

plot4 = ggplot(RESULTATS.FAM.tmp, aes(x=reorder(FAM,COR, FUN=mean), y=COR, fill=DATA))+
  geom_boxplot()+ facet_grid(.~POURCENT)+
  theme(axis.text.x = element_text(angle=90))+
  xlab("families")+ylab("predictive ability")
plot4

plot5a = ggplot(RESULTATS.FAM.tmp[RESULTATS.FAM.tmp$DATA=="real",], aes(x=reorder(FAM,COR, FUN=mean), y=COR, fill=DATA))+geom_boxplot()+ theme(axis.text.x = element_text(angle=90))+xlab("families")+ylab("predictive ability")+theme(legend.position="none")+ facet_grid(.~POURCENT)
plot5b = ggplot(RESULTATS.FAM.tmp[RESULTATS.FAM.tmp$DATA=="simul",], aes(x=reorder(FAM,COR, FUN=mean), y=COR, fill=DATA))+geom_boxplot(fill="chartreuse3")+ theme(axis.text.x = element_text(angle=90))+xlab("families")+ylab("predictive ability")+ facet_grid(.~POURCENT)
plot5c = ggplot(RESULTATS.FAM.tmp[RESULTATS.FAM.tmp$DATA=="simulGV",], aes(x=reorder(FAM,COR, FUN=mean), y=COR, fill=DATA))+geom_boxplot(fill="cadetblue3")+ theme(axis.text.x = element_text(angle=90))+xlab("families")+ylab("predictive ability")+ facet_grid(.~POURCENT)
plot5 = plot_grid(plot5a, plot5b, plot5c, nrow=1)
plot5

plot6a = ggplot(RESULTATS.FAM[RESULTATS.FAM$DATA=="real",], aes(x=reorder(FAM,COR, FUN=mean), y=COR, fill=DATA))+geom_boxplot()+ theme(axis.text.x = element_text(angle=90))+xlab("families")+ylab("predictive ability")+theme(legend.position="none")+ facet_grid(.~POURCENT)
plot6b = ggplot(RESULTATS.FAM[RESULTATS.FAM$DATA=="simul",], aes(x=reorder(FAM,COR, FUN=mean), y=COR, fill=DATA))+geom_boxplot(fill="chartreuse3")+ theme(axis.text.x = element_text(angle=90))+xlab("families")+ylab("predictive ability")+ facet_grid(.~POURCENT)
plot6c = ggplot(RESULTATS.FAM[RESULTATS.FAM$DATA=="simulGV",], aes(x=reorder(FAM,COR, FUN=mean), y=COR, fill=DATA))+geom_boxplot(fill="cadetblue3")+ theme(axis.text.x = element_text(angle=90))+xlab("families")+ylab("predictive ability")+ facet_grid(.~POURCENT)
plot6 = plot_grid(plot6a, plot6b, plot6c, nrow=1)
plot6

ggsave("CheckAccuracyP2_intrafam_1.tiff", plot4, width= 20, height= 20, units= "cm", scale = 1, dpi= 600)
ggsave("CheckAccuracyP2_intrafam_2.tiff", plot5, width= 40, height= 20, units= "cm", scale = 1, dpi= 600)
ggsave("CheckAccuracyP2_intrafam_3.tiff", plot6, width= 60, height= 20, units= "cm", scale = 1, dpi= 600)

