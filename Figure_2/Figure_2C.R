library("ggplot2")
library("dplyr")
library("ggpubr")

setwd("/Users/megangrace/Library/CloudStorage/Box-Box/Behringer_Lab_Box_Drive/Manuscripts/In_Progress/Antimicrobial_Resistance/Manuscript/GitHub/Figure_2")
AMR_Clones<- read.table("Evolved_Clones_AMR.txt",sep="\t",header=TRUE)
AMR_Clones

AMR_Clones$Background<-factor(AMR_Clones$Background,levels=c("Ancestor","P103","P108","P403","P408","P503","P508"))
AMR_Clones_Ans_only<-AMR_Clones[which(AMR_Clones$Background=="Ancestor"),]

AMR_Clones_Ans_Drug<-AMR_Clones_Ans_only %>% group_by(Drug) %>% summarise(mean.drug=mean(Difference))



ggplot(data=AMR_Clones, aes(x=Background,y=Difference,,group=Background,fill=Background))+
  #geom_point(alpha=0.4)+
  #geom_boxplot(color="black", alpha=0.4)+
  #stat_summary(fun = "mean", geom = "point", size=3, pch=21, color="black")+
  stat_summary(fun = "mean", geom = "bar", color="black", alpha=1)+
  stat_summary(fun.data = "mean_se", geom = "errorbar", width=0.2, color="black")+
  geom_hline(data=AMR_Clones_Ans_Drug, aes(yintercept=mean.drug),color="black",size=0.7, linetype="dashed")+
  #scale_color_manual(values=c("#888888","#E69F00","#E69F00","#56b4e9","#56b4e9","#009e73","#009e73"),name="Strain")+
  scale_fill_manual(values=c("#888888","#DCBC77","#DCBC77","#4e74aa","#4e74aa","#bb6a8e","#bb6a8e"),name="Strain")+
  #xlab("Genotypes")+
  ylab("Radius of Clearence (mm)")+
  scale_x_discrete(breaks=c("Ancestor","P103","P108","P403","P408","P503","P508"),labels=c("WT","1-day Pop A","1-day Pop B","10-day Pop A","10-day Pop B","100-day Pop A","100-day Pop B"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1), panel.grid.minor = element_blank(),panel.grid.major = element_blank(),axis.title.x = element_blank(),legend.position = "none")+
  facet_grid(~Drug)

###Stats
ERY<-AMR_Clones[which(AMR_Clones$Drug=="ERY"),]
NAL<-AMR_Clones[which(AMR_Clones$Drug=="NAL"),]
PEN<-AMR_Clones[which(AMR_Clones$Drug=="PEN"),]
STR<-AMR_Clones[which(AMR_Clones$Drug=="STR"),]
TET<-AMR_Clones[which(AMR_Clones$Drug=="TET"),]
TMS<-AMR_Clones[which(AMR_Clones$Drug=="TMS"),]

compare_means(data=ERY,Difference~Evo_Round_1, ref.group = "Ancestor", method="t.test")
compare_means(data=NAL,Difference~Evo_Round_1, ref.group = "Ancestor", method="t.test")
compare_means(data=PEN,Difference~Evo_Round_1, ref.group = "Ancestor",method="t.test")
compare_means(data=STR,Difference~Evo_Round_1, ref.group = "Ancestor",method="t.test")
compare_means(data=TET,Difference~Evo_Round_1, ref.group = "Ancestor",method="t.test")
compare_means(data=TMS,Difference~Evo_Round_1, ref.group = "Ancestor",method="t.test")


