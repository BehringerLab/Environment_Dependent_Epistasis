###Figure_3

setwd("/Users/megangrace/Library/CloudStorage/Box-Box/Behringer_Lab_Box_Drive/Manuscripts/In_Progress/Antimicrobial_Resistance/Manuscript/GitHub/Figure_3")


###Figure 3A Plot OmpF and AcrRAB Mutations

#BiocManager::install("trackViewer")
library(trackViewer)
library(cowplot)

##Draw acrRAB and mutations

features <- GRanges("chr1", IRanges(c(481254,484426,485761), 
                                    width=c(3150,1194,648),
                                    names=c("acrB","acrA","acrR")),
                    fill = c("#AA0000","#AA0000","#888888"),
                    height = c(0.2))
SNP <- c(485710,485630,485406,485341,485037,484367,484234,484228,484157,483406,483376,483368,483367,483352,483343,483332,483206,483183,483162,482834,482771,482720,482362,482035,481928,481927,481618,481615,481504)
SNPName<-c("Mut_1","Mut_2","Mut_3","Mut_4","Mut_5","Mut_6","Mut_7","Mut_8","Mut_9","Mut_10","Mut_11","Mut_12","Mut_13","Mut_14","Mut_15","Mut_16","Mut_17","Mut_18","Mut_19","Mut_20","Mut_21","Mut_22","Mut_23","Mut_24","Mut_25","Mut_26","Mut_27","Mut_28","Mut_29")
SNPCount<-c(1,1,1,3,1,2,1,1,1,1,1,1,1,1,1,1,3,1,1,1,2,1,1,1,1,2,1,1,2)
SNPType<-c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey")
Populations <- GRanges("chr1", IRanges(SNP, width=1, names=(SNPName)),
                       color = (SNPType),
                       score = SNPCount)
ranges <- GRanges("chr1", IRanges(481000, 487000))


Populations$label.parameter.gp <- gpar(fontsize=10)

ranges$label.parameter.gp <- gpar(fontsize=10)

acrRAB_locus<-lolliplot(Populations, features,ranges,ylab="", yaxis=FALSE, xaxis.gp = gpar(fontsize = 10))

#Draw ompF and Mutations

ompF_features <- GRanges("chr1", IRanges(c(985894), 
                                         width=c(1089),
                                         names=c("ompF")),
                         fill = c("#0000AA"),
                         height = c(0.2))
ompF_SNP <- c(987098,986580,986579,986567,986549,986550)
ompF_SNPName<-c("Mut_1","Mut_2","Mut_3","Mut_4","Mut_5","Mut_6")
ompF_SNPCount<-c(2,12,4,16,1,1)
ompF_SNPType<-c("grey","grey","grey","grey","grey","grey")
ompF_Populations <- GRanges("chr1", IRanges(ompF_SNP, width=1, names=(ompF_SNPName)),
                            color = (ompF_SNPType),
                            score = ompF_SNPCount)
ompF_ranges <- GRanges("chr1", IRanges(985550, 987350))


ompF_Populations$label.parameter.gp <- gpar(fontsize=10)

ompF_ranges$label.parameter.gp <- gpar(fontsize=10)

ompF_locus<-lolliplot(ompF_Populations, ompF_features,ompF_ranges,ylab="", yaxis=FALSE, xaxis.gp = gpar(fontsize = 10))


#Figure 3B Plot Selection Rate for Reconstructed Mutants. 

library(readr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(cowplot)
library(dbplyr)
library(scales)

acrB_ompF_comp_data<-read.table("Mutant_Competition_Data.txt", sep="\t", header=TRUE)
acrB_ompF_comp_data

acrB_ompF_comp_data$Mutant <- factor(acrB_ompF_comp_data$Mutant, levels = c("ompF", "acrB", "ompF/acrB"))


ggplot(data=acrB_ompF_comp_data, aes(x=Mutant, y=s, fill=Mutant))+
  stat_summary(fun=mean, geom="bar", color="black")+
  stat_summary(fun.data = "mean_se", geom = "errorbar", color="black")+
  scale_fill_manual(values=c("#0F99AA","#F68B26","#aaaacd"))+
  geom_point()+
  ylab("Selection Rate (s)")+
  theme_bw()+
  facet_wrap(~DayPair)

#### Figure 3C 

acrB_ompF_mut_order<-read.table("Mutant_Order_Data.txt", sep="\t",header=TRUE)
acrB_ompF_mut_order

acrB_ompF_mut_order$Order <- factor(acrB_ompF_mut_order$Order, levels = c("OmpF First", "OmpF/AcrAB same time", "AcrAB First","AcrAB Only","None"))


ggplot(data=acrB_ompF_mut_order, aes(x=Order, fill=Order))+
  geom_bar(stat="count", color="black")+
  scale_fill_manual(values=c("#0F99AA","#aaaacd","#aaaacd","#F68B26","#aaaaaa"))+
  facet_wrap(Treatment~.,nrow=3)+
  theme_bw()

#Figure 3D AcrB and OmpF Antimicrobial Susceptibility

acrB_ompF_ABS<-read.table("Mutants_10_day_AMR.txt",sep="\t",header=TRUE)

acrB_ompF_ABS

acrB_ompF_ABS$Background <- factor(acrB_ompF_ABS$Background, levels = c("Ancestor", "P403", "P408","ompF","AcrB","AcrB/OmpF"))


ggplot(data=acrB_ompF_ABS, aes(x=Background, y=Difference, fill=Background))+
  stat_summary(fun=mean, geom="bar", color="black")+
  stat_summary(fun.data = "mean_se", geom = "errorbar", color="black", width=0.5)+
  scale_fill_manual(values=c("#888888","#0000cc","#0000cc","#0F99AA","#F68B26","#aaaacd"))+
  facet_wrap(~Drug, nrow=1)+
  ylab("Radius of Clearence (mm)")+
  theme_bw()+
  theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle=45,hjust=1))

  



