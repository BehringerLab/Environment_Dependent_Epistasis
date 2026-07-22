####Figure_4 and Supplemental Figures S3 and S4. 

##Figure_4A Competition Assays for Reconstruced Mutants in 1-day Spent Media

library(ggplot2)

setwd("/Users/megangrace/Library/CloudStorage/Box-Box/Behringer_Lab_Box_Drive/Manuscripts/In_Progress/Antimicrobial_Resistance/Manuscript/GitHub/Figure_4")

OmpF_AcrB_Spent_Fit<-read.table("OmpF_AcrB_SpentMedia_Comps.txt", sep="\t",header=TRUE)
OmpF_AcrB_Spent_Fit 

OmpF_AcrB_Spent_Fit$Comp2<-factor(OmpF_AcrB_Spent_Fit$Comp2,levels=c("AcrB","OmpF","DoubleMut"))


ggplot(data=OmpF_AcrB_Spent_Fit, aes(x=Comp2,y=s,col=Comp2,fill=Comp2))+
  #geom_point(alpha=1)+
  stat_summary(fun.y=mean, geom="bar", color="black")+
  stat_summary(fun.data="mean_se",geom="errorbar", size=1, width=0.2,color="black")+
  scale_fill_manual(values=c("#f68b26","#0999aa","#005700"))+
  theme_bw()

data.aov<-aov(data=OmpF_AcrB_Spent_Fit, s~Comp2)
TukeyHSD(data.aov)

OmpF_AcrB_Spent_Fit %>% group_by(Comp2) %>% summarise(s.mean=mean(s))

####Figure 4B, S3, and S4 Carbon Utilization Assays

library(ggplot2)
library(growthcurver)
library(dplyr)
library(lubridate)
library(readr)
library(ggpubr)
library(reshape2)
library(cowplot)


#Read in data from each biolog plate replicate

setwd("/Users/megangrace/Library/CloudStorage/Box-Box/Behringer_Lab_Box_Drive/Manuscripts/In_Progress/Antimicrobial_Resistance/Manuscript/GitHub/Figure_4/Biolog_Data")

WT1GC <- 'Biolog_WT_1.csv'
WT2GC <- 'Biolog_WT_2.csv'
WT3GC <- 'Biolog_WT_3.csv'
WT4GC <- 'Biolog_WT_4.csv'
OmpF1GC <- 'Biolog_OMPF_1.csv'
OmpF2GC <- 'Biolog_OMPF_2.csv'
OmpF3GC <- 'Biolog_OMPF_3.csv'

WT1 <- read_csv(WT1GC)
WT2 <- read_csv(WT2GC)
WT3 <- read_csv(WT3GC)
WT4 <- read_csv(WT4GC)
OmpF1 <- read_csv(OmpF1GC)
OmpF2 <- read_csv(OmpF2GC)
OmpF3 <- read_csv(OmpF3GC)

#Convert all Biolog data to long form (melt)
melt_data <- function(data) {
  data |>
    melt(id = c("Time", "T° Read 2:600")) |>
    setNames(c("Time", "T° Read 2:600", "Well", "OD600"))
}

WT1melt <- melt_data(WT1)
WT2melt <- melt_data(WT2)
WT3melt <- melt_data(WT3)
WT4melt <- melt_data(WT4)
OmpF1melt <- melt_data(OmpF1)
OmpF2melt <- melt_data(OmpF2)
OmpF3melt <- melt_data(OmpF3)

#Combine all Biolog replicates into one dataset for plotting growth curves
all_data <- rbind(
  WT1melt |> mutate(strain = "WT", replicate = 1),
  WT2melt |> mutate(strain = "WT", replicate = 2),
  WT3melt |> mutate(strain = "WT", replicate = 3),
  WT4melt |> mutate(strain = "WT", replicate = 4),
  OmpF1melt |> mutate(strain = "OmpF", replicate = 1),
  OmpF2melt |> mutate(strain = "OmpF", replicate = 2), 
  OmpF3melt |> mutate(strain = "OmpF", replicate = 3)
)

all_data

#Changing time format from hh:mm:ss to seconds for plotting growth curves

all_data <- all_data |>
  mutate(Time = hms::as_hms(Time)) |>
  group_by(Well) |>
  mutate(Time_elapsed = as.numeric(Time - first(Time))) |>
  ungroup()

all_data

### Remove temperature column from replicate biolog plate data for calculating growth curve parameters

WT1	<-	WT1[-2]
WT2	<-	WT2[-2]
WT3	<-	WT3[-2]
WT4	<-	WT4[-2]
OmpF1	<-	OmpF1[-2]
OmpF2	<-	OmpF2[-2]
OmpF3	<-	OmpF3[-2]

##Format time column into minutes for estimating growth curve parametes
WT1$Time <- as.numeric(lubridate::hms(WT1$Time), "minutes")
WT2$Time <- as.numeric(lubridate::hms(WT2$Time), "minutes")
WT3$Time <- as.numeric(lubridate::hms(WT3$Time), "minutes")
WT4$Time <- as.numeric(lubridate::hms(WT4$Time), "minutes")
OmpF1$Time <- as.numeric(lubridate::hms(OmpF1$Time), "minutes")
OmpF2$Time <- as.numeric(lubridate::hms(OmpF2$Time), "minutes")
OmpF3$Time <- as.numeric(lubridate::hms(OmpF3$Time), "minutes")

# Trim data to first 15 hours(901 minutes)
WT1	<-	WT1[which(WT1$Time<901),]
WT2	<-	WT2[which(WT2$Time<901),]
WT3	<-	WT3[which(WT3$Time<901),]
WT4	<-	WT4[which(WT4$Time<901),]
OmpF1	<-	OmpF1[which(OmpF1$Time<901),]
OmpF2	<-	OmpF2[which(OmpF2$Time<901),]
OmpF3	<-	OmpF3[which(OmpF3$Time<901),]

#Run growthcurver to estimate growth curve parameters
WT1_values<-SummarizeGrowthByPlate(WT1)
WT2_values<-SummarizeGrowthByPlate(WT2)
WT3_values<-SummarizeGrowthByPlate(WT3)
WT4_values<-SummarizeGrowthByPlate(WT4)
OmpF1_values<-SummarizeGrowthByPlate(OmpF1)
OmpF2_values<-SummarizeGrowthByPlate(OmpF2)
OmpF3_values<-SummarizeGrowthByPlate(OmpF3)

#Add Metadata to Growth Curve Parameter Data
WT1_values$Background <- "WT"
WT2_values$Background <- "WT"
WT3_values$Background <- "WT"
WT4_values$Background <- "WT"
OmpF1_values$Background <- "OmpF"
OmpF2_values$Background <- "OmpF"
OmpF3_values$Background <- "OmpF"

WT1_values$Replicate <- "WT1"
WT2_values$Replicate <- "WT2"
WT3_values$Replicate <- "WT3"
WT4_values$Replicate <- "WT4"
OmpF1_values$Replicate <- "OmpF1"
OmpF2_values$Replicate <- "OmpF2"
OmpF3_values$Replicate <- "OmpF3"

#Combine all growth curve parameter datasets
All_GC_Values<-rbind(WT1_values,WT2_values,WT3_values,WT4_values,OmpF1_values,OmpF2_values,OmpF3_values)
All_GC_Values

#Look at AUCe Values across all resources to determine cutoff for informative growth
ggplot(All_GC_Values,aes(x=sample,y=auc_e,col=Background))+
  stat_summary(fun.x=mean, geom="point")+
  stat_summary(fun.data="mean_se",geom="errorbar", size=1, width=0.2)+
  theme(axis.text.x=element_text(size=8,angle=90,hjust=1))

#Calculate Max AUCe value for each resource and create list to remove resources where Max AUCe is <250
Max_Values<-All_GC_Values %>% group_by(sample) %>% summarise(max.aucE=max(auc_e))

print(Max_Values[which(Max_Values$max.aucE<250),],n=50)

#Remove resources where max AUCe <250
All_GC_Values_clean<-All_GC_Values[which(All_GC_Values$sample!="A1"&All_GC_Values$sample!="A8"&All_GC_Values$sample!="B10"&All_GC_Values$sample!="B12"&All_GC_Values$sample!="C5"&All_GC_Values$sample!="D1"&All_GC_Values$sample!="D10"&All_GC_Values$sample!="D11"&All_GC_Values$sample!="D2"&All_GC_Values$sample!="D3"&All_GC_Values$sample!="D4"&All_GC_Values$sample!="D5"&All_GC_Values$sample!="E1"&All_GC_Values$sample!="E5"&All_GC_Values$sample!="E6"&All_GC_Values$sample!="E9"&All_GC_Values$sample!="F11"&All_GC_Values$sample!="F2"&All_GC_Values$sample!="F3"&All_GC_Values$sample!="F4"&All_GC_Values$sample!="F7"&All_GC_Values$sample!="G2"&All_GC_Values$sample!="G7"&All_GC_Values$sample!="H11"&All_GC_Values$sample!="H12"&All_GC_Values$sample!="H2"&All_GC_Values$sample!="H3"&All_GC_Values$sample!="H4"&All_GC_Values$sample!="H5"&All_GC_Values$sample!="H6"&All_GC_Values$sample!="H7"),]

#Identify resources where OmpF has significant difference in AUCe

Compare_Samples_AUC<-compare_means(data=All_GC_Values_clean, auc_e~Background, group.by = "sample",method="t.test")
print(Compare_Samples_AUC, n=80)

my_comparisons <- list( c("WT", "OmpF"))

##Plot Supplemental Figure S3 - AUC across resources
ggplot(All_GC_Values_clean,aes(x=Background,y=auc_e,col=Background))+
  geom_jitter(alpha=0.5)+
  stat_summary(fun.x=mean, geom="point")+
  stat_summary(fun.data="mean_se",geom="errorbar", size=1, width=0.2)+
  stat_compare_means(comparisons = my_comparisons,method="t.test")+
  scale_y_continuous(limits=c(0,1500))+
  scale_color_manual(values=c("#0999AA","#888888"))+
  theme_bw()+
  theme(legend.position="top",axis.text.x=element_text(size=8,angle=90,hjust=1))+
  facet_wrap(~sample, nrow=3)

#Plot Supplemental Data 4 - Growth Curves across Resources

Resources<-c("A10","A11","A12","A2","A3","A4","A5","A6","A7","A9","B1","B11","B2","B3","B4","B5","B6","B7","B8","B9","C1","C10","C11","C12","C2","C3","C4","C6","C7","C8","C9","D12","D6","D7","D8","D9","E10","E11","E12","E2","E3","E4","E7","E8","F1","F10","F12","F5","F6","F8","F9","G1","G10","G11","G12","G3","G4","G5","G6","G8","G9","H1","H10","H8","H9")

all_data |>
  filter(Well %in% Resources, Time_elapsed <= 54000) |>
  group_by(replicate, strain, Time_elapsed, Well) |>
  summarise(mean_OD600 = mean(OD600), .groups = 'drop') |>
  ggplot(aes(x = Time_elapsed/3600, y = mean_OD600, col = strain,)) +
  stat_summary(fun=mean,geom="line")+
  ylim(0, 2) +
  scale_color_manual(values=c("#0999AA","#888888"))+
  theme_bw() +
  theme(legend.position = "top") +
  #geom_hline(yintercept = 0.25, linetype = "dashed", color = "gray") +
  labs(color = "Strain", y = "Mean OD600", x = "Time (hours)")+
  facet_wrap(~Well,nrow=3,scales="free_x")

#Plot Figure 4B - Growth Curves of Selected Representative Resources

Rep_Resources<-c("A5","A7","C3","B9","D6","F5","F9","C12","H1","G3","G5","F12")

all_data |>
  filter(Well %in% Rep_Resources, Time_elapsed <= 54000) |>
  group_by(replicate, strain, Time_elapsed, Well) |>
  summarise(mean_OD600 = mean(OD600), .groups = 'drop') |>
  ggplot(aes(x = Time_elapsed/3600, y = mean_OD600, col = strain,)) +
  stat_summary(fun=mean,geom="line")+
  ylim(0, 2) +
  scale_color_manual(values=c("#0999AA","#888888"))+
  theme_bw() +
  theme(legend.position = "top") +
  #geom_hline(yintercept = 0.25, linetype = "dashed", color = "gray") +
  labs(color = "Strain", y = "Mean OD600", x = "Time (hours)")+
  facet_wrap(~Well,nrow=3,scales="free_x")





