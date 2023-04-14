######################################
#
#BLOOD SMEAR GRAPHING 
#
# Created by KW October 2022 
# modified by K.Weisgerber on 4/13/23
#
############################################
#load libraries
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(emmeans)
library(rstatix)
library(performance)
library(car)

setwd("/Users/kaitl/Desktop/LS_ploidy/raw_data/")

#this is the master data file, will probably have to subset 
Master<-read.csv("comb_ploidy_data.csv")

Cells<-Master[Master$Celltype=="avg_cell", c("Smear_Pop", "Smr_samp", "Diam_um", "avg_volume_um3","avg_volume_uL", "avg_surface_um2","SAVratio")]
#okay "cells" has all the data for the cell type 

#lets make a file for nuclei while were at it 
Nucleis<-Master[Master$Celltype=="avg_nuclei", c("Smear_Pop", "Smr_samp", "Diam_um", "avg_volume_um3","avg_volume_uL", "avg_surface_um2","SAVratio")]

##remove the TEN, that was a standard not a sturgeon 

Cells1<-Cells[!(Cells$Smear_Pop=="TEN"),]

Nuc<-Nucleis[!(Nucleis$Smear_Pop=="TEN"),]
#############################################################

#starting with the cell volume 

##making a violoin plot with mean as a crossbar 
cellv<-Cells1 %>% 
  filter(!is.na(Smear_Pop))%>% 
  mutate(CC_Pop = fct_relevel(Smear_Pop,"BR", "LR",
                              "NR","PDB", "UMFeb", "UMMay"))%>%
  ggplot(five, mapping = aes(x = Smear_Pop, y= avg_volume_um3)) +
  geom_violin(lwd = 0.5)+ 
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.5,
               color = "grey")+
  geom_point(shape = 16) + 
  theme_classic() +
  ylab(bquote('Average Cell Volume'~ ~ (um^3))) +
  xlab('Sample Population')+ 
  theme(text = element_text(size = 15))
cellv

#making the axis labels full population names 
CV<- cellv + scale_x_discrete("Sample Population", labels = c("BR" = "Birthday Rapids", "LR" = "Landing River", "NR" = "Nelson River", "PDB" = "Pointe du Bois", "UMFeb" = "UM 268", "UMMay" = "UM 344" ))
CV

#add the statistical significance 
CV_stat<- CV + ylim(4000, 13500)+
  annotate("text",
                         x = 1, y = 10000,
                         label = "ac",
                         size = 6)+
  annotate("text",
           x = 2, y = 10500,
           label = "ac",
           size = 6)+
  annotate("text",
           x = 3, y = 8500,
           label = "b",
           size = 6)+
  annotate("text",
           x = 4, y = 13500,
           label = "ab",
           size = 6)+
  annotate("text",
           x = 5, y = 10000,
           label = "ac",
           size = 6)+
  annotate("text",
           x = 6, y = 13000,
           label = "ac",
           size = 6) 
CV_stat
ggsave(filename = "Cellvolume.tiff", plot = CV_stat, dpi = 600, height = 10, width = 12)


##################################

#cell surface area 
cellsa<-Cells1 %>% 
  filter(!is.na(Smear_Pop))%>% 
  mutate(CC_Pop = fct_relevel(Smear_Pop,"BR", "LR",
                              "NR","PDB", "UMFeb", "UMMay"))%>%
  ggplot(five, mapping = aes(x = Smear_Pop, y= avg_surface_um2)) +
  geom_violin(lwd = 0.5)+ 
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.5,
               color = "grey")+
  geom_point(shape = 16) + 
  theme_classic() +
  ylab(bquote('Average Cell Surface Area'~ ~ (um^2))) +
  xlab('Sample Population')+ 
  theme(text = element_text(size = 15))
cellsa

#making the axis labels full population names 
CSA<- cellsa + scale_x_discrete("Sample Population", labels = c("BR" = "Birthday Rapids", "LR" = "Landing River", "NR" = "Nelson River", "PDB" = "Pointe du Bois", "UMFeb" = "UM 268", "UMMay" = "UM 344" ))
CSA

#add the statistical significance 
CSA_stat<- CSA + ylim(350, 800)+
  annotate("text",
           x = 1, y = 630,
           label = "a",
           size = 6)+
  annotate("text",
           x = 2, y = 650,
           label = "ab",
           size = 6)+
  annotate("text",
           x = 3, y = 575,
           label = "c",
           size = 6)+
  annotate("text",
           x = 4, y = 800,
           label = "abc",
           size = 6)+
  annotate("text",
           x = 5, y = 660,
           label = "ab",
           size = 6)+
  annotate("text",
           x = 6, y = 800,
           label = "a",
           size = 6) 
CSA_stat
ggsave(filename = "Cellsurface.tiff", plot = CSA_stat, dpi = 600, height = 10, width = 12)

#########################################

##cell surface area to volume ratio 

cellsav<-Cells1 %>% 
  filter(!is.na(Smear_Pop))%>% 
  mutate(CC_Pop = fct_relevel(Smear_Pop,"BR", "LR",
                              "NR","PDB", "UMFeb", "UMMay"))%>%
  ggplot(five, mapping = aes(x = Smear_Pop, y= SAVratio)) +
  geom_violin(lwd = 0.5)+ 
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.5,
               color = "grey")+
  geom_point(shape = 16) + 
  theme_classic() +
  ylab(bquote('Average Cell SA:V Ratio'~ ~ (um))) +
  xlab('Sample Population')+ 
  theme(text = element_text(size = 15))
cellsav

#making the axis labels full population names 
CSAV<- cellsav + scale_x_discrete("Sample Population", labels = c("BR" = "Birthday Rapids", "LR" = "Landing River", "NR" = "Nelson River", "PDB" = "Pointe du Bois", "UMFeb" = "UM 268", "UMMay" = "UM 344" ))
CSAV

#add the statistical significance 
CSAV_stat<- CSAV + ylim(0.055, 0.09)+
  annotate("text",
           x = 1, y = 0.077,
           label = "a",
           size = 6)+
  annotate("text",
           x = 2, y = 0.08,
           label = "ab",
           size = 6)+
  annotate("text",
           x = 3, y = 0.089,
           label = "c",
           size = 6)+
  annotate("text",
           x = 4, y = 0.087,
           label = "ac",
           size = 6)+
  annotate("text",
           x = 5, y = 0.082,
           label = "ab",
           size = 6)+
  annotate("text",
           x = 6, y = 0.089,
           label = "ab",
           size = 6) 
CSAV_stat
ggsave(filename = "CellSAV.tiff", plot = CSAV_stat, dpi = 600, height = 10, width = 12)


#############################################################

##nuclei volume 

nucleiv<-Nuc %>% 
  filter(!is.na(Smear_Pop))%>% 
  mutate(CC_Pop = fct_relevel(Smear_Pop,"BR", "LR",
                              "NR","PDB", "UMFeb", "UMMay"))%>%
  ggplot(five, mapping = aes(x = Smear_Pop, y= avg_volume_um3)) +
  geom_violin(lwd = 0.5)+ 
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.5,
               color = "grey")+
  geom_point(shape = 16) + 
  theme_classic() +
  ylab(bquote('Average Nuclei Volume'~ ~ (um^3))) +
  xlab('Sample Population')+ 
  theme(text = element_text(size = 15))
nucleiv

#making the axis labels full population names 
NV<- nucleiv + scale_x_discrete("Sample Population", labels = c("BR" = "Birthday Rapids", "LR" = "Landing River", "NR" = "Nelson River", "PDB" = "Pointe du Bois", "UMFeb" = "UM 268", "UMMay" = "UM 344" ))
NV

#add the statistical significance 
NV_stat<- NV + ylim(300, 1900)+
  annotate("text",
           x = 1, y = 900,
           label = "a",
           size = 6)+
  annotate("text",
           x = 2, y = 1000,
           label = "ab",
           size = 6)+
  annotate("text",
           x = 3, y = 700,
           label = "c",
           size = 6)+
  annotate("text",
           x = 4, y = 1000,
           label = "bc",
           size = 6)+
  annotate("text",
           x = 5, y = 1900,
           label = "abd",
           size = 6)+
  annotate("text",
           x = 6, y = 1500,
           label = "abd",
           size = 6) 
NV_stat
ggsave(filename = "NucleiV.tiff", plot = NV_stat, dpi = 600, height = 10, width = 12)

#####################################

##nuclei surface area 

nucleisa<-Nuc %>% 
  filter(!is.na(Smear_Pop))%>% 
  mutate(CC_Pop = fct_relevel(Smear_Pop,"BR", "LR",
                              "NR","PDB", "UMFeb", "UMMay"))%>%
  ggplot(five, mapping = aes(x = Smear_Pop, y= avg_surface_um2)) +
  geom_violin(lwd = 0.5)+ 
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.5,
               color = "grey")+
  geom_point(shape = 16) + 
  theme_classic() +
  ylab(bquote('Average Nuclei Surface Area'~ ~ (um^2))) +
  xlab('Sample Population')+ 
  theme(text = element_text(size = 15))
nucleisa

#making the axis labels full population names 
NSA<- nucleisa + scale_x_discrete("Sample Population", labels = c("BR" = "Birthday Rapids", "LR" = "Landing River", "NR" = "Nelson River", "PDB" = "Pointe du Bois", "UMFeb" = "UM 268", "UMMay" = "UM 344" ))
NSA

#add the statistical significance 
NSA_stat<- NSA + ylim(50, 225)+
  annotate("text",
           x = 1, y = 130,
           label = "a",
           size = 6)+
  annotate("text",
           x = 2, y = 140,
           label = "ab",
           size = 6)+
  annotate("text",
           x = 3, y = 120,
           label = "c",
           size = 6)+
  annotate("text",
           x = 4, y = 150,
           label = "c",
           size = 6)+
  annotate("text",
           x = 5, y = 210,
           label = "abd",
           size = 6)+
  annotate("text",
           x = 6, y = 190,
           label = "abd",
           size = 6) 
NSA_stat
ggsave(filename = "NucleiSA.tiff", plot = NSA_stat, dpi = 600, height = 10, width = 12)

##########################################

##Nuclei surface area to volume ratio 

Nucsav<-Nuc %>% 
  filter(!is.na(Smear_Pop))%>% 
  mutate(CC_Pop = fct_relevel(Smear_Pop,"BR", "LR",
                              "NR","PDB", "UMFeb", "UMMay"))%>%
  ggplot(five, mapping = aes(x = Smear_Pop, y= SAVratio)) +
  geom_violin(lwd = 0.5)+ 
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.5,
               color = "grey")+
  geom_point(shape = 16) + 
  theme_classic() +
  ylab(bquote('Average Nuclei SA:V Ratio'~ ~ (um))) +
  xlab('Sample Population')+ 
  theme(text = element_text(size = 15))
Nucsav

#making the axis labels full population names 
NSAV<- Nucsav + scale_x_discrete("Sample Population", labels = c("BR" = "Birthday Rapids", "LR" = "Landing River", "NR" = "Nelson River", "PDB" = "Pointe du Bois", "UMFeb" = "UM 268", "UMMay" = "UM 344" ))
NSAV

#add the statistical significance 
NSAV_stat<- NSAV + ylim(0.1, 0.23)+
  scale_y_continuous()+
  annotate("text",
           x = 1, y = 0.189,
           label = "a",
           size = 6)+
  annotate("text",
           x = 2, y = 0.195,
           label = "ab",
           size = 6)+
  annotate("text",
           x = 3, y = 0.22,
           label = "c",
           size = 6)+
  annotate("text",
           x = 4, y = 0.22,
           label = "bd",
           size = 6)+
  annotate("text",
           x = 5, y = 0.197,
           label = "ab",
           size = 6)+
  annotate("text",
           x = 6, y = 0.2,
           label = "abd",
           size = 6) 
NSAV_stat
ggsave(filename = "NucSAV.tiff", plot = NSAV_stat, dpi = 600, height = 10, width = 12)

###########################################

#combine all 6 into a single plot 

AllSmear<-(CV_stat | NV_stat) /
  (CSA_stat| NSA_stat) /
  (CSAV_stat| NSAV_stat)+
  plot_annotation(tag_levels = 'A')
AllSmear
ggsave(filename = "Combsmears2.tiff", plot = AllSmear, dpi = 600, height = 16, width = 18)
