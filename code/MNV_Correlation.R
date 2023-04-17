######################################
#
#RELATIONSHIP BTWN MNV AND BIOMETRICS 
#
# Created by KW October 2022 
# modified by K.Weisgerber on 4/17/23
#
############################################
#load libraries
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(emmeans)
library(rstatix)
library(performance)
library(ggpubr)

setwd("/Users/kaitl/Desktop/LS_ploidy/raw_data/")

#this is the master data file, will probably have to subset 
Master<-read.csv("comb_ploidy_data.csv")

#this makes it just the data for the five main populations 
main<- Master[c(1:598,865:1020), ]

UMonly<- Master[c(599:976),]

#want to see if there is a correlation between length/ mass and MNV 
#first need to know if this data is normal or not 

######################################
##starting with main 5 population groups 

shapiro.test(main$avg_nuclei)
shapiro.test(main$TLeng_mm)
shapiro.test(main$mass_g)
# these are not showing normal distribution, very small p 

ggqqplot(main$avg_nuclei)
#this doesnt look bad 

ggqqplot(main$TLeng_mm)
#the length is def not normal 

ggqqplot(main$mass_g)
#also not normal 

#could use spearmand or kendall correlation 

lcor<-cor.test(main$avg_nuclei, main$TLeng_mm, method = "spearman")
lcor
#	Spearman's rank correlation rho
#data:  main$avg_nuclei and main$TLeng_mm
#S = 111406419, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#-0.7112369 

mcor<-cor.test(main$avg_nuclei, main$mass_g, method = "spearman")
mcor
#data:  main$avg_nuclei and main$mass_g
#S = 115352899, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#-0.7219141 

##this is telling us that when volume increases, mass and length decrease 
## or that smaller animals have larger nuclei volumes 

############################################
##UM population 

shapiro.test(UMonly$avg_nuclei)
shapiro.test(UMonly$TLeng_mm)
shapiro.test(UMonly$mass_g)
#these are also not normal 

UMlcor<-cor.test(UMonly$avg_nuclei, UMonly
                 $TLeng_mm, method = "spearman")
UMlcor

#data:  UMonly$avg_nuclei and UMonly$TLeng_mm
#S = 10474538, p-value = 1.357e-14
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#-0.392963 

UMmcor<-cor.test(UMonly$avg_nuclei, UMonly$mass_g, method = "spearman")
UMmcor
#data:  UMonly$avg_nuclei and UMonly$mass_g
#S = 11316935, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#-0.4195918 

###for both the main 5 and the  UM, larger volume are in smaller length and mass 
##so as they grow, the MNV decreases 
