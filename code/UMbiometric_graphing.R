######################################
#
#BIOMETRICS GRAPHING BETWEEN UM TIME POINTS   
#
# Created by KW October 2022 
# modified by K.Weisgerber on 4/15/23
#
#code for dates to DPH 
#175 dph = Nov, 200= Dec, 225 = Jan, 240 = Jan_, 268 = Feb, 
## 325 = april, 344 = may 
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
library(patchwork)

##total length 


UMTL<-UMonly %>% 
  filter(!is.na(CC_Pop))%>% 
  mutate(CC_Pop = fct_relevel(CC_Pop,"UMNov_", "UMDec",
                              "UMJan","UMJan_","UMFeb","UMApril","UMMay"))%>%
  ggplot(UMonly, mapping = aes(x = CC_Pop, y= TLeng_mm)) +
  geom_violin(lwd = 0.5)+ 
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.5,
               color = "grey")+
  geom_point(shape = 16) +
  theme_classic() +
  ylab(bquote('Total length'~ ~ (mm))) +
  xlab('Days post hatch')+ 
  theme(text = element_text(size = 15))
UMTL

#relabelling x axis 
UM_TL<-UMTL + scale_x_discrete("Days post hatch", labels = c("UMNov_" = "175", "UMDec" = "200", "UMJan" = "225", "UMJan_" = "240", "UMFeb" = "268", "UMApril" = "325", "UMMay" = "344"))
UM_TL

#adding significance values 
leng_um<- UM_TL + annotate("text",
                         x = 1, y = 160,
                         label = "a",
                         size = 6)+
  annotate("text",
           x = 2, y = 175,
           label = "ab",
           size = 6)+
  annotate("text",
           x = 3, y = 190,
           label = "c",
           size = 6)+
  annotate("text",
           x = 4, y = 235,
           label = "bc",
           size = 6)+
  annotate("text",
           x = 5, y = 235,
           label = "d",
           size = 6)+
  annotate("text",
           x = 6, y = 320,
           label = "d",
           size = 6)+
  annotate("text",
           x = 7, y = 275,
           label = "d",
           size = 6)
leng_um
#save graph 

ggsave(filename = "LengthUM.tiff", plot = leng_um, dpi = 600, height = 10, width = 12)


##############################################

##wet mass 

UMMASS<-UMonly %>% 
  filter(!is.na(CC_Pop))%>% 
  mutate(CC_Pop = fct_relevel(CC_Pop,"UMNov_", "UMDec",
                              "UMJan","UMJan_","UMFeb","UMApril","UMMay"))%>%
  ggplot(UMonly, mapping = aes(x = CC_Pop, y= mass_g)) +
  geom_violin(lwd = 0.5)+ 
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.5,
               color = "grey")+
  geom_point(shape = 16) +
  theme_classic() +
  ylab(bquote('Wet mass'~ ~ (g))) +
  xlab('Days post hatch')+ 
  theme(text = element_text(size = 15))
UMMASS

#relabelling x axis 
UM_MASS<-UMMASS + scale_x_discrete("Days post hatch", labels = c("UMNov_" = "175", "UMDec" = "200", "UMJan" = "225", "UMJan_" = "240", "UMFeb" = "268", "UMApril" = "325", "UMMay" = "344"))
UM_MASS

#adding significance values 
mass_um<- UM_MASS + annotate("text",
                           x = 1, y = 18,
                           label = "a",
                           size = 6)+
  annotate("text",
           x = 2, y = 20,
           label = "a",
           size = 6)+
  annotate("text",
           x = 3, y = 22,
           label = "a",
           size = 6)+
  annotate("text",
           x = 4, y = 37,
           label = "ab",
           size = 6)+
  annotate("text",
           x = 5, y = 35,
           label = "bc",
           size = 6)+
  annotate("text",
           x = 6, y = 47,
           label = "c",
           size = 6)

mass_um
#save graph 

ggsave(filename = "wetmassUM.tiff", plot = mass_um, dpi = 600, height = 10, width = 12)


###################################

##condition factor 
UMK<-data_K1 %>% 
  filter(!is.na(CC_Pop))%>% 
  mutate(CC_Pop = fct_relevel(CC_Pop,"UMNov_", "UMDec",
                              "UMJan","UMJan_","UMFeb","UMApril","UMMay"))%>%
  ggplot(UMonly, mapping = aes(x = CC_Pop, y= K)) +
  geom_violin(lwd = 0.5)+ 
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.5,
               color = "grey")+
  geom_point(shape = 16) +
  theme_classic() +
  ylab(bquote('Condition factor'~ ~ (K))) +
  xlab('Days post hatch')+ 
  theme(text = element_text(size = 15))
UMK

#relabelling x axis
UM_K<-UMK + scale_x_discrete("Days post hatch", labels = c("UMNov_" = "175", "UMDec" = "200", "UMJan" = "225", "UMJan_" = "240", "UMFeb" = "268", "UMApril" = "325", "UMMay" = "344"))
UM_K

#adding significance values 
cond_um<- UM_K + annotate("text",
                           x = 1, y = 0.43,
                           label = "a",
                           size = 6)+
  annotate("text",
           x = 2, y = 0.4,
           label = "a",
           size = 6)+
  annotate("text",
           x = 3, y = 0.42,
           label = "a",
           size = 6)+
  annotate("text",
           x = 4, y = 0.4,
           label = "a",
           size = 6)+
  annotate("text",
           x = 5, y = 0.4,
           label = "a",
           size = 6)+
  annotate("text",
           x = 6, y = 0.43,
           label = "a",
           size = 6)
cond_um
#save graph 

ggsave(filename = "conditionUM.tiff", plot = cond_um, dpi = 600, height = 10, width = 12)

################################################
##standard growth rate 
##filter out november 
UMSGRonly<- Master[c(646:969), c(1:13)]

UMSGR<-UMSGRonly %>% 
  filter(!is.na(CC_Pop))%>% 
  mutate(CC_Pop = fct_relevel(CC_Pop,"UMDec",
                              "UMJan","UMJan_","UMFeb","UMApril","UMMay"))%>%
  ggplot(UMSGRonly, mapping = aes(x = CC_Pop, y= SGR)) +
  geom_violin(lwd = 0.5)+ 
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.5,
               color = "grey")+
  geom_point(shape = 16) +
  theme_classic() +
  ylab(bquote('Standard growth rate'~ ~ (g/day))) +
  xlab('Days post hatch')+ 
  theme(text = element_text(size = 15))
UMSGR

#relabelling x axis 
UM_SGR<-UMSGR + scale_x_discrete("Days post hatch", labels = c("UMNov_" = "175", "UMDec" = "200", "UMJan" = "225", "UMJan_" = "240", "UMFeb" = "268", "UMApril" = "325", "UMMay" = "344"))
UM_SGR

#adding significance values 
SGR_um<- UM_SGR + annotate("text",
                           x = 1, y = 9.5,
                           label = "ab",
                           size = 6)+
  annotate("text",
           x = 2, y = 12,
           label = "ab",
           size = 6)+
  annotate("text",
           x = 3, y = 11,
           label = "ab",
           size = 6)+
  annotate("text",
           x = 4, y = 12,
           label = "a",
           size = 6)+
  annotate("text",
           x = 5, y = 4,
           label = "b",
           size = 6)+
  annotate("text",
           x = 6, y = 12,
           label = "b",
           size = 6)
SGR_um
#save graph 

ggsave(filename = "SGRUM.tiff", plot = SGR_um, dpi = 600, height = 10, width = 12)


##########################
##combine the 4 graphs into one 

UMBiomet<-(leng_um | mass_um) /
  (cond_um| SGR_um) +
  plot_annotation(tag_levels = 'A')
UMBiomet

##the dots are a bit large.. but hopefully saving will help 
ggsave(filename = "Comb_UMbiomet.tiff", plot = UMBiomet, dpi = 600, height = 12, width = 14)
