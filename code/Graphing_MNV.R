###############################################
#
#GRAPHING MNV COMPARISON IN BOTH ANAYLSES   
#
# Created by KW October 2022 
# modified by K.Weisgerber on 4/15/23
#
##############################################

library(ggbeeswarm)
library(patchwork)

#starting with the main 5 populations 
##making a violoin plot with mean as a crossbar 
CCfive<-five %>% 
  filter(!is.na(CC_Pop))%>% 
  mutate(CC_Pop = fct_relevel(CC_Pop,"GRH_BR", "GRH_LR",
                              "NR","PDB", "UMMay"))%>%
  ggplot(five, mapping = aes(x = CC_Pop, y= avg_nuclei)) +
  geom_violin(lwd = 0.5)+ 
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.5,
               color = "grey")+
  geom_point(shape = 16) + 
  theme_classic() +
  ylab(bquote('Average modal nuclei volume'~ ~ (fL))) +
  xlab('Sample population')+ 
  theme(text = element_text(size = 15))
CCfive

#making the axis labels full population names 
cc_five<- CCfive + scale_x_discrete("Sample population", labels = c("GRH_BR" = "Birthday Rapids", "GRH_LR" = "Landing River", "NR" = "Nelson River", "PDB" = "Pointe du Bois", "UMMay" = "University" ))
cc_five

#adding significance labels 
gg_five<- cc_five +   annotate("text",
                               x = 1, y = 70,
                               label = "a",
                               size = 6)+
  annotate("text",
           x = 2, y = 78,
           label = "a",
           size = 6)

gg_five
#final graph for 5 populations 

#save file 
ggsave(filename = "MNV_five.tiff", plot = gg_five, dpi = 600, height = 10, width = 14)


################################################
#graphin UM only over time 

UMCC<-UMonly %>% 
  filter(!is.na(CC_Pop))%>% 
  mutate(CC_Pop = fct_relevel(CC_Pop,"UMNov_", "UMDec",
                              "UMJan","UMJan_","UMFeb","UMApril","UMMay"))%>%
  ggplot(UMonly, mapping = aes(x = CC_Pop, y= avg_nuclei)) +
  geom_violin(lwd = 0.5)+ 
  stat_summary(fun = "mean",
               geom = "crossbar",
               width = 0.5,
               color = "grey")+
  geom_point(shape = 16) +
  theme_classic() +
  ylab(bquote('Average modal nuclei volume'~ ~ (fL))) +
  xlab('Days post hatch')+ 
  theme(text = element_text(size = 15))
UMCC

#relabelling x axis 
UM_CC<-UMCC + scale_x_discrete("Days post hatch", labels = c("UMNov_" = "175", "UMDec" = "200", "UMJan" = "225", "UMJan_" = "240", "UMFeb" = "268", "UMApril" = "325", "UMMay" = "344"))
UM_CC

#adding significance values 
gg_um<- UM_CC + annotate("text",
                           x = 1, y = 67,
                           label = "ab",
                           size = 6)+
  annotate("text",
           x = 2, y = 84,
           label = "b",
           size = 6)+
  annotate("text",
           x = 3, y = 85,
           label = "b",
           size = 6)+
  annotate("text",
           x = 4, y = 86,
           label = "b",
           size = 6)+
  annotate("text",
           x = 5, y = 87,
           label = "ab",
           size = 6)+
  annotate("text",
           x = 6, y = 80,
           label = "ab",
           size = 6)+
  annotate("text",
           x = 7, y = 81,
           label = "a",
           size = 6)
gg_um
#save graph 

ggsave(filename = "MNV_UM.tiff", plot = gg_um, dpi = 600, height = 10, width = 14)

##################################
##combine the two plots into a single figure 

CCpatch<-(gg_five) /
  (gg_um) +
  plot_annotation(tag_levels = 'A')
CCpatch
ggsave(filename = "MNV_comb_final.tiff", plot = CCpatch, dpi = 600, height = 12, width = 14)
