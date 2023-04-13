######################################
#
#REPEATED MEASURES OF NUCLEI VOLUME OVER 14 DAYS 
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
library(afex)
library(PMCMRplus)
library(rstatix)

setwd("/Users/kaitl/Desktop/LS_ploidy/raw_data/")

#read in the repeated sample data 
rep<-read.csv("rep_CC_data.csv")

#convert into long format and change id and time into factors 
rep1<-rep%>%
  gather(key = "time", value = "score", t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,
         t11,t12,t13,t14)%>%
  convert_as_factor(id, time)

#check the new data format 
head(rep1, 4)

#compute the summary stats by time 
rep1%>%
  group_by(time)%>%
  get_summary_stats(score, type = "mean_sd")

#check the normality assumptions 

##starting with outliers 

out<-rep1%>%
  group_by(time)%>%
  identify_outliers(score)
print(out)

##normal distribution 

rep1%>%
  group_by(time)%>%
  shapiro_test(score)
#large p value is normal distribution 

##try qqplot 

ggqqplot(rep1, "score", facet.by = "time")

##the qq plots look mostly normal, but all of the p values are very small 

#use the non parametric equivalent of the repeated ANOVA - Friedmans test 

##equation for this is a~b|c ; a is dependent, b is within subject and c is the identifier 

friedman_test(score ~ time | id, data = rep1)

## A tibble: 1 × 6
#.y.       n statistic    df        p method       
#* <chr> <int>     <dbl> <dbl>    <dbl> <chr>        
#  1 score    10      34.9    13 0.000890 Friedman test

#this is suggesting a significant difference 

#look at pairwise comparisons with Nemenyi-Wilcoxon-Wilcox all-pairs test for a two-way balanced complete block design 

rep_comp = frdAllPairsNemenyiTest(score ~ time | id, data = rep1)
rep_comp

#   t1    t10   t11   t12   t13   t14   t2    t3    t4    t5    t6    t7    t8   
#t10 0.328 -     -     -     -     -     -     -     -     -     -     -     -    
#t11 0.931 1.000 -     -     -     -     -     -     -     -     -     -     -    
#t12 0.328 1.000 1.000 -     -     -     -     -     -     -     -     -     -    
#t13 0.782 1.000 1.000 1.000 -     -     -     -     -     -     -     -     -    
#t14 0.843 1.000 1.000 1.000 1.000 -     -     -     -     -     -     -     -    
#t2  0.931 1.000 1.000 1.000 1.000 1.000 -     -     -     -     -     -     -    
#t3  1.000 0.892 1.000 0.892 0.998 0.999 1.000 -     -     -     -     -     -    
#t4  0.913 1.000 1.000 1.000 1.000 1.000 1.000 1.000 -     -     -     -     -    
#t5  0.437 1.000 1.000 1.000 1.000 1.000 1.000 0.946 1.000 -     -     -     -    
#t6  0.536 1.000 1.000 1.000 1.000 1.000 1.000 0.973 1.000 1.000 -     -     -    
#t7  0.536 1.000 1.000 1.000 1.000 1.000 1.000 0.973 1.000 1.000 1.000 -     -    
#t8  0.091 1.000 0.969 1.000 0.996 0.992 0.969 0.556 0.977 1.000 1.000 1.000 -    
#t9  0.437 1.000 1.000 1.000 1.000 1.000 1.000 0.946 1.000 1.000 1.000 1.000 1.000


###accepting significance level 0.05, non of these comparisons are significant. 

#check the effectsize with Kendall's W 

friedman_effsize(score ~ time | id, data = rep1)
## A tibble: 1 × 5
#.y.       n effsize method    magnitude
#* <chr> <int>   <dbl> <chr>     <ord>    
#  1 score    10   0.268 Kendall W small 

#small effect size 


#just to be sure testing pairwise with wilcox test and bonferroni correction 

pwc<- rep1 %>% 
  wilcox_test(score ~ time, paired = TRUE, p.adjust.method = "bonferroni")
pwc

#none of these are significant either 

###the MNV did not change over the 14 days 
