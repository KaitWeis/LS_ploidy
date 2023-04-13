###############################################
#
#MNV COMPARISON BTWN UM POPULATION OVER TIME  
#
# Created by KW October 2022 
# modified by K.Weisgerber on 4/13/23
#
##############################################
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

UMonly<- Master[c(599:976), c(1:2)]

#check for outliers 
UMout<-UMonly %>% 
  group_by(CC_Pop) %>% 
  identify_outliers(avg_nuclei)

#normality for each group of this measure 
UMonly %>% 
  group_by(CC_Pop) %>% 
  shapiro_test(avg_nuclei)
#groups are not normal 

ggqqplot(UMonly, "avg_nuclei", facet.by = "CC_Pop")
#plots look good 

#levens for homogentiy of variance 
leveneTest(avg_nuclei ~ CC_Pop, data = UMonly)
#       Df F value Pr(>F)
#group   6  0.5781 0.7478

UMmod<- lm(avg_nuclei ~ CC_Pop, data = UMonly)
UMmod

#Coefficients:
#(Intercept)   CC_PopUMDec   CC_PopUMFeb   CC_PopUMJan  CC_PopUMJan_   CC_PopUMMay  
#54.6745        1.7275        0.3611        1.1895        2.2705       -1.2395  
#CC_PopUMNov_  
#0.5335 

summary(UMmod)

#Residuals:
#Min      1Q  Median      3Q     Max 
#-4.9509 -1.9770 -0.7658  0.6781 29.2960 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   54.6745     0.6498  84.134   <2e-16 ***
#  CC_PopUMDec    1.7275     0.8834   1.956   0.0513 .  
#CC_PopUMFeb    0.3611     0.8925   0.405   0.6861    
#CC_PopUMJan    1.1895     0.8879   1.340   0.1812    
#CC_PopUMJan_   2.2705     0.8879   2.557   0.0110 *  
#CC_PopUMMay   -1.2395     0.7695  -1.611   0.1081    
#CC_PopUMNov_   0.5335     0.8790   0.607   0.5443    
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 4.058 on 356 degrees of freedom
#(15 observations deleted due to missingness)
#Multiple R-squared:  0.08349,	Adjusted R-squared:  0.06804 
#F-statistic: 5.405 on 6 and 356 DF,  p-value: 2.313e-05

##suggesting there is a difference btwn these times 

ggqqplot(residuals(UMmod))
#bit of a tail at the end 

shapiro_test(residuals(UMmod))

#test homogentiy of variance of the model as a whole 
plot(UMmod, 1)
#  variable         statistic  p.value
#<chr>                <dbl>    <dbl>
#  1 residuals(UMmod)     0.619 1.48e-27

check_model(UMmod, check = "all")
#looks pretty good 

anova(UMmod)
#Response: avg_nuclei
#Df Sum Sq Mean Sq F value    Pr(>F)    
#CC_Pop      6  534.1  89.016  5.4048 2.313e-05 ***
#  Residuals 356 5863.2  16.470  

AIC(UMmod)
#2056.034

emmeans(UMmod, pairwise ~ CC_Pop, adjust = "tukey")
#$emmeans
#CC_Pop  emmean    SE  df lower.CL upper.CL
#UMApril   54.7 0.650 356     53.4     56.0
#UMDec     56.4 0.598 356     55.2     57.6
#UMFeb     55.0 0.612 356     53.8     56.2
#UMJan     55.9 0.605 356     54.7     57.1
#UMJan_    56.9 0.605 356     55.8     58.1
#UMMay     53.4 0.412 356     52.6     54.2
#UMNov_    55.2 0.592 356     54.0     56.4
#Confidence level used: 0.95 
#$contrasts
#contrast         estimate    SE  df t.ratio p.value
#UMApril - UMDec    -1.727 0.883 356  -1.956  0.4451
#UMApril - UMFeb    -0.361 0.893 356  -0.405  0.9997
#UMApril - UMJan    -1.190 0.888 356  -1.340  0.8327
#UMApril - UMJan_   -2.270 0.888 356  -2.557  0.1426
#UMApril - UMMay     1.239 0.769 356   1.611  0.6755
#UMApril - UMNov_   -0.534 0.879 356  -0.607  0.9966
#UMDec - UMFeb       1.366 0.856 356   1.597  0.6846
#UMDec - UMJan       0.538 0.851 356   0.632  0.9957
#UMDec - UMJan_     -0.543 0.851 356  -0.638  0.9955
#UMDec - UMMay       2.967 0.727 356   4.084  0.0011
#UMDec - UMNov_      1.194 0.842 356   1.418  0.7916
#UMFeb - UMJan      -0.828 0.860 356  -0.963  0.9615
#UMFeb - UMJan_     -1.909 0.860 356  -2.219  0.2878
#UMFeb - UMMay       1.601 0.738 356   2.170  0.3147
#UMFeb - UMNov_     -0.172 0.851 356  -0.203  1.0000
#UMJan - UMJan_     -1.081 0.856 356  -1.263  0.8681
#UMJan - UMMay       2.429 0.732 356   3.318  0.0172
#UMJan - UMNov_      0.656 0.846 356   0.775  0.9872
#UMJan_ - UMMay      3.510 0.732 356   4.795  <.0001
#UMJan_ - UMNov_     1.737 0.846 356   2.052  0.3838
#UMMay - UMNov_     -1.773 0.721 356  -2.458  0.1780


sjstats::anova_stats(car::Anova(UMmod, type = 3)) %>% dplyr::select(1:7)
#term      |    sumsq | meansq |  df | statistic | p.value | etasq
#-----------------------------------------------------------------
#CC_Pop    |  534.096 | 89.016 |   6 |     5.405 |  < .001 | 0.083
#Residuals | 5863.217 | 16.470 | 356 |           |         |  

model_performance(UMmod)
#AIC      |     AICc |      BIC |    R2 | R2 (adj.) |  RMSE | Sigma
#------------------------------------------------------------------
#2056.034 | 2056.441 | 2087.189 | 0.083 |     0.068 | 4.019 | 4.058