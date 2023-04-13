######################################
#
#MNV COMPARISON BTWN 5 MAIN POPULATIONS  
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

#this makes it just the coulter counter data for the five main populations 
five<- Master[c(1:598,865:1020), c(1:2)]

#lets check the normality of this data 

#check for outliers 
fiveout<-five %>% 
  group_by(CC_Pop) %>% 
  identify_outliers(avg_nuclei)
###there are some extreme outliers 

#normality for each group of this measure 
five %>% 
  group_by(CC_Pop) %>% 
  shapiro_test(avg_nuclei)

#### none are normal 

ggqqplot(five, "avg_nuclei", facet.by = "CC_Pop")

#levens for homogentiy of variance 
leveneTest(avg_nuclei ~ CC_Pop, data = five)
#group   4  5.7701 0.0001417 ***

#evaulate as a linear model 

fivemod<- lm(avg_nuclei ~ CC_Pop, data = five)
fivemod
#Coefficients:
#(Intercept)  CC_PopGRH_LR      CC_PopNR     CC_PopPDB   CC_PopUMMay  
#51.7152       -0.1331       -1.6892       -3.7639        1.7198  
summary(fivemod)

#Residuals:
#Min       1Q   Median       3Q      Max 
#-29.3575  -0.8727  -0.1533   0.4556  24.1988 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   51.7152     0.1944 266.069  < 2e-16 ***
#  CC_PopGRH_LR  -0.1331     0.2749  -0.484    0.628    
#CC_PopNR      -1.6892     0.4081  -4.139 3.89e-05 ***
#  CC_PopPDB     -3.7639     0.2383 -15.794  < 2e-16 ***
#  CC_PopUMMay    1.7198     0.3102   5.545 4.11e-08 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 2.381 on 734 degrees of freedom
#(15 observations deleted due to missingness)
#Multiple R-squared:  0.4308,	Adjusted R-squared:  0.4277 
#F-statistic: 138.9 on 4 and 734 DF,  p-value: < 2.2e-16


ggqqplot(residuals(fivemod))
#the model has tails at the start and end 

shapiro_test(residuals(fivemod))
## A tibble: 1 × 3
#variable           statistic  p.value
#<chr>                  <dbl>    <dbl>
#  1 residuals(fivemod)     0.569 3.81e-39

#test homogentiy of variance of the model as a whole 
plot(fivemod, 1)

check_model(fivemod, check = "all")

AIC(fivemod)
#3386.067

kruskal.test(avg_nuclei ~ CC_Pop, data = five)

#Kruskal-Wallis chi-squared = 513.21, df = 4, p-value < 2.2e-16


pairwise.wilcox.test(five$avg_nuclei, five$CC_Pop, p.adjust.method = "BH")
#       GRH_BR  GRH_LR  NR      PDB    
#  GRH_LR 0.33    -       -       -      
#  NR     1.2e-13 3.4e-14 -       -      
#  PDB    < 2e-16 < 2e-16 < 2e-16 -      
#  UMMay  2.0e-08 1.6e-07 2.8e-15 < 2e-16

##all groups have different MNV except LR and BR based on non-parametric 

######overall the check model isn't bad, the sample size is very large so a parametric test may be okay to use; will also be more useful because of differences in variability 

anova(fivemod)

#Response: avg_nuclei
#Df Sum Sq Mean Sq F value    Pr(>F)    
#CC_Pop      4 3147.5  786.88  138.86 < 2.2e-16 ***
# Residuals 734 4159.5    5.67                      


emmeans(fivemod, pairwise ~ CC_Pop, adjust = "tukey")

#$emmeans
#CC_Pop emmean    SE  df lower.CL upper.CL
#GRH_BR   51.7 0.194 734     51.3     52.1
#GRH_LR   51.6 0.194 734     51.2     52.0
#NR       50.0 0.359 734     49.3     50.7
#PDB      48.0 0.138 734     47.7     48.2
#UMMay    53.4 0.242 734     53.0     53.9
#Confidence level used: 0.95 
#$contrasts
#contrast        estimate    SE  df t.ratio p.value
#GRH_BR - GRH_LR    0.133 0.275 734   0.484  0.9888
#GRH_BR - NR        1.689 0.408 734   4.139  0.0004
#GRH_BR - PDB       3.764 0.238 734  15.794  <.0001
#GRH_BR - UMMay    -1.720 0.310 734  -5.545  <.0001
#GRH_LR - NR        1.556 0.408 734   3.813  0.0014
#GRH_LR - PDB       3.631 0.238 734  15.235  <.0001
#GRH_LR - UMMay    -1.853 0.310 734  -5.974  <.0001
#NR - PDB           2.075 0.384 734   5.396  <.0001
#NR - UMMay        -3.409 0.433 734  -7.879  <.0001
#PDB - UMMay       -5.484 0.278 734 -19.706  <.0001

sjstats::anova_stats(car::Anova(fivemod, type = 3)) %>% dplyr::select(1:7)
#term      |    sumsq |  meansq |  df | statistic | p.value | etasq
#------------------------------------------------------------------
#CC_Pop    | 3147.500 | 786.875 |   4 |   138.856 |  < .001 | 0.431
#Residuals | 4159.461 |   5.667 | 734 |           |         |   

#using the parametric data; results are the same 

model_performance(fivemod)
#AIC      |     AICc |      BIC |    R2 | R2 (adj.) |  RMSE | Sigma
#------------------------------------------------------------------
#3386.067 | 3386.182 | 3413.699 | 0.431 |     0.428 | 2.372 | 2.381