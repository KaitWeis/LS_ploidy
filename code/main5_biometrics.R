######################################
#
#BIOMETRICS ANALYSIS BETWEEN FIVE MAIN POPULATIONS  
#
# Created by KW October 2022 
# modified by K.Weisgerber on 4/15/23
#
############################################
#load libraries
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(emmeans)
library(rstatix)
library(performance)

setwd("/Users/kaitl/Desktop/LS_ploidy/raw_data/")

#this is the master data file, will probably have to subset 
Master<-read.csv("comb_ploidy_data.csv")

#this makes it just the data for the five main populations 
main<- Master[c(1:598,865:1020), ]

###running normality, linear models, and anovas on the models 

###############################

### length 

#check for outliers 
TLout<-main %>% 
  group_by(CC_Pop) %>% 
  identify_outliers(TLeng_mm)

#normality for each group of this measure 
main %>% 
  group_by(CC_Pop) %>% 
  shapiro_test(TLeng_mm)

ggqqplot(main, "TLeng_mm", facet.by = "CC_Pop")
#slight tail in NR 

TLmod<- lm(TLeng_mm ~ CC_Pop, data = main)
TLmod
#Coefficients:
#(Intercept)  CC_PopGRH_LR      CC_PopNR     CC_PopPDB   CC_PopUMMay  
#208.45          9.02       1098.88        833.22        -36.87 

summary(TLmod)
#Residuals:
#Min      1Q  Median      3Q     Max 
#-479.67  -43.33   -1.45   26.42 3092.67 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    208.45      15.86  13.145   <2e-16 ***
#  CC_PopGRH_LR     9.02      22.43   0.402    0.688    
#CC_PopNR      1098.88      33.60  32.709   <2e-16 ***
#  CC_PopPDB      833.22      19.44  42.853   <2e-16 ***
#  CC_PopUMMay    -36.87      24.71  -1.492    0.136    
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 194.2 on 741 degrees of freedom
#(8 observations deleted due to missingness)
#Multiple R-squared:  0.8373,	Adjusted R-squared:  0.8364 
#F-statistic: 953.5 on 4 and 741 DF,  p-value: < 2.2e-16

ggqqplot(residuals(TLmod))

shapiro_test(residuals(TLmod))
#residuals are non normal overall; because grouping widely different ages? 

#test homogentiy of variance of the model as a whole 
plot(TLmod, 1)

check_model(TLmod, check = "all")

anova(TLmod)
#Response: TLeng_mm
#Df    Sum Sq  Mean Sq F value    Pr(>F)    
#CC_Pop      4 143862062 35965516  953.47 < 2.2e-16 ***
#Residuals 741  27951064    37721      

AIC(TLmod)
#9985.362

emmeans(TLmod, pairwise ~ CC_Pop, adjust = "tukey")
#$emmeans
#CC_Pop emmean   SE  df lower.CL upper.CL
#GRH_BR    208 15.9 741      177      240
#GRH_LR    217 15.9 741      186      249
#NR       1307 29.6 741     1249     1365
#PDB      1042 11.3 741     1020     1064
#UMMay     172 19.0 741      134      209
#Confidence level used: 0.95 
#$contrasts
#contrast        estimate   SE  df t.ratio p.value
#GRH_BR - GRH_LR    -9.02 22.4 741  -0.402  0.9945
#GRH_BR - NR     -1098.88 33.6 741 -32.709  <.0001
#GRH_BR - PDB     -833.22 19.4 741 -42.853  <.0001
#GRH_BR - UMMay     36.87 24.7 741   1.492  0.5682
#GRH_LR - NR     -1089.86 33.6 741 -32.440  <.0001
#GRH_LR - PDB     -824.20 19.4 741 -42.389  <.0001
#GRH_LR - UMMay     45.89 24.7 741   1.857  0.3419
#NR - PDB          265.66 31.7 741   8.385  <.0001
#NR - UMMay       1135.75 35.2 741  32.299  <.0001
#PDB - UMMay       870.09 22.0 741  39.475  <.0001

sjstats::anova_stats(car::Anova(TLmod, type = 3)) %>% dplyr::select(1:7)
#term      |     sumsq |    meansq |  df | statistic | p.value | etasq
#---------------------------------------------------------------------
#CC_Pop    | 1.439e+08 | 3.597e+07 |   4 |   953.468 |  < .001 | 0.837
#Residuals | 2.795e+07 | 37720.735 | 741 |           |         |  

##need to run a non parametric as well 

kruskal.test(TLeng_mm ~ CC_Pop, data = main)

#data:  TLeng_mm by CC_Pop
#Kruskal-Wallis chi-squared = 591.02, df = 4, p-value < 2.2e-16

kruskal_effsize(TLeng_mm ~ CC_Pop, data = main)
#  .y.          n effsize method  magnitude
#* <chr>    <int>   <dbl> <chr>   <ord>    
#  1 TLeng_mm   754   0.784 eta2[H] large 

pairwise.wilcox.test(main$TLeng_mm, main$CC_Pop, p.adjust.method = "BH") 

#       GRH_BR  GRH_LR  NR      PDB    
#GRH_LR 8.9e-07 -       -       -      
#NR     < 2e-16 < 2e-16 -       -      
#PDB    < 2e-16 < 2e-16 4e-04   -      
#UMMay  5.9e-15 < 2e-16 < 2e-16 < 2e-16

##be interested to compare just the hatchery to each other, and just the wild 

###############################

### wet mass  

#check for outliers 
massout<-main %>% 
  group_by(CC_Pop) %>% 
  identify_outliers(mass_g)

#normality for each group of this measure 
main %>% 
  group_by(CC_Pop) %>% 
  shapiro_test(mass_g)

ggqqplot(main, "mass_g", facet.by = "CC_Pop")
#NR and PDB have tails both ends 

Massmod<- lm(mass_g ~ CC_Pop, data = main)
Massmod
#Coefficients:
#(Intercept)  CC_PopGRH_LR      CC_PopNR     CC_PopPDB   CC_PopUMMay  
#34.3494        0.9309     7576.8225     7358.5197      -12.3840  

summary(Massmod)
#Residuals:
#Min      1Q  Median      3Q     Max 
#-6941.2  -892.9    -3.5    10.7 17555.1 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    34.3494   240.6407   0.143    0.887    
#CC_PopGRH_LR    0.9309   340.3174   0.003    0.998    
#CC_PopNR     7576.8225   509.8161  14.862   <2e-16 ***
#  CC_PopPDB    7358.5197   295.0530  24.940   <2e-16 ***
#  CC_PopUMMay   -12.3840   368.0533  -0.034    0.973    
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 2947 on 748 degrees of freedom
#(1 observation deleted due to missingness)
#Multiple R-squared:  0.6106,	Adjusted R-squared:  0.6085 
#F-statistic: 293.2 on 4 and 748 DF,  p-value: < 2.2e-16

ggqqplot(residuals(Massmod))

shapiro_test(residuals(Massmod))
##not normal 

#test homogentiy of variance of the model as a whole 
plot(Massmod, 1)

check_model(Massmod, check = "all")

anova(Massmod)
#           Df     Sum Sq    Mean Sq F value    Pr(>F)    
#CC_Pop      4 1.0189e+10 2547145337  293.24 < 2.2e-16 ***
#Residuals 748 6.4973e+09    8686193   

AIC(Massmod)
#14174.77

emmeans(Massmod, pairwise ~ CC_Pop, adjust = "tukey")
#similar to the length, it is similar between hatchery groups and different otherwise; also similar btwn the two wild 

#$emmeans
#CC_Pop emmean  SE  df lower.CL upper.CL
#GRH_BR   34.3 241 748     -438      507
#GRH_LR   35.3 241 748     -437      508
#NR     7611.2 449 748     6729     8494
#PDB    7392.9 171 748     7058     7728
#UMMay    22.0 278 748     -525      569
#Confidence level used: 0.95 
#$contrasts
#contrast         estimate  SE  df t.ratio p.value
#GRH_BR - GRH_LR    -0.931 340 748  -0.003  1.0000
#GRH_BR - NR     -7576.822 510 748 -14.862  <.0001
#GRH_BR - PDB    -7358.520 295 748 -24.940  <.0001
#GRH_BR - UMMay     12.384 368 748   0.034  1.0000
#GRH_LR - NR     -7575.892 510 748 -14.860  <.0001
#GRH_LR - PDB    -7357.589 295 748 -24.937  <.0001
#GRH_LR - UMMay     13.315 368 748   0.036  1.0000
#NR - PDB          218.303 481 748   0.454  0.9912
#NR - UMMay       7589.207 529 748  14.354  <.0001
#PDB - UMMay      7370.904 327 748  22.565  <.0001

sjstats::anova_stats(car::Anova(Massmod, type = 3)) %>% dplyr::select(1:7)

#term      |     sumsq |    meansq |  df | statistic | p.value | etasq
#---------------------------------------------------------------------
#  CC_Pop    | 1.019e+10 | 2.547e+09 |   4 |   293.241 |  < .001 | 0.611
#Residuals | 6.497e+09 | 8.686e+06 | 748 |           |         |  

##try a non parametric 
kruskal.test(mass_g ~ CC_Pop, data = main)
#data:  mass_g by CC_Pop
#Kruskal-Wallis chi-squared = 587.36, df = 4, p-value < 2.2e-16

kruskal_effsize(mass_g ~ CC_Pop, data = main)
#  .y.        n effsize method  magnitude
#* <chr>  <int>   <dbl> <chr>   <ord>    
#  1 mass_g   754   0.779 eta2[H] large 

pairwise.wilcox.test(main$mass_g, main$CC_Pop, p.adjust.method = "BH")

#       GRH_BR GRH_LR NR     PDB   
#GRH_LR 0.081  -      -      -     
#NR     <2e-16 <2e-16 -      -     
#PDB    <2e-16 <2e-16 0.531  -     
#UMMay  <2e-16 <2e-16 <2e-16 <2e-16

###############################

### condition factor  
##needs to be converted to a number 
data_K<- transform(main, K = as.numeric(K))

#check for outliers 
Kout<-data_K %>% 
  group_by(CC_Pop) %>% 
  identify_outliers(K)

#normality for each group of this measure 
data_K %>% 
  group_by(CC_Pop) %>% 
  shapiro_test(K)

ggqqplot(data_K, "K", facet.by = "CC_Pop")
##NR bad tail on negative side 

Kmod<- lm(K ~ CC_Pop, data = data_K)
Kmod
#Coefficients:
#(Intercept)  CC_PopGRH_LR      CC_PopNR     CC_PopPDB   CC_PopUMMay  
#0.37289      -0.03424       0.08465       0.22830      -0.03290  

summary(Kmod)
#Residuals:
#Min       1Q   Median       3Q      Max 
#-0.53527 -0.03374 -0.00560  0.02865  0.43212 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.372889   0.007436  50.146  < 2e-16 ***
#  CC_PopGRH_LR -0.034235   0.010516  -3.256  0.00118 ** 
#  CC_PopNR      0.084652   0.015754   5.373 1.04e-07 ***
#  CC_PopPDB     0.228304   0.009117  25.040  < 2e-16 ***
#  CC_PopUMMay  -0.032899   0.011588  -2.839  0.00465 ** 
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 0.09107 on 741 degrees of freedom
#(8 observations deleted due to missingness)
#Multiple R-squared:  0.6372,	Adjusted R-squared:  0.6353 
#F-statistic: 325.4 on 4 and 741 DF,  p-value: < 2.2e-16

ggqqplot(residuals(Kmod))
#large tails 

shapiro_test(residuals(Kmod))
#not normal 

#test homogentiy of variance of the model as a whole 
plot(Kmod, 1)

check_model(Kmod, check = "all")

anova(Kmod)
#Response: K
#Df Sum Sq Mean Sq F value    Pr(>F)    
#CC_Pop      4 10.795 2.69885  325.39 < 2.2e-16 ***
#Residuals 741  6.146 0.00829       

AIC(Kmod)
#-1450.939

emmeans(Kmod, pairwise ~ CC_Pop, adjust = "tukey")
#$emmeans
#CC_Pop emmean      SE  df lower.CL upper.CL
#GRH_BR  0.373 0.00744 741    0.358    0.387
#GRH_LR  0.339 0.00744 741    0.324    0.353
#NR      0.458 0.01389 741    0.430    0.485
#PDB     0.601 0.00528 741    0.591    0.612
#UMMay   0.340 0.00889 741    0.323    0.357
#Confidence level used: 0.95 
#$contrasts
#contrast        estimate      SE  df t.ratio p.value
#GRH_BR - GRH_LR  0.03424 0.01052 741   3.256  0.0104
#GRH_BR - NR     -0.08465 0.01575 741  -5.373  <.0001
#GRH_BR - PDB    -0.22830 0.00912 741 -25.040  <.0001
#GRH_BR - UMMay   0.03290 0.01159 741   2.839  0.0374
#GRH_LR - NR     -0.11889 0.01575 741  -7.547  <.0001
#GRH_LR - PDB    -0.26254 0.00912 741 -28.795  <.0001
#GRH_LR - UMMay  -0.00134 0.01159 741  -0.115  1.0000
#NR - PDB        -0.14365 0.01486 741  -9.669  <.0001
#NR - UMMay       0.11755 0.01649 741   7.129  <.0001
#PDB - UMMay      0.26120 0.01034 741  25.272  <.0001

sjstats::anova_stats(car::Anova(Kmod, type = 3)) %>% dplyr::select(1:7)
#term      |  sumsq | meansq |  df | statistic | p.value | etasq
#---------------------------------------------------------------
#CC_Pop    | 10.795 |  2.699 |   4 |   325.390 |  < .001 | 0.637
#Residuals |  6.146 |  0.008 | 741 |           |         |    


###non parametric
kruskal.test(K ~ CC_Pop, data = data_K)
#data:  K by CC_Pop
#Kruskal-Wallis chi-squared = 515.72, df = 4, p-value < 2.2e-16

kruskal_effsize(K ~ CC_Pop, data = data_K)
#  .y.       n effsize method  magnitude
#* <chr> <int>   <dbl> <chr>   <ord>    
#  1 K       754   0.683 eta2[H] large 

pairwise.wilcox.test(data_K$K, data_K$CC_Pop, p.adjust.method = "BH")
#       GRH_BR  GRH_LR  NR      PDB    
#GRH_LR < 2e-16 -       -       -      
#NR     2.6e-08 1.1e-08 -       -      
#PDB    < 2e-16 < 2e-16 3.7e-07 -      
#UMMay  8.8e-16 0.016   3.0e-07 < 2e-16

#all significantly different K 


#######################################################
##analyze MNV in 5 populations accounting for mass/ length 

#already know this data is not normal, can run interaction in kruskal? 

#length first 

FLCC<-lm(avg_nuclei ~ TLeng_mm*CC_Pop, data = main)
summary(FLCC)

#Residuals:
#Min       1Q   Median       3Q      Max 
#-29.1704  -0.9024  -0.1073   0.5287  24.4149 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           53.106920   2.209246  24.038   <2e-16 ***
#  TLeng_mm              -0.006676   0.010558  -0.632   0.5274    
#CC_PopGRH_LR          -4.788506   3.334105  -1.436   0.1514    
#CC_PopNR              -2.528631   2.349466  -1.076   0.2822    
#CC_PopPDB             -4.782361   2.361973  -2.025   0.0433 *  
#  CC_PopUMMay            4.236748   2.550625   1.661   0.0971 .  
#TLeng_mm:CC_PopGRH_LR  0.021684   0.015574   1.392   0.1642    
#TLeng_mm:CC_PopNR      0.006261   0.010572   0.592   0.5539    
#TLeng_mm:CC_PopPDB     0.006318   0.010588   0.597   0.5509    
#TLeng_mm:CC_PopUMMay  -0.013947   0.012657  -1.102   0.2709    
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 2.365 on 721 degrees of freedom
#(23 observations deleted due to missingness)
#Multiple R-squared:  0.4478,	Adjusted R-squared:  0.4409 
#F-statistic: 64.96 on 9 and 721 DF,  p-value: < 2.2e-16

anova(FLCC)
#                 Df Sum Sq Mean Sq  F value  Pr(>F)    
#TLeng_mm          1 2237.7 2237.66 400.1825 < 2e-16 ***
#  CC_Pop            4  972.4  243.11  43.4772 < 2e-16 ***
#  TLeng_mm:CC_Pop   4   58.9   14.73   2.6348 0.03311 *  
#  Residuals       721 4031.5    5.59   

#there is a significant interaction between length and population 
##as well as the significant difference within population and within lengths that were already analyzed 

AIC(FLCC)
#3344.663
#this is high.. but LOWER than the individual models without interaction 

FMCC<-lm(avg_nuclei ~ mass_g*CC_Pop, data = main)
summary(FMCC)
#Residuals:
#Min       1Q   Median       3Q      Max 
#-29.2144  -0.9069  -0.1421   0.4893  24.2643 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)         51.551629   0.747767  68.941  < 2e-16 ***
#  mass_g               0.004763   0.021032   0.226 0.820910    
#CC_PopGRH_LR        -0.985508   1.174983  -0.839 0.401889    
#CC_PopNR            -1.192229   1.034063  -1.153 0.249307    
#CC_PopPDB           -3.509941   0.794942  -4.415 1.16e-05 ***
#  CC_PopUMMay          3.287703   0.864206   3.804 0.000154 ***
#  mass_g:CC_PopGRH_LR  0.024035   0.032748   0.734 0.463223    
#mass_g:CC_PopNR     -0.004805   0.021033  -0.228 0.819339    
#mass_g:CC_PopPDB    -0.004775   0.021032  -0.227 0.820459    
#mass_g:CC_PopUMMay  -0.062629   0.025755  -2.432 0.015267 *  
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 2.363 on 728 degrees of freedom
#(16 observations deleted due to missingness)
#Multiple R-squared:  0.4437,	Adjusted R-squared:  0.4368 
#F-statistic: 64.51 on 9 and 728 DF,  p-value: < 2.2e-16

anova(FMCC)
#Response: avg_nuclei
#Df Sum Sq Mean Sq  F value    Pr(>F)    
#mass_g          1 1725.0 1724.99 308.9530 < 2.2e-16 ***
#  CC_Pop          4 1424.0  355.99  63.7599 < 2.2e-16 ***
#  mass_g:CC_Pop   4   92.9   23.22   4.1592  0.002447 ** 
#  Residuals     728 4064.7    5.58  

##again significant interaction btwn and within each 

AIC(FMCC)
#3375.487


FKCC<-lm(avg_nuclei ~ K*CC_Pop, data = data_K)
summary(FKCC)
#Residuals:
#Min       1Q   Median       3Q      Max 
#-29.2736  -0.8938  -0.1365   0.4822  23.8052 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)      45.995      1.816  25.325  < 2e-16 ***
#  K                15.340      4.843   3.167 0.001603 ** 
#  CC_PopGRH_LR      7.745      3.413   2.269 0.023555 *  
#  CC_PopNR          3.836      2.041   1.880 0.060570 .  
#CC_PopPDB         1.796      1.961   0.916 0.360065    
#CC_PopUMMay       9.915      2.227   4.451 9.88e-06 ***
#  K:CC_PopGRH_LR  -21.711      9.795  -2.217 0.026964 *  
#  K:CC_PopNR      -14.894      5.194  -2.868 0.004258 ** 
#  K:CC_PopPDB     -15.074      4.992  -3.020 0.002621 ** 
#  K:CC_PopUMMay   -21.897      6.077  -3.603 0.000336 ***
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 2.361 on 721 degrees of freedom
#(23 observations deleted due to missingness)
#Multiple R-squared:  0.4495,	Adjusted R-squared:  0.4426 
#F-statistic: 65.41 on 9 and 721 DF,  p-value: < 2.2e-16

anova(FKCC)
#Response: avg_nuclei
#Df Sum Sq Mean Sq  F value    Pr(>F)    
#K           1 1807.9 1807.92 324.3241 < 2.2e-16 ***
#  CC_Pop      4 1396.8  349.21  62.6442 < 2.2e-16 ***
#  K:CC_Pop    4   76.7   19.16   3.4378  0.008535 ** 
#  Residuals 721 4019.2    5.57

AIC(FKCC)
#3342.415

