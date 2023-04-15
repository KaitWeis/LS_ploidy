######################################
#
#BIOMETRICS ANALYSIS BETWEEN UM TIME POINTS   
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

setwd("/Users/kaitl/Desktop/LS_ploidy/raw_data/")

#this is the master data file, will probably have to subset 
Master<-read.csv("comb_ploidy_data.csv")

#subset data for UM time points 
UMonly<- Master[c(599:976), c(1:13)]

#check for outliers 
TLout<-UMonly %>% 
  group_by(CC_Pop) %>% 
  identify_outliers(TLeng_mm)

#normality for each group of this measure 
UMonly %>% 
  group_by(CC_Pop) %>% 
  shapiro_test(TLeng_mm)
#all normal except UM may slightly 

ggqqplot(UMonly, "TLeng_mm", facet.by = "CC_Pop")
##looks good 

TLUMmod<- lm(TLeng_mm ~ CC_Pop, data = UMonly)
TLUMmod
#Coefficients:
#(Intercept)   CC_PopUMDec   CC_PopUMFeb   CC_PopUMJan  CC_PopUMJan_   CC_PopUMMay  
#169.462       -49.809        -9.052       -42.239       -30.706         2.119  
#CC_PopUMNov_  
#-62.376  

summary(TLUMmod)
#Residuals:
#Min      1Q  Median      3Q     Max 
#-72.581 -17.617   1.538  18.580 129.538 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   169.462      4.716  35.933  < 2e-16 ***
#  CC_PopUMDec   -49.809      6.411  -7.770 8.11e-14 ***
#  CC_PopUMFeb    -9.052      6.477  -1.398    0.163    
#CC_PopUMJan   -42.239      6.443  -6.555 1.90e-10 ***
#  CC_PopUMJan_  -30.706      6.443  -4.766 2.73e-06 ***
#  CC_PopUMMay     2.119      5.523   0.384    0.701    
#CC_PopUMNov_  -62.376      6.379  -9.778  < 2e-16 ***
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 29.45 on 364 degrees of freedom
#(7 observations deleted due to missingness)
#Multiple R-squared:  0.4114,	Adjusted R-squared:  0.4017 
#F-statistic:  42.4 on 6 and 364 DF,  p-value: < 2.2e-16

ggqqplot(residuals(TLUMmod))
#nice

shapiro_test(residuals(TLUMmod))
#  variable           statistic p.value
#<chr>                  <dbl>   <dbl>
#  1 residuals(TLUMmod)     0.987 0.00244

#test homogentiy of variance of the model as a whole 
plot(TLUMmod, 1)

check_model(TLUMmod, check = "all")

anova(TLUMmod)
#Response: TLeng_mm
#Df Sum Sq Mean Sq F value    Pr(>F)    
#CC_Pop      6 220650   36775  42.396 < 2.2e-16 ***
#Residuals 364 315738     867       

AIC(TLUMmod)
#3571.791

emmeans(TLUMmod, pairwise ~ CC_Pop, adjust = "tukey")
#$emmeans
#CC_Pop  emmean   SE  df lower.CL upper.CL
#UMApril    169 4.72 364    160.2      179
#UMDec      120 4.34 364    111.1      128
#UMFeb      160 4.44 364    151.7      169
#UMJan      127 4.39 364    118.6      136
#UMJan_     139 4.39 364    130.1      147
#UMMay      172 2.87 364    165.9      177
#UMNov_     107 4.30 364     98.6      116
#Confidence level used: 0.95 
#$contrasts
#contrast         estimate   SE  df t.ratio p.value
#UMApril - UMDec     49.81 6.41 364   7.770  <.0001
#UMApril - UMFeb      9.05 6.48 364   1.398  0.8030
#UMApril - UMJan     42.24 6.44 364   6.555  <.0001
#UMApril - UMJan_    30.71 6.44 364   4.766  0.0001
#UMApril - UMMay     -2.12 5.52 364  -0.384  0.9997
#UMApril - UMNov_    62.38 6.38 364   9.778  <.0001
#UMDec - UMFeb      -40.76 6.21 364  -6.563  <.0001
#UMDec - UMJan       -7.57 6.18 364  -1.226  0.8838
#UMDec - UMJan_     -19.10 6.18 364  -3.094  0.0344
#UMDec - UMMay      -51.93 5.21 364  -9.972  <.0001
#UMDec - UMNov_      12.57 6.11 364   2.057  0.3806
#UMFeb - UMJan       33.19 6.24 364   5.315  <.0001
#UMFeb - UMJan_      21.65 6.24 364   3.468  0.0104
#UMFeb - UMMay      -11.17 5.29 364  -2.112  0.3476
#UMFeb - UMNov_      53.32 6.18 364   8.631  <.0001
#UMJan - UMJan_     -11.53 6.21 364  -1.858  0.5100
#UMJan - UMMay      -44.36 5.25 364  -8.453  <.0001
#UMJan - UMNov_      20.14 6.14 364   3.278  0.0195
#UMJan_ - UMMay     -32.83 5.25 364  -6.255  <.0001
#UMJan_ - UMNov_     31.67 6.14 364   5.156  <.0001
#UMMay - UMNov_      64.50 5.17 364  12.478  <.0001

sjstats::anova_stats(car::Anova(TLUMmod, type = 3)) %>% dplyr::select(1:7)
#term      |     sumsq |    meansq |  df | statistic | p.value | etasq
#---------------------------------------------------------------------
#  CC_Pop    | 2.206e+05 | 36774.917 |   6 |    42.396 |  < .001 | 0.411
#Residuals | 3.157e+05 |   867.412 | 364 |           |         |      


####mass 

#check for outliers 
Massout<-UMonly %>% 
  group_by(CC_Pop) %>% 
  identify_outliers(mass_g)

#normality for each group of this measure 
UMonly %>% 
  group_by(CC_Pop) %>% 
  shapiro_test(mass_g)
#januarys and may non

ggqqplot(UMonly, "mass_g", facet.by = "CC_Pop")
#may has more variability because it has more samples; but not bad 

M_UMmod<- lm(mass_g ~ CC_Pop, data = UMonly)
M_UMmod
#Coefficients:
#(Intercept)   CC_PopUMDec   CC_PopUMFeb   CC_PopUMJan  CC_PopUMJan_   CC_PopUMMay  
#16.179       -10.722        -2.530        -9.683        -7.485         5.787  
#CC_PopUMNov_  
#-12.150  

summary(M_UMmod)
#Residuals:
#Min      1Q  Median      3Q     Max 
#-18.865  -3.720  -0.667   2.463  61.645 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    16.179      1.582  10.224  < 2e-16 ***
#  CC_PopUMDec   -10.722      2.151  -4.985 9.54e-07 ***
#  CC_PopUMFeb    -2.530      2.173  -1.164 0.245113    
#CC_PopUMJan    -9.683      2.162  -4.479 1.00e-05 ***
#  CC_PopUMJan_   -7.485      2.162  -3.462 0.000598 ***
#  CC_PopUMMay     5.787      1.837   3.149 0.001769 ** 
#  CC_PopUMNov_  -12.150      2.140  -5.676 2.78e-08 ***
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 9.882 on 371 degrees of freedom
#Multiple R-squared:  0.3384,	Adjusted R-squared:  0.3277 
#F-statistic: 31.62 on 6 and 371 DF,  p-value: < 2.2e-16

ggqqplot(residuals(M_UMmod))

shapiro_test(residuals(M_UMmod))

#test homogentiy of variance of the model as a whole 
plot(M_UMmod, 1)

check_model(M_UMmod, check = "all")
#model as a whole is not bad 

anova(M_UMmod)
#Response: mass_g
#Df Sum Sq Mean Sq F value    Pr(>F)    
#CC_Pop      6  18528 3088.01  31.622 < 2.2e-16 ***
#Residuals 371  36230   97.65        

AIC(M_UMmod)
#2813.435

emmeans(M_UMmod, pairwise ~ CC_Pop, adjust = "tukey")
#$emmeans
#CC_Pop  emmean    SE  df lower.CL upper.CL
#UMApril  16.18 1.582 371    13.07    19.29
#UMDec     5.46 1.457 371     2.59     8.32
#UMFeb    13.65 1.490 371    10.72    16.58
#UMJan     6.50 1.473 371     3.60     9.39
#UMJan_    8.69 1.473 371     5.80    11.59
#UMMay    21.97 0.934 371    20.13    23.80
#UMNov_    4.03 1.441 371     1.19     6.86
#Confidence level used: 0.95 
#$contrasts
#contrast         estimate   SE  df t.ratio p.value
#UMApril - UMDec     10.72 2.15 371   4.985  <.0001
#UMApril - UMFeb      2.53 2.17 371   1.164  0.9071
#UMApril - UMJan      9.68 2.16 371   4.479  0.0002
#UMApril - UMJan_     7.49 2.16 371   3.462  0.0106
#UMApril - UMMay     -5.79 1.84 371  -3.149  0.0291
#UMApril - UMNov_    12.15 2.14 371   5.676  <.0001
#UMDec - UMFeb       -8.19 2.08 371  -3.931  0.0019
#UMDec - UMJan       -1.04 2.07 371  -0.502  0.9988
#UMDec - UMJan_      -3.24 2.07 371  -1.562  0.7065
#UMDec - UMMay      -16.51 1.73 371  -9.540  <.0001
#UMDec - UMNov_       1.43 2.05 371   0.697  0.9927
#UMFeb - UMJan        7.15 2.10 371   3.414  0.0125
#UMFeb - UMJan_       4.96 2.10 371   2.365  0.2165
#UMFeb - UMMay       -8.32 1.76 371  -4.730  0.0001
#UMFeb - UMNov_       9.62 2.07 371   4.641  0.0001
#UMJan - UMJan_      -2.20 2.08 371  -1.055  0.9406
#UMJan - UMMay      -15.47 1.74 371  -8.869  <.0001
#UMJan - UMNov_       2.47 2.06 371   1.197  0.8950
#UMJan_ - UMMay     -13.27 1.74 371  -7.610  <.0001
#UMJan_ - UMNov_      4.66 2.06 371   2.263  0.2648
#UMMay - UMNov_      17.94 1.72 371  10.444  <.0001

sjstats::anova_stats(car::Anova(M_UMmod, type = 3)) %>% dplyr::select(1:7)

#term      |     sumsq |   meansq |  df | statistic | p.value | etasq
#--------------------------------------------------------------------
#CC_Pop    | 18528.083 | 3088.014 |   6 |    31.622 |  < .001 | 0.338
#Residuals | 36229.875 |   97.655 | 371 |           |         |  

###condition factor 
##convert to character 
data_K1<- transform(UMonly, K = as.numeric(K))

#check for outliers 
HBout<-data_K1 %>% 
  group_by(CC_Pop) %>% 
  identify_outliers(K)

#normality for each group of this measure 
data_K1 %>% 
  group_by(CC_Pop) %>% 
  shapiro_test(K)

ggqqplot(data_K1, "K", facet.by = "CC_Pop")
#all look good, UMMay with some tail 

UM_Kmod<- lm(K ~ CC_Pop, data = data_K1)
UM_Kmod
#Coefficients:
#(Intercept)   CC_PopUMDec   CC_PopUMFeb   CC_PopUMJan  CC_PopUMJan_   CC_PopUMMay  
#0.302028      0.004508      0.003530     -0.011327     -0.012578      0.037962  
#CC_PopUMNov_  
#0.009020 

summary(UM_Kmod)
#Residuals:
#Min        1Q    Median        3Q       Max 
#-0.207194 -0.021950 -0.004522  0.015947  0.279886 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   0.302028   0.007018  43.038  < 2e-16 ***
#  CC_PopUMDec   0.004508   0.009539   0.473    0.637    
#CC_PopUMFeb   0.003530   0.009638   0.366    0.714    
#CC_PopUMJan  -0.011327   0.009588  -1.181    0.238    
#CC_PopUMJan_ -0.012578   0.009588  -1.312    0.190    
#CC_PopUMMay   0.037962   0.008218   4.619 5.35e-06 ***
#  CC_PopUMNov_  0.009020   0.009493   0.950    0.343    
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 0.04383 on 364 degrees of freedom
#(7 observations deleted due to missingness)
#Multiple R-squared:  0.1595,	Adjusted R-squared:  0.1457 
#F-statistic: 11.51 on 6 and 364 DF,  p-value: 8.353e-12

ggqqplot(residuals(UM_Kmod))

shapiro_test(residuals(UM_Kmod))

#test homogentiy of variance of the model as a whole 
plot(UM_Kmod, 1)

check_model(UM_Kmod, check = "all")
#looks good 

anova(UM_Kmod)
#Response: K
#Df  Sum Sq   Mean Sq F value    Pr(>F)    
#CC_Pop      6 0.13268 0.0221129  11.513 8.353e-12 ***
#Residuals 364 0.69911 0.0019206                      
  
AIC(UM_Kmod)
#-1258.856

emmeans(UM_Kmod, pairwise ~ CC_Pop, adjust = "tukey")
#$emmeans
#CC_Pop  emmean      SE  df lower.CL upper.CL
#UMApril  0.302 0.00702 364    0.288    0.316
#UMDec    0.307 0.00646 364    0.294    0.319
#UMFeb    0.306 0.00661 364    0.293    0.319
#UMJan    0.291 0.00653 364    0.278    0.304
#UMJan_   0.289 0.00653 364    0.277    0.302
#UMMay    0.340 0.00428 364    0.332    0.348
#UMNov_   0.311 0.00639 364    0.298    0.324
#Confidence level used: 0.95 
#$contrasts
#contrast          estimate      SE  df t.ratio p.value
#UMApril - UMDec  -0.004508 0.00954 364  -0.473  0.9992
#UMApril - UMFeb  -0.003530 0.00964 364  -0.366  0.9998
#UMApril - UMJan   0.011327 0.00959 364   1.181  0.9009
#UMApril - UMJan_  0.012578 0.00959 364   1.312  0.8462
#UMApril - UMMay  -0.037962 0.00822 364  -4.619  0.0001
#UMApril - UMNov_ -0.009020 0.00949 364  -0.950  0.9639
#UMDec - UMFeb     0.000978 0.00924 364   0.106  1.0000
#UMDec - UMJan     0.015835 0.00919 364   1.723  0.6008
#UMDec - UMJan_    0.017086 0.00919 364   1.859  0.5088
#UMDec - UMMay    -0.033454 0.00775 364  -4.317  0.0004
#UMDec - UMNov_   -0.004512 0.00909 364  -0.496  0.9989
#UMFeb - UMJan     0.014857 0.00929 364   1.599  0.6831
#UMFeb - UMJan_    0.016108 0.00929 364   1.734  0.5939
#UMFeb - UMMay    -0.034432 0.00787 364  -4.375  0.0003
#UMFeb - UMNov_   -0.005490 0.00919 364  -0.597  0.9969
#UMJan - UMJan_    0.001251 0.00924 364   0.135  1.0000
#UMJan - UMMay    -0.049289 0.00781 364  -6.312  <.0001
#UMJan - UMNov_   -0.020347 0.00914 364  -2.226  0.2841
#UMJan_ - UMMay   -0.050540 0.00781 364  -6.472  <.0001
#UMJan_ - UMNov_  -0.021598 0.00914 364  -2.363  0.2176
#UMMay - UMNov_    0.028942 0.00769 364   3.763  0.0037

sjstats::anova_stats(car::Anova(UM_Kmod, type = 3)) %>%dplyr::select(1:7)
#term      | sumsq | meansq |  df | statistic | p.value | etasq
#--------------------------------------------------------------
#CC_Pop    | 0.133 |  0.022 |   6 |    11.513 |  < .001 | 0.160
#Residuals | 0.699 |  0.002 | 364 |           |         |  

###standard growth rate 
##this is the average of each group, not a direct comparison of SGR per individual over time because sampples were random 
##there is no value for november here because nothing to calculate from previous to this date measure 

#check for outliers 
SGRout<-UMonly %>% 
  group_by(CC_Pop) %>% 
  identify_outliers(SGR)

#normality for each group of this measure 
UMonly %>% 
  group_by(CC_Pop) %>% 
  shapiro_test(SGR)

ggqqplot(UMonly, "SGR", facet.by = "CC_Pop")

SGRmod<- lm(SGR ~ CC_Pop, data = UMonly)
SGRmod
#Coefficients:
#(Intercept)   CC_PopUMDec   CC_PopUMFeb   CC_PopUMJan  CC_PopUMJan_   CC_PopUMMay  
#-0.09109       1.65917       2.14007       0.52396       1.20632       0.14961 
summary(SGRmod)
#Residuals:
#Min       1Q   Median       3Q      Max 
#-10.3005  -1.5561   0.0418   1.4274  10.4507 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)   
#(Intercept)  -0.09109    0.53251  -0.171  0.86432   
#CC_PopUMDec   1.65917    0.68447   2.424  0.01611 * 
#  CC_PopUMFeb   2.14007    0.69059   3.099  0.00218 **
#  CC_PopUMJan   0.52396    0.68747   0.762  0.44674   
#CC_PopUMJan_  1.20632    0.68747   1.755  0.08061 . 
#CC_PopUMMay   0.14961    0.75309   0.199  0.84269   
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 
#Residual standard error: 2.917 on 234 degrees of freedom
#(138 observations deleted due to missingness)
#Multiple R-squared:  0.06584,	Adjusted R-squared:  0.04588 
#F-statistic: 3.298 on 5 and 234 DF,  p-value: 0.006736

ggqqplot(residuals(SGRmod))
#not bad 

shapiro_test(residuals(SGRmod))

#test homogentiy of variance of the model as a whole 
plot(SGRmod, 1)

check_model(SGRmod, check = "all")
#overall looks good 

anova(SGRmod)
#Response: SGR
#Df Sum Sq Mean Sq F value   Pr(>F)   
#CC_Pop      5  140.3 28.0597  3.2984 0.006736 **
#Residuals 234 1990.7  8.5071                    

AIC(SGRmod)
#1202.83

emmeans(SGRmod, pairwise ~ CC_Pop, adjust = "tukey")
#$emmeans
#CC_Pop   emmean    SE  df lower.CL upper.CL
#UMApril -0.0911 0.533 234   -1.140    0.958
#UMDec    1.5681 0.430 234    0.721    2.415
#UMFeb    2.0490 0.440 234    1.183    2.915
#UMJan    0.4329 0.435 234   -0.424    1.289
#UMJan_   1.1152 0.435 234    0.259    1.972
#UMMay    0.0585 0.533 234   -0.991    1.108
#Confidence level used: 0.95 
#$contrasts
#contrast         estimate    SE  df t.ratio p.value
#UMApril - UMDec    -1.659 0.684 234  -2.424  0.1522
#UMApril - UMFeb    -2.140 0.691 234  -3.099  0.0262
#UMApril - UMJan    -0.524 0.687 234  -0.762  0.9735
#UMApril - UMJan_   -1.206 0.687 234  -1.755  0.4972
#UMApril - UMMay    -0.150 0.753 234  -0.199  1.0000
#UMDec - UMFeb      -0.481 0.615 234  -0.782  0.9704
#UMDec - UMJan       1.135 0.612 234   1.856  0.4319
#UMDec - UMJan_      0.453 0.612 234   0.741  0.9766
#UMDec - UMMay       1.510 0.684 234   2.205  0.2392
#UMFeb - UMJan       1.616 0.618 234   2.613  0.0980
#UMFeb - UMJan_      0.934 0.618 234   1.510  0.6582
#UMFeb - UMMay       1.990 0.691 234   2.882  0.0488
#UMJan - UMJan_     -0.682 0.615 234  -1.110  0.8770
#UMJan - UMMay       0.374 0.687 234   0.545  0.9942
#UMJan_ - UMMay      1.057 0.687 234   1.537  0.6407

sjstats::anova_stats(car::Anova(SGRmod, type = 3)) %>% dplyr::select(1:7)
#term      |    sumsq | meansq |  df | statistic | p.value | etasq
#-----------------------------------------------------------------
#CC_Pop    |  140.299 | 28.060 |   5 |     3.298 |   0.007 | 0.066
#Residuals | 1990.656 |  8.507 | 234 |           |         | 


################################################
##account for interaction of time point as factor 

UM_LCC<-lm(avg_nuclei ~ TLeng_mm*CC_Pop, data = UMonly)
summary(UM_LCC)

#Residuals:
#  Min      1Q  Median      3Q     Max 
#-5.3168 -1.8603 -0.7492  0.9265 29.6312 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           56.508956   2.800951  20.175   <2e-16 ***
#  TLeng_mm              -0.010825   0.016071  -0.674    0.501    
#CC_PopUMDec           -0.417300   6.318470  -0.066    0.947    
#CC_PopUMFeb            3.602074   4.901821   0.735    0.463    
#CC_PopUMJan            0.485691   4.630200   0.105    0.917    
#CC_PopUMJan_           1.245024   4.600535   0.271    0.787    
#CC_PopUMMay            0.834712   3.563449   0.234    0.815    
#CC_PopUMNov_          -3.048057   4.919253  -0.620    0.536    
#TLeng_mm:CC_PopUMDec   0.013419   0.049734   0.270    0.787    
#TLeng_mm:CC_PopUMFeb  -0.020816   0.029537  -0.705    0.481    
#TLeng_mm:CC_PopUMJan   0.001938   0.032790   0.059    0.953    
#TLeng_mm:CC_PopUMJan_  0.004995   0.030509   0.164    0.870    
#TLeng_mm:CC_PopUMMay  -0.009798   0.020095  -0.488    0.626    
#TLeng_mm:CC_PopUMNov_  0.027140   0.040662   0.667    0.505    
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 4.086 on 342 degrees of freedom
#(22 observations deleted due to missingness)
#Multiple R-squared:  0.08649,	Adjusted R-squared:  0.05177 
#F-statistic: 2.491 on 13 and 342 DF,  p-value: 0.002917

anova(UM_LCC)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#TLeng_mm          1  268.7 268.658 16.0879 7.43e-05 ***
#CC_Pop            6  242.9  40.477  2.4238  0.02621 *  
#TLeng_mm:CC_Pop   6   29.2   4.867  0.2915  0.94083    
#Residuals       342 5711.2  16.699 

AIC(UM_LCC)
#2028.274


UM_MCC<-lm(avg_nuclei ~ mass_g*CC_Pop, data = UMonly)
summary(UM_MCC)
#Residuals:
#Min      1Q  Median      3Q     Max 
#-5.2564 -1.9096 -0.7802  0.8286 29.5116 
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)         55.511775   1.328574  41.783   <2e-16 ***
#  mass_g              -0.051753   0.071622  -0.723    0.470    
#CC_PopUMDec          0.453384   2.307328   0.196    0.844    
#CC_PopUMFeb          0.848065   2.036866   0.416    0.677    
#CC_PopUMJan          0.430595   1.890901   0.228    0.820    
#CC_PopUMJan_         1.709440   1.761777   0.970    0.333    
#CC_PopUMMay         -0.672444   1.522804  -0.442    0.659    
#CC_PopUMNov_        -0.963914   2.011784  -0.479    0.632    
#mass_g:CC_PopUMDec   0.131803   0.335604   0.393    0.695    
#mass_g:CC_PopUMFeb  -0.045274   0.126158  -0.359    0.720    
#mass_g:CC_PopUMJan   0.039686   0.198389   0.200    0.842    
#mass_g:CC_PopUMJan_  0.019972   0.134164   0.149    0.882    
#mass_g:CC_PopUMMay  -0.006114   0.076037  -0.080    0.936    
#mass_g:CC_PopUMNov_  0.215619   0.352352   0.612    0.541    
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 4.059 on 349 degrees of freedom
#(15 observations deleted due to missingness)
#Multiple R-squared:  0.1013,	Adjusted R-squared:  0.06778 
#F-statistic: 3.025 on 13 and 349 DF,  p-value: 0.0003059

anova(UM_MCC)
#Response: avg_nuclei
#Df Sum Sq Mean Sq F value    Pr(>F)    
#mass_g          1  416.8  416.80 25.3001 7.866e-07 ***
#CC_Pop          6  217.0   36.17  2.1956   0.04296 *  
#mass_g:CC_Pop   6   13.9    2.32  0.1409   0.99070    
#Residuals     349 5749.6   16.47  

AIC(UM_MCC)
#2062.928

UM_KCC<-lm(avg_nuclei ~ K*CC_Pop, data = data_K1)
summary(UM_KCC)
#Residuals:
#Min      1Q  Median      3Q     Max 
#-5.0798 -1.9055 -0.7066  0.7605 29.1847 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)      53.490      4.647  11.512   <2e-16 ***
#  K                 3.921     15.232   0.257    0.797    
#CC_PopUMDec     -11.378      9.354  -1.216    0.225    
#CC_PopUMFeb      -1.331      8.567  -0.155    0.877    
#CC_PopUMJan      -4.710      7.547  -0.624    0.533    
#CC_PopUMJan_     -2.205      6.931  -0.318    0.751    
#CC_PopUMMay       2.420      5.152   0.470    0.639    
#CC_PopUMNov_      1.687      8.403   0.201    0.841    
#K:CC_PopUMDec    42.695     30.488   1.400    0.162    
#K:CC_PopUMFeb     5.493     27.980   0.196    0.844    
#K:CC_PopUMJan    20.447     25.420   0.804    0.422    
#K:CC_PopUMJan_   15.631     23.309   0.671    0.503    
#K:CC_PopUMMay   -10.477     16.498  -0.635    0.526    
#K:CC_PopUMNov_   -3.821     27.111  -0.141    0.888    
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 4.076 on 342 degrees of freedom
#(22 observations deleted due to missingness)
#Multiple R-squared:  0.09101,	Adjusted R-squared:  0.05646 
#F-statistic: 2.634 on 13 and 342 DF,  p-value: 0.001617

anova(UM_KCC)
#Response: avg_nuclei
#Df Sum Sq Mean Sq F value    Pr(>F)    
#K           1   44.8  44.810  2.6967 0.1014784    
#CC_Pop      6  408.1  68.012  4.0930 0.0005556 ***
#K:CC_Pop    6  116.1  19.349  1.1644 0.3248813    
#Residuals 342 5682.9  16.617 

AIC(UM_KCC)
#2026.508

UM_SGRCC<-lm(avg_nuclei ~ SGR*CC_Pop, data = UMonly)
summary(UM_SGRCC)
#Residuals:
#Min      1Q  Median      3Q     Max 
#-5.6112 -1.9704 -0.9032  0.7749 29.2435 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       54.3987     0.8373  64.966   <2e-16 ***
#  SGR               -1.4130     1.3991  -1.010   0.3136    
#CC_PopUMDec        2.1874     1.1578   1.889   0.0601 .  
#CC_PopUMFeb        1.0967     1.2452   0.881   0.3794    
#CC_PopUMJan        1.4875     1.0807   1.376   0.1700    
#CC_PopUMJan_       2.5128     1.1163   2.251   0.0253 *  
#  CC_PopUMMay       -0.6283     1.2016  -0.523   0.6016    
#SGR:CC_PopUMDec    1.2956     1.4268   0.908   0.3648    
#SGR:CC_PopUMFeb    1.1886     1.4312   0.830   0.4072    
#SGR:CC_PopUMJan    1.3617     1.4182   0.960   0.3380    
#SGR:CC_PopUMJan_   1.4430     1.4242   1.013   0.3121    
#SGR:CC_PopUMMay    1.4272     1.4096   1.012   0.3124    
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 4.533 on 226 degrees of freedom
#(140 observations deleted due to missingness)
#Multiple R-squared:  0.05933,	Adjusted R-squared:  0.01355 
#F-statistic: 1.296 on 11 and 226 DF,  p-value: 0.2277

anova(UM_SGRCC)
#Response: avg_nuclei
#Df Sum Sq Mean Sq F value  Pr(>F)  
#SGR          1    1.0   0.984  0.0479 0.82697  
#CC_Pop       5  259.7  51.932  2.5275 0.02996 *
#  SGR:CC_Pop   5   32.3   6.453  0.3141 0.90422  
#Residuals  226 4643.6  20.547 

AIC(UM_SGRCC)
#1408.507

