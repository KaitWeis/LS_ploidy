######################################
#
#BLOOD SMEAR ANALYSIS   
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

###################
##starting with cell values 

###cell volume um3####

#check for outliers 
CVout<-Cells1 %>% 
  group_by(Smear_Pop) %>% 
  identify_outliers(avg_volume_um3)

#normality for each group of this measure 
Cells1 %>% 
  group_by(Smear_Pop) %>% 
  shapiro_test(avg_volume_um3)
##normal 

ggqqplot(Cells1, "avg_volume_um3", facet.by = "Smear_Pop")

#levens for homogentiy of variance 
Cells1 %>% levene_test(avg_volume_um3 ~ Smear_Pop)
#homogenous 

CVmod<- lm(avg_volume_um3 ~ Smear_Pop, data = Cells1)
CVmod
#Coefficients:
#(Intercept)     Smear_PopLR     Smear_PopNR    Smear_PopPDB  Smear_PopUMFeb  
#7451.19           35.75        -1380.78         -880.86          -38.41  
#Smear_PopUMMay  
#741.49  

summary(CVmod)
#Residuals:
#Min      1Q  Median      3Q     Max 
#-3040.8  -668.6  -157.0   648.7  6327.3 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     7451.19     359.67  20.717  < 2e-16 ***
#  Smear_PopLR       35.75     473.53   0.075  0.93993    
#Smear_PopNR    -1380.78     420.47  -3.284  0.00131 ** 
#  Smear_PopPDB    -880.86     418.65  -2.104  0.03727 *  
#  Smear_PopUMFeb   -38.41     440.50  -0.087  0.93065    
#Smear_PopUMMay   741.49     422.41   1.755  0.08151 .  
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 1193 on 132 degrees of freedom
#Multiple R-squared:  0.301,	Adjusted R-squared:  0.2745 
#F-statistic: 11.37 on 5 and 132 DF,  p-value: 3.935e-09

ggqqplot(residuals(CVmod))
#looks good 

shapiro_test(residuals(CVmod))

#test homogentiy of variance of the model as a whole 
plot(CVmod, 1)

check_model(CVmod, check = "all")

anova(CVmod)
#Response: avg_volume_um3
#Df    Sum Sq  Mean Sq F value    Pr(>F)    
#Smear_Pop   5  80871859 16174372  11.367 3.935e-09 ***
#Residuals 132 187834912  1422992    

AIC(CVmod)
#2354.714

emmeans(CVmod, pairwise ~ Smear_Pop, adjust = "tukey")
#$emmeans
#Smear_Pop emmean  SE  df lower.CL upper.CL
#BR          7451 360 132     6740     8163
#LR          7487 308 132     6878     8096
#NR          6070 218 132     5640     6501
#PDB         6570 214 132     6147     6994
#UMFeb       7413 254 132     6910     7916
#UMMay       8193 222 132     7755     8631
#Confidence level used: 0.95 
#$contrasts
#contrast      estimate  SE  df t.ratio p.value
#BR - LR          -35.8 474 132  -0.075  1.0000
#BR - NR         1380.8 420 132   3.284  0.0162
#BR - PDB         880.9 419 132   2.104  0.2916
#BR - UMFeb        38.4 441 132   0.087  1.0000
#BR - UMMay      -741.5 422 132  -1.755  0.4981
#LR - NR         1416.5 377 132   3.755  0.0035
#LR - PDB         916.6 375 132   2.443  0.1493
#LR - UMFeb        74.2 399 132   0.186  1.0000
#LR - UMMay      -705.7 379 132  -1.860  0.4313
#NR - PDB        -499.9 306 132  -1.636  0.5761
#NR - UMFeb     -1342.4 335 132  -4.009  0.0014
#NR - UMMay     -2122.3 311 132  -6.832  <.0001
#PDB - UMFeb     -842.5 333 132  -2.533  0.1220
#PDB - UMMay    -1622.4 308 132  -5.264  <.0001
#UMFeb - UMMay   -779.9 337 132  -2.312  0.1965

sjstats::anova_stats(car::Anova(CVmod, type = 3)) %>% dplyr::select(1:7)
#term      |     sumsq |    meansq |  df | statistic | p.value | etasq
#---------------------------------------------------------------------
#  Smear_Pop | 8.087e+07 | 1.617e+07 |   5 |    11.366 |  < .001 | 0.301
#Residuals | 1.878e+08 | 1.423e+06 | 132 |           |         | 

########################################
### Cell surface area### 
#check for outliers 
CSout<-Cells1 %>% 
  group_by(Smear_Pop) %>% 
  identify_outliers(avg_surface_um2)

#normality for each group of this measure 
Cells1 %>% 
  group_by(Smear_Pop) %>% 
  shapiro_test(avg_surface_um2)
#normal 

ggqqplot(Cells1, "avg_surface_um2", facet.by = "Smear_Pop")

#levens for homogentiy of variance 
Cells1 %>% levene_test(avg_surface_um2 ~ Smear_Pop)

CSmod<- lm(avg_surface_um2 ~ Smear_Pop, data = Cells1)
CSmod
#Coefficients:
#(Intercept)     Smear_PopLR     Smear_PopNR    Smear_PopPDB  Smear_PopUMFeb  
#513.8601          0.7678        -59.8408        -40.4624          6.3430  
#Smear_PopUMMay  
#52.6226  

summary(CSmod)
#Residuals:
#Min       1Q   Median       3Q      Max 
#-141.557  -34.139   -4.646   26.929  270.664 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    513.8601    16.6126  30.932  < 2e-16 ***
#  Smear_PopLR      0.7678    21.8714   0.035  0.97205    
#Smear_PopNR    -59.8408    19.4208  -3.081  0.00251 ** 
#  Smear_PopPDB   -40.4624    19.3366  -2.093  0.03831 *  
#  Smear_PopUMFeb   6.3430    20.3461   0.312  0.75572    
#Smear_PopUMMay  52.6226    19.5105   2.697  0.00791 ** 
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 55.1 on 132 degrees of freedom
#Multiple R-squared:  0.3603,	Adjusted R-squared:  0.3361 
#F-statistic: 14.87 on 5 and 132 DF,  p-value: 1.456e-11

ggqqplot(residuals(CSmod))

shapiro_test(residuals(CSmod))

#test homogentiy of variance of the model as a whole 
plot(CSmod, 1)

check_model(CSmod, check = "all")

anova(CSmod)
#Response: avg_surface_um2
#Df Sum Sq Mean Sq F value    Pr(>F)    
#Smear_Pop   5 225731   45146  14.871 1.456e-11 ***
#  Residuals 132 400719    3036   
AIC(CSmod)
#1506.006

emmeans(CSmod, pairwise ~ Smear_Pop, adjust = "tukey")
#$emmeans
#Smear_Pop emmean   SE  df lower.CL upper.CL
#BR           514 16.6 132      481      547
#LR           515 14.2 132      486      543
#NR           454 10.1 132      434      474
#PDB          473  9.9 132      454      493
#UMFeb        520 11.7 132      497      543
#UMMay        566 10.2 132      546      587
#Confidence level used: 0.95 
#$contrasts
#contrast      estimate   SE  df t.ratio p.value
#BR - LR         -0.768 21.9 132  -0.035  1.0000
#BR - NR         59.841 19.4 132   3.081  0.0296
#BR - PDB        40.462 19.3 132   2.093  0.2976
#BR - UMFeb      -6.343 20.3 132  -0.312  0.9996
#BR - UMMay     -52.623 19.5 132  -2.697  0.0826
#LR - NR         60.609 17.4 132   3.479  0.0088
#LR - PDB        41.230 17.3 132   2.379  0.1712
#LR - UMFeb      -5.575 18.4 132  -0.302  0.9997
#LR - UMMay     -51.855 17.5 132  -2.959  0.0417
#NR - PDB       -19.378 14.1 132  -1.373  0.7429
#NR - UMFeb     -66.184 15.5 132  -4.279  0.0005
#NR - UMMay    -112.463 14.3 132  -7.838  <.0001
#PDB - UMFeb    -46.805 15.4 132  -3.047  0.0326
#PDB - UMMay    -93.085 14.2 132  -6.540  <.0001
#UMFeb - UMMay  -46.280 15.6 132  -2.971  0.0403

sjstats::anova_stats(car::Anova(CSmod, type = 3)) %>% dplyr::select(1:7)

#term      |     sumsq |    meansq |  df | statistic | p.value | etasq
#---------------------------------------------------------------------
# Smear_Pop | 2.257e+05 | 45146.176 |   5 |    14.872 |  < .001 | 0.360
#Residuals | 4.007e+05 |  3035.747 | 132 |           |         | 

#########################################
##Cell surface area to volume ratio####

#check for outliers 
HBout<-Cells1 %>% 
  group_by(Smear_Pop) %>% 
  identify_outliers(SAVratio)

#normality for each group of this measure 
Cells1 %>% 
  group_by(Smear_Pop) %>% 
  shapiro_test(SAVratio)

ggqqplot(Cells1, "SAVratio", facet.by = "Smear_Pop")

#levens for homogentiy of variance 
Cells1 %>% levene_test(SAVratio ~ Smear_Pop)

CSAVmod<- lm(SAVratio ~ Smear_Pop, data = Cells1)
CSAVmod
#Coefficients:
#(Intercept)     Smear_PopLR     Smear_PopNR    Smear_PopPDB  Smear_PopUMFeb  
#0.0692697      -0.0001117       0.0059969       0.0036267       0.0012792  
#Smear_PopUMMay  
#0.0006975 

summary(CSAVmod)
#Residuals:
#Min         1Q     Median         3Q        Max 
#-0.0152068 -0.0025540 -0.0002141  0.0023367  0.0134094 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     0.0692697  0.0011967  57.884  < 2e-16 ***
#  Smear_PopLR    -0.0001117  0.0015755  -0.071   0.9436    
#Smear_PopNR     0.0059969  0.0013990   4.287 3.47e-05 ***
#  Smear_PopPDB    0.0036267  0.0013929   2.604   0.0103 *  
#  Smear_PopUMFeb  0.0012792  0.0014656   0.873   0.3844    
#Smear_PopUMMay  0.0006975  0.0014054   0.496   0.6205    
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 0.003969 on 132 degrees of freedom
#Multiple R-squared:  0.2529,	Adjusted R-squared:  0.2246 
#F-statistic: 8.938 on 5 and 132 DF,  p-value: 2.468e-07

ggqqplot(residuals(CSAVmod))

shapiro_test(residuals(CSAVmod))

#test homogentiy of variance of the model as a whole 
plot(CSAVmod, 1)

check_model(CSAVmod, check = "all")

anova(CSAVmod)
#Response: SAVratio
#Df     Sum Sq    Mean Sq F value    Pr(>F)    
#Smear_Pop   5 0.00070402 1.4080e-04  8.9383 2.468e-07 ***
#Residuals 132 0.00207937 1.5753e-05  

AIC(CSAVmod)

emmeans(CSAVmod, pairwise ~ Smear_Pop, adjust = "tukey")
#$emmeans
#Smear_Pop emmean       SE  df lower.CL upper.CL
#BR        0.0693 0.001197 132   0.0669   0.0716
#LR        0.0692 0.001025 132   0.0671   0.0712
#NR        0.0753 0.000725 132   0.0738   0.0767
#PDB       0.0729 0.000713 132   0.0715   0.0743
#UMFeb     0.0705 0.000846 132   0.0689   0.0722
#UMMay     0.0700 0.000737 132   0.0685   0.0714
#Confidence level used: 0.95 
#$contrasts
#contrast       estimate      SE  df t.ratio p.value
#BR - LR        0.000112 0.00158 132   0.071  1.0000
#BR - NR       -0.005997 0.00140 132  -4.287  0.0005
#BR - PDB      -0.003627 0.00139 132  -2.604  0.1036
#BR - UMFeb    -0.001279 0.00147 132  -0.873  0.9524
#BR - UMMay    -0.000697 0.00141 132  -0.496  0.9962
#LR - NR       -0.006109 0.00126 132  -4.867  <.0001
#LR - PDB      -0.003738 0.00125 132  -2.995  0.0378
#LR - UMFeb    -0.001391 0.00133 132  -1.047  0.9011
#LR - UMMay    -0.000809 0.00126 132  -0.641  0.9876
#NR - PDB       0.002370 0.00102 132   2.332  0.1889
#NR - UMFeb     0.004718 0.00111 132   4.235  0.0006
#NR - UMMay     0.005299 0.00103 132   5.127  <.0001
#PDB - UMFeb    0.002348 0.00111 132   2.122  0.2826
#PDB - UMMay    0.002929 0.00103 132   2.857  0.0549
#UMFeb - UMMay  0.000582 0.00112 132   0.518  0.9954

sjstats::anova_stats(car::Anova(CSAVmod, type = 3)) %>% dplyr::select(1:7)
#term      | sumsq | meansq |  df | statistic | p.value | etasq
#--------------------------------------------------------------
#Smear_Pop | 0.001 |      0 |   5 |     8.938 |  < .001 | 0.253
#Residuals | 0.002 |      0 | 132 |           |         |    


##################################################################
##Nuclei volume ##
#check for outliers 
HBout<-Nuc %>% 
  group_by(Smear_Pop) %>% 
  identify_outliers(avg_volume_um3)

#normality for each group of this measure 
Nuc %>% 
  group_by(Smear_Pop) %>% 
  shapiro_test(avg_volume_um3)

ggqqplot(Nuc, "avg_volume_um3", facet.by = "Smear_Pop")

#levens for homogentiy of variance 
Nuc %>% levene_test(avg_volume_um3 ~ Smear_Pop)

NVmod<- lm(avg_volume_um3 ~ Smear_Pop, data = Nuc)
NVmod
#   (Intercept)     Smear_PopLR     Smear_PopNR    Smear_PopPDB  Smear_PopUMFeb  
#655.69          -57.35         -228.94         -184.56           63.93  
#Smear_PopUMMay  
#-2.82 

summary(NVmod)
#Residuals:
#Min      1Q  Median      3Q     Max 
#-249.08  -73.18  -11.28   48.01 1073.65 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)      655.70      46.33  14.154  < 2e-16 ***
#  Smear_PopLR      -57.35      60.99  -0.940 0.348811    
#Smear_PopNR     -228.94      54.16  -4.227 4.38e-05 ***
#  Smear_PopPDB    -184.56      53.92  -3.423 0.000826 ***
#  Smear_PopUMFeb    63.92      56.74   1.127 0.261929    
#Smear_PopUMMay    -2.82      54.41  -0.052 0.958735    
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 153.6 on 132 degrees of freedom
#Multiple R-squared:  0.3533,	Adjusted R-squared:  0.3288 
#F-statistic: 14.42 on 5 and 132 DF,  p-value: 2.911e-11

ggqqplot(residuals(NVmod))

shapiro_test(residuals(NVmod))

#test homogentiy of variance of the model as a whole 
plot(NVmod, 1)

check_model(NVmod, check = "all")

anova(NVmod)
#Response: avg_volume_um3
#Df  Sum Sq Mean Sq F value    Pr(>F)    
#Smear_Pop   5 1702485  340497  14.423 2.911e-11 ***
#Residuals 132 3116264   23608         
AIC(NVmod)
#1789.062

emmeans(NVmod, pairwise ~ Smear_Pop, adjust = "tukey")

#$emmeans
#Smear_Pop emmean   SE  df lower.CL upper.CL
#BR           656 46.3 132      564      747
#LR           598 39.7 132      520      677
#NR           427 28.1 132      371      482
#PDB          471 27.6 132      417      526
#UMFeb        720 32.8 132      655      784
#UMMay        653 28.5 132      596      709
#Confidence level used: 0.95 
#$contrasts
#contrast      estimate   SE  df t.ratio p.value
#BR - LR          57.35 61.0 132   0.940  0.9353
#BR - NR         228.94 54.2 132   4.227  0.0006
#BR - PDB        184.56 53.9 132   3.423  0.0105
#BR - UMFeb      -63.93 56.7 132  -1.127  0.8695
#BR - UMMay        2.82 54.4 132   0.052  1.0000
#LR - NR         171.60 48.6 132   3.532  0.0074
#LR - PDB        127.21 48.3 132   2.632  0.0967
#LR - UMFeb     -121.27 51.4 132  -2.357  0.1792
#LR - UMMay      -54.53 48.9 132  -1.116  0.8741
#NR - PDB        -44.38 39.4 132  -1.128  0.8690
#NR - UMFeb     -292.87 43.1 132  -6.791  <.0001
#NR - UMMay     -226.12 40.0 132  -5.651  <.0001
#PDB - UMFeb    -248.49 42.8 132  -5.801  <.0001
#PDB - UMMay    -181.74 39.7 132  -4.579  0.0002
#UMFeb - UMMay    66.75 43.4 132   1.536  0.6415


sjstats::anova_stats(car::Anova(NVmod, type = 3)) %>% dplyr::select(1:7)
#term      |     sumsq |    meansq |  df | statistic | p.value | etasq
#---------------------------------------------------------------------
#Smear_Pop | 1.702e+06 | 3.405e+05 |   5 |    14.423 |  < .001 | 0.353
#Residuals | 3.116e+06 | 23608.059 | 132 |           |         |  


########################################
##Nuclei surface area##

#check for outliers 
HBout<-Nuc %>% 
  group_by(Smear_Pop) %>% 
  identify_outliers(avg_surface_um2)

#normality for each group of this measure 
Nuc %>% 
  group_by(Smear_Pop) %>% 
  shapiro_test(avg_surface_um2)

ggqqplot(Nuc, "avg_surface_um2", facet.by = "Smear_Pop")

#levens for homogentiy of variance 
Nuc %>% levene_test(avg_surface_um2 ~ Smear_Pop)

NSmod<- lm(avg_surface_um2 ~ Smear_Pop, data = Nuc)
NSmod
#Coefficients:
#(Intercept)     Smear_PopLR     Smear_PopNR    Smear_PopPDB  Smear_PopUMFeb  
#102.015          -5.918         -24.848         -23.056           7.247  
#Smear_PopUMMay  
#4.818  

summary(NSmod)
#Residuals:
#Min      1Q  Median      3Q     Max 
#-26.425  -7.779  -1.319   5.956  78.830 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     102.015      4.345  23.481  < 2e-16 ***
#  Smear_PopLR      -5.918      5.720  -1.035    0.303    
#Smear_PopNR     -24.848      5.079  -4.892 2.85e-06 ***
#  Smear_PopPDB    -23.056      5.057  -4.559 1.16e-05 ***
#  Smear_PopUMFeb    7.247      5.321   1.362    0.176    
#Smear_PopUMMay    4.818      5.102   0.944    0.347    
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 14.41 on 132 degrees of freedom
#Multiple R-squared:  0.4872,	Adjusted R-squared:  0.4678 
#F-statistic: 25.08 on 5 and 132 DF,  p-value: < 2.2e-16

ggqqplot(residuals(NSmod))

shapiro_test(residuals(NSmod))

#test homogentiy of variance of the model as a whole 
plot(NSmod, 1)

check_model(NSmod, check = "all")

anova(NSmod)
#Response: avg_surface_um2
#Df Sum Sq Mean Sq F value    Pr(>F)    
#Smear_Pop   5  26036  5207.3   25.08 < 2.2e-16 ***
#Residuals 132  27406   207.6 

AIC(NSmod)
#1135.823

emmeans(NSmod, pairwise ~ Smear_Pop, adjust = "tukey")
#$emmeans
#Smear_Pop emmean   SE  df lower.CL upper.CL
#BR         102.0 4.34 132     93.4    110.6
#LR          96.1 3.72 132     88.7    103.5
#NR          77.2 2.63 132     72.0     82.4
#PDB         79.0 2.59 132     73.8     84.1
#UMFeb      109.3 3.07 132    103.2    115.3
#UMMay      106.8 2.68 132    101.5    112.1
#Confidence level used: 0.95 
#$contrasts
#contrast      estimate   SE  df t.ratio p.value
#BR - LR           5.92 5.72 132   1.035  0.9054
#BR - NR          24.85 5.08 132   4.892  <.0001
#BR - PDB         23.06 5.06 132   4.559  0.0002
#BR - UMFeb       -7.25 5.32 132  -1.362  0.7496
#BR - UMMay       -4.82 5.10 132  -0.944  0.9342
#LR - NR          18.93 4.56 132   4.154  0.0008
#LR - PDB         17.14 4.53 132   3.782  0.0032
#LR - UMFeb      -13.16 4.82 132  -2.728  0.0764
#LR - UMMay      -10.74 4.58 132  -2.343  0.1847
#NR - PDB         -1.79 3.69 132  -0.485  0.9966
#NR - UMFeb      -32.09 4.04 132  -7.935  <.0001
#NR - UMMay      -29.67 3.75 132  -7.906  <.0001
#PDB - UMFeb     -30.30 4.02 132  -7.544  <.0001
#PDB - UMMay     -27.87 3.72 132  -7.488  <.0001
#UMFeb - UMMay     2.43 4.07 132   0.596  0.9911

sjstats::anova_stats(car::Anova(NSmod, type = 3)) %>% dplyr::select(1:7)
#term      |     sumsq |   meansq |  df | statistic | p.value | etasq
#--------------------------------------------------------------------
#Smear_Pop | 26036.297 | 5207.259 |   5 |    25.080 |  < .001 | 0.487
#Residuals | 27406.286 |  207.623 | 132 |           |   


###########################################
##Nuclei surface area to volume ratio##

#check for outliers 
HBout<-Nuc %>% 
  group_by(Smear_Pop) %>% 
  identify_outliers(SAVratio)

#normality for each group of this measure 
Nuc %>% 
  group_by(Smear_Pop) %>% 
  shapiro_test(SAVratio)

ggqqplot(Nuc, "SAVratio", facet.by = "Smear_Pop")

#levens for homogentiy of variance 
Nuc %>% levene_test(SAVratio ~ Smear_Pop)

NSAVmod<- lm(SAVratio ~ Smear_Pop, data = Nuc)
NSAVmod
#Coefficients:
#(Intercept)     Smear_PopLR     Smear_PopNR    Smear_PopPDB  Smear_PopUMFeb  
#0.156324        0.005824        0.025763        0.015052        0.001370  
#Smear_PopUMMay  
#0.009950 

summary(NSAVmod)
#Residuals:
#Min        1Q    Median        3Q       Max 
#-0.052806 -0.006236 -0.000202  0.007915  0.032435 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    0.156324   0.003886  40.226  < 2e-16 ***
#  Smear_PopLR    0.005824   0.005116   1.138  0.25707    
#Smear_PopNR    0.025763   0.004543   5.671 8.57e-08 ***
#  Smear_PopPDB   0.015052   0.004523   3.328  0.00113 ** 
#  Smear_PopUMFeb 0.001370   0.004759   0.288  0.77393    
#Smear_PopUMMay 0.009950   0.004564   2.180  0.03102 *  
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 0.01289 on 132 degrees of freedom
#Multiple R-squared:  0.3286,	Adjusted R-squared:  0.3032 
#F-statistic: 12.92 on 5 and 132 DF,  p-value: 3.1e-10

ggqqplot(residuals(NSAVmod))

shapiro_test(residuals(NSAVmod))

#test homogentiy of variance of the model as a whole 
plot(NSAVmod, 1)

check_model(NSAVmod, check = "all")

anova(NSAVmod)
#Response: SAVratio
#Df   Sum Sq    Mean Sq F value  Pr(>F)    
#Smear_Pop   5 0.010734 0.00214690  12.924 3.1e-10 ***
#Residuals 132 0.021928 0.00016612   

AIC(NSAVmod)
#-801.4931

emmeans(NSAVmod, pairwise ~ Smear_Pop, adjust = "tukey")
#$emmeans
#Smear_Pop emmean      SE  df lower.CL upper.CL
#BR         0.156 0.00389 132    0.149    0.164
#LR         0.162 0.00333 132    0.156    0.169
#NR         0.182 0.00235 132    0.177    0.187
#PDB        0.171 0.00231 132    0.167    0.176
#UMFeb      0.158 0.00275 132    0.152    0.163
#UMMay      0.166 0.00239 132    0.162    0.171
#Confidence level used: 0.95 
#$contrasts
#contrast      estimate      SE  df t.ratio p.value
#BR - LR       -0.00582 0.00512 132  -1.138  0.8645
#BR - NR       -0.02576 0.00454 132  -5.671  <.0001
#BR - PDB      -0.01505 0.00452 132  -3.328  0.0141
#BR - UMFeb    -0.00137 0.00476 132  -0.288  0.9997
#BR - UMMay    -0.00995 0.00456 132  -2.180  0.2541
#LR - NR       -0.01994 0.00408 132  -4.892  <.0001
#LR - PDB      -0.00923 0.00405 132  -2.276  0.2112
#LR - UMFeb     0.00445 0.00432 132   1.032  0.9064
#LR - UMMay    -0.00413 0.00410 132  -1.007  0.9150
#NR - PDB       0.01071 0.00330 132   3.245  0.0182
#NR - UMFeb     0.02439 0.00362 132   6.743  <.0001
#NR - UMMay     0.01581 0.00336 132   4.711  0.0001
#PDB - UMFeb    0.01368 0.00359 132   3.808  0.0029
#PDB - UMMay    0.00510 0.00333 132   1.532  0.6442
#UMFeb - UMMay -0.00858 0.00364 132  -2.355  0.1802

sjstats::anova_stats(car::Anova(NSAVmod, type = 3)) %>% dplyr::select(1:7)
#term      | sumsq | meansq |  df | statistic | p.value | etasq
#--------------------------------------------------------------
#  Smear_Pop | 0.011 |  0.002 |   5 |    12.924 |  < .001 | 0.329
#Residuals | 0.022 |  0.000 | 132 |           |         | 

