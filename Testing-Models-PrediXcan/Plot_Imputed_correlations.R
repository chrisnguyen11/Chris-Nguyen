cat("Compare Imputed Spearman Correlation")

library(data.table)
library(dplyr)
library(ggplot2)
#library(viridis)
library(tidyr)

"%&%" = function(a,b) paste(a,b,sep="")

## Global Variables
tiss="Tcell"
dir="/home/chris/topmed/expression/baseline_elasticnet/imputed_results/"

## Read in spearman correlations

#Correlations calculated with `genetic_pi1.R` script

### AFA
AFAmodelAFAExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_AFA_R2_-1_model_AFA_expression_spearman_correlation.txt",header=T)
ALLmodelAFAExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_ALL_R2_-1_model_AFA_expression_spearman_correlation.txt",header=T)
#CHNmodelAFAExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_CHN_R2_-1_model_AFA_expression_spearman_correlation.txt",header=T)
CAUmodelAFAExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_CAU_R2_-1_model_AFA_expression_spearman_correlation.txt",header=T)
HISmodelAFAExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_HIS_R2_-1_model_AFA_expression_spearman_correlation.txt",header=T)

### ALL
AFAmodelALLExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_AFA_R2_-1_model_ALL_expression_spearman_correlation.txt",header=T)
ALLmodelALLExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_ALL_R2_-1_model_ALL_expression_spearman_correlation.txt",header=T)
#CHNmodelALLExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_CHN_R2_-1_model_ALL_expression_spearman_correlation.txt",header=T)
CAUmodelALLExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_CAU_R2_-1_model_ALL_expression_spearman_correlation.txt",header=T)
HISmodelALLExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_HIS_R2_-1_model_ALL_expression_spearman_correlation.txt",header=T)

### CAU
AFAmodelCAUExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_AFA_R2_-1_model_CAU_expression_spearman_correlation.txt",header=T)
ALLmodelCAUExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_ALL_R2_-1_model_CAU_expression_spearman_correlation.txt",header=T)
#CHNmodelCAUExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_CHN_R2_-1_model_CAU_expression_spearman_correlation.txt",header=T)
CAUmodelCAUExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_CAU_R2_-1_model_CAU_expression_spearman_correlation.txt",header=T)
HISmodelCAUExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_HIS_R2_-1_model_CAU_expression_spearman_correlation.txt",header=T)

### CHN
#AFAmodelCHNExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_AFA_R2_-1_model_CHN_expression_spearman_correlation.txt",header=T)
#ALLmodelCHNExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_ALL_R2_-1_model_CHN_expression_spearman_correlation.txt",header=T)
#CHNmodelCHNExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_CHN_R2_-1_model_CHN_expression_spearman_correlation.txt",header=T)
#CAUmodelCHNExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_CAU_R2_-1_model_CHN_expression_spearman_correlation.txt",header=T)
#HISmodelCHNExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_HIS_R2_-1_model_CHN_expression_spearman_correlation.txt",header=T)

### HIS
AFAmodelHISExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_AFA_R2_-1_model_HIS_expression_spearman_correlation.txt",header=T)
ALLmodelHISExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_ALL_R2_-1_model_HIS_expression_spearman_correlation.txt",header=T)
#CHNmodelHISExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_CHN_R2_-1_model_HIS_expression_spearman_correlation.txt",header=T)
CAUmodelHISExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_CAU_R2_-1_model_HIS_expression_spearman_correlation.txt",header=T)
HISmodelHISExp<-fread("" %&% dir %&% "" %&% tiss %&% "/" %&% tiss %&% "_HIS_R2_-1_model_HIS_expression_spearman_correlation.txt",header=T)

## Concat all models for

### AFA
AFAmodelAFAExp<-AFAmodelAFAExp %>% mutate(model_pop="AFA",ID="AFAfilt")
ALLmodelAFAExp<-ALLmodelAFAExp %>% mutate(model_pop="ALL",ID="ALLfilt")
#CHNmodelAFAExp<-CHNmodelAFAExp %>% mutate(model_pop="CHN",ID="CHNfilt")
CAUmodelAFAExp<-CAUmodelAFAExp %>% mutate(model_pop="CAU",ID="CAUfilt")
HISmodelAFAExp<-HISmodelAFAExp %>% mutate(model_pop="HIS",ID="HISfilt")

AFA_Exp_Correlations<-bind_rows(AFAmodelAFAExp,
                                ALLmodelAFAExp,
#                                CHNmodelAFAExp,
                                CAUmodelAFAExp,
                                HISmodelAFAExp)

### ALL 
AFAmodelALLExp<-AFAmodelALLExp %>% mutate(model_pop="AFA",ID="AFAfilt")
ALLmodelALLExp<-ALLmodelALLExp %>% mutate(model_pop="ALL",ID="ALLfilt")
#CHNmodelALLExp<-CHNmodelALLExp %>% mutate(model_pop="CHN",ID="CHNfilt")
CAUmodelALLExp<-CAUmodelALLExp %>% mutate(model_pop="CAU",ID="CAUfilt")
HISmodelALLExp<-HISmodelALLExp %>% mutate(model_pop="HIS",ID="HISfilt")

ALL_Exp_Correlations<-bind_rows(AFAmodelALLExp,
                                ALLmodelALLExp,
#                                CHNmodelALLExp,
                                CAUmodelALLExp,
                                HISmodelALLExp)
### CHN
#AFAmodelCHNExp<-AFAmodelCHNExp %>% mutate(model_pop="AFA",ID="AFAfilt")
#ALLmodelCHNExp<-ALLmodelCHNExp %>% mutate(model_pop="ALL",ID="ALLfilt")
#CHNmodelCHNExp<-CHNmodelCHNExp %>% mutate(model_pop="CHN",ID="CHNfilt")
#CAUmodelCHNExp<-CAUmodelCHNExp %>% mutate(model_pop="CAU",ID="CAUfilt")
#HISmodelCHNExp<-HISmodelCHNExp %>% mutate(model_pop="HIS",ID="HISfilt")

#CHN_Exp_Correlations<-bind_rows(AFAmodelCHNExp,
#                                ALLmodelCHNExp,
#                                CHNmodelCHNExp,
#                                CAUmodelCHNExp,
#                                HISmodelCHNExp)

### CAU
AFAmodelCAUExp<-AFAmodelCAUExp %>% mutate(model_pop="AFA",ID="AFAfilt")
ALLmodelCAUExp<-ALLmodelCAUExp %>% mutate(model_pop="ALL",ID="ALLfilt")
#CHNmodelCAUExp<-CHNmodelCAUExp %>% mutate(model_pop="CHN",ID="CHNfilt")
CAUmodelCAUExp<-CAUmodelCAUExp %>% mutate(model_pop="CAU",ID="CAUfilt")
HISmodelCAUExp<-HISmodelCAUExp %>% mutate(model_pop="HIS",ID="HISfilt")

CAU_Exp_Correlations<-bind_rows(AFAmodelCAUExp,
                                ALLmodelCAUExp,
#                                CHNmodelCAUExp,
                                CAUmodelCAUExp,
                                HISmodelCAUExp)

### HIS
AFAmodelHISExp<-AFAmodelHISExp %>% mutate(model_pop="AFA",ID="AFAfilt")
ALLmodelHISExp<-ALLmodelHISExp %>% mutate(model_pop="ALL",ID="ALLfilt")
#CHNmodelHISExp<-CHNmodelHISExp %>% mutate(model_pop="CHN",ID="CHNfilt")
CAUmodelHISExp<-CAUmodelHISExp %>% mutate(model_pop="CAU",ID="CAUfilt")
HISmodelHISExp<-HISmodelHISExp %>% mutate(model_pop="HIS",ID="HISfilt")

HIS_Exp_Correlations<-bind_rows(AFAmodelHISExp,
                                ALLmodelHISExp,
#                                CHNmodelHISExp,
                                CAUmodelHISExp,
                                HISmodelHISExp)

## Plot the distributions of concat models

xlimit<-c(-1.2,1.2)
ylimit<-c(0,2500)
bin_width=0.05
### AFA
ggplot(data=AFA_Exp_Correlations, aes(x=estimate,fill=ID)) + geom_boxplot() + scale_x_continuous(limits=xlimit) + 
  scale_y_continuous(limits=ylimit) + facet_wrap(~ID) + ggtitle("" %&% tiss %&% " AFA Estimated Spearman Correlation Coefficient Frequecies") + 
  xlab("Estimated Spearman Correlation Coefficient") + ylab("Frequency")
ggsave("" %&% tiss %&% "_AFA_estimated_spearman_correlations.png")

### ALL
ggplot(data=ALL_Exp_Correlations, aes(x=estimate,fill=ID)) + geom_boxplot() + scale_x_continuous(limits=xlimit) + 
  scale_y_continuous(limits=ylimit) + facet_wrap(~ID) + ggtitle("" %&% tiss %&% " ALL Estimated Spearman Correlation Coefficient Frequecies") + 
  xlab("Estimated Spearman Correlation Coefficient") + ylab("Frequency")
ggsave("" %&% tiss %&% "_ALL_estimated_spearman_correlations.png")

### CHN
#ggplot(data=CHN_Exp_Correlations, aes(x=estimate,fill=ID)) + geom_boxplot() + scale_x_continuous(limits=xlimit) + 
#  scale_y_continuous(limits=ylimit) + facet_wrap(~ID) + ggtitle("" %&% tiss %&% " CHN Estimated Spearman Correlation Coefficient Frequecies") + 
#  xlab("Estimated Spearman Correlation Coefficient") + ylab("Frequency")
#ggsave("" %&% tiss %&% "_CHN_estimated_spearman_correlations.png")

### CAU
ggplot(data=CAU_Exp_Correlations, aes(x=estimate,fill=ID)) + geom_boxplot() + scale_x_continuous(limits=xlimit) + 
  scale_y_continuous(limits=ylimit) + facet_wrap(~ID) + ggtitle("" %&% tiss %&% " CAU Estimated Spearman Correlation Coefficient Frequecies") + 
  xlab("Estimated Spearman Correlation Coefficient") + ylab("Frequency")
ggsave("" %&% tiss %&% "_CAU_estimated_spearman_correlations.png")

### HIS
ggplot(data=HIS_Exp_Correlations, aes(x=estimate,fill=ID)) + geom_boxplot() + scale_x_continuous(limits=xlimit) + 
  scale_y_continuous(limits=ylimit) + facet_wrap(~ID) + ggtitle("" %&% tiss %&% " HIS Estimated Spearman Correlation Coefficient Frequecies") + 
  xlab("Estimated Spearman Correlation Coefficient") + ylab("Frequency")
ggsave("" %&% tiss %&% "_HIS_estimated_spearman_correlations.png")
