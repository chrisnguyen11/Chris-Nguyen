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

### PBMC AFA models
PBMC_AFA_model_PBMC_AFA_Exp<-fread("" %&% dir %&% "PBMC/PBMC_AFA_R2_-1_model_AFA_expression_spearman_correlation.txt",header=T)
PBMC_AFA_model_PBMC_ALL_Exp<-fread("" %&% dir %&% "PBMC/PBMC_AFA_R2_-1_model_ALL_expression_spearman_correlation.txt",header=T)
PBMC_AFA_model_PBMC_CAU_Exp<-fread("" %&% dir %&% "PBMC/PBMC_AFA_R2_-1_model_CAU_expression_spearman_correlation.txt",header=T)
PBMC_AFA_model_PBMC_CHN_Exp<-fread("" %&% dir %&% "PBMC/PBMC_AFA_R2_-1_model_CHN_expression_spearman_correlation.txt",header=T)
PBMC_AFA_model_PBMC_HIS_Exp<-fread("" %&% dir %&% "PBMC/PBMC_AFA_R2_-1_model_HIS_expression_spearman_correlation.txt",header=T)

PBMC_AFA_model_Mono_AFA_Exp<-fread("" %&% dir %&% "/cross_tissue/PBMC_AFA_R2_-1_model_Mono_AFA_expression_spearman_correlation.txt",header=T)
PBMC_AFA_model_Mono_ALL_Exp<-fread("" %&% dir %&% "/cross_tissue/PBMC_AFA_R2_-1_model_Mono_ALL_expression_spearman_correlation.txt",header=T)
PBMC_AFA_model_Mono_CAU_Exp<-fread("" %&% dir %&% "/cross_tissue/PBMC_AFA_R2_-1_model_Mono_CAU_expression_spearman_correlation.txt",header=T)
PBMC_AFA_model_Mono_HIS_Exp<-fread("" %&% dir %&% "/cross_tissue/PBMC_AFA_R2_-1_model_Mono_HIS_expression_spearman_correlation.txt",header=T)

PBMC_AFA_model_Tcell_AFA_Exp<-fread("" %&% dir %&% "/cross_tissue/PBMC_AFA_R2_-1_model_Tcell_AFA_expression_spearman_correlation.txt",header=T)
PBMC_AFA_model_Tcell_ALL_Exp<-fread("" %&% dir %&% "/cross_tissue/PBMC_AFA_R2_-1_model_Tcell_ALL_expression_spearman_correlation.txt",header=T)
PBMC_AFA_model_Tcell_CAU_Exp<-fread("" %&% dir %&% "/cross_tissue/PBMC_AFA_R2_-1_model_Tcell_CAU_expression_spearman_correlation.txt",header=T)
PBMC_AFA_model_Tcell_HIS_Exp<-fread("" %&% dir %&% "/cross_tissue/PBMC_AFA_R2_-1_model_Tcell_HIS_expression_spearman_correlation.txt",header=T)

PBMC_AFA_model_PBMC_AFA_Exp<-PBMC_AFA_model_PBMC_AFA_Exp %>% mutate(model_pop="PBMC AFA",ID="PBMC AFA")
PBMC_AFA_model_PBMC_ALL_Exp<-PBMC_AFA_model_PBMC_ALL_Exp %>% mutate(model_pop="PBMC AFA",ID="PBMC ALL")
PBMC_AFA_model_PBMC_CAU_Exp<-PBMC_AFA_model_PBMC_CAU_Exp %>% mutate(model_pop="PBMC AFA",ID="PBMC CAU")
PBMC_AFA_model_PBMC_CHN_Exp<-PBMC_AFA_model_PBMC_CHN_Exp %>% mutate(model_pop="PBMC AFA",ID="PBMC CHN")
PBMC_AFA_model_PBMC_HIS_Exp<-PBMC_AFA_model_PBMC_HIS_Exp %>% mutate(model_pop="PBMC AFA",ID="PBMC HIS")

PBMC_AFA_model_Mono_AFA_Exp<-PBMC_AFA_model_Mono_AFA_Exp %>% mutate(model_pop="PBMC AFA",ID="Mono AFA")
PBMC_AFA_model_Mono_ALL_Exp<-PBMC_AFA_model_Mono_ALL_Exp %>% mutate(model_pop="PBMC AFA",ID="Mono ALL")
PBMC_AFA_model_Mono_CAU_Exp<-PBMC_AFA_model_Mono_CAU_Exp %>% mutate(model_pop="PBMC AFA",ID="Mono CAU")
PBMC_AFA_model_Mono_HIS_Exp<-PBMC_AFA_model_Mono_HIS_Exp %>% mutate(model_pop="PBMC AFA",ID="Mono HIS")

PBMC_AFA_model_Tcell_AFA_Exp<-PBMC_AFA_model_Tcell_AFA_Exp %>% mutate(model_pop="PBMC AFA",ID="Tcell AFA")
PBMC_AFA_model_Tcell_ALL_Exp<-PBMC_AFA_model_Tcell_ALL_Exp %>% mutate(model_pop="PBMC AFA",ID="Tcell ALL")
PBMC_AFA_model_Tcell_CAU_Exp<-PBMC_AFA_model_Tcell_CAU_Exp %>% mutate(model_pop="PBMC AFA",ID="Tcell CAU")
PBMC_AFA_model_Tcell_HIS_Exp<-PBMC_AFA_model_Tcell_HIS_Exp %>% mutate(model_pop="PBMC AFA",ID="Tcell HIS")

AFA_Model_Exp_Correlations<-bind_rows(PBMC_AFA_model_PBMC_AFA_Exp,
#                                      PBMC_AFA_model_PBMC_ALL_Exp,
                                      PBMC_AFA_model_PBMC_CAU_Exp,
                                      PBMC_AFA_model_PBMC_CHN_Exp,
                                      PBMC_AFA_model_PBMC_HIS_Exp,
                                      PBMC_AFA_model_Mono_AFA_Exp,
#                                      PBMC_AFA_model_Mono_ALL_Exp,
                                      PBMC_AFA_model_Mono_CAU_Exp,
                                      PBMC_AFA_model_Mono_HIS_Exp,
                                      PBMC_AFA_model_Tcell_AFA_Exp,
#                                      PBMC_AFA_model_Tcell_ALL_Exp,
                                      PBMC_AFA_model_Tcell_ALL_Exp,
                                      PBMC_AFA_model_Tcell_HIS_Exp)

## Plot the distributions of concat models

ylimit<-c(-1,1)
boxplot_width=0.1

### AFA model
ggplot(data=AFA_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_violin() + geom_boxplot(width=boxplot_width)  +
  scale_y_continuous(limits=ylimit) + ggtitle("PBMC AFA Model - Estimated Spearman Correlation Coefficient - Cross Tissue") + 
  guides(fill = guide_legend(title="Populations")) + xlab("Populations") + ylab("Estimated Spearman Correlation Coefficient") +
  scale_fill_manual(values = c('salmon1','salmon2','salmon3','plum','plum1','plum2','plum3','steelblue1','steelblue2','steelblue3'))
ggsave("cross_tissue_PBMC_AFA_model_estimated_spearman_correlations.png",width=12,height=7)
