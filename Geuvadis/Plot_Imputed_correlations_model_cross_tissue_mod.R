cat("Compare Imputed Spearman Correlation")

library(data.table)
library(dplyr)
library(ggplot2)
#library(viridis)
library(tidyr)

"%&%" = function(a,b) paste(a,b,sep="")

## Global Variables
tiss="Tcell"
dir="/home/chris/geuvadis/baseline_elasticnet/imputed_results/PBMC/spearman/gene_model_overlap/overlap_"

## Read in spearman correlations

#Correlations calculated with `genetic_pi1.R` script

### PBMC AFA models
PBMC_AFA_model_geuvadis_ALL_Exp<-fread("" %&% dir %&% "PBMC_AFA_R2_-1_model_geuvadis_ALL_expression_spearman_correlation.txt",header=T)
PBMC_AFA_model_geuvadis_CEU_Exp<-fread("" %&% dir %&% "PBMC_AFA_R2_-1_model_geuvadis_CEU_expression_spearman_correlation.txt",header=T)
PBMC_AFA_model_geuvadis_FIN_Exp<-fread("" %&% dir %&% "PBMC_AFA_R2_-1_model_geuvadis_FIN_expression_spearman_correlation.txt",header=T)
PBMC_AFA_model_geuvadis_GBR_Exp<-fread("" %&% dir %&% "PBMC_AFA_R2_-1_model_geuvadis_GBR_expression_spearman_correlation.txt",header=T)
PBMC_AFA_model_geuvadis_TSI_Exp<-fread("" %&% dir %&% "PBMC_AFA_R2_-1_model_geuvadis_TSI_expression_spearman_correlation.txt",header=T)
PBMC_AFA_model_geuvadis_YRI_Exp<-fread("" %&% dir %&% "PBMC_AFA_R2_-1_model_geuvadis_YRI_expression_spearman_correlation.txt",header=T)

PBMC_ALL_model_geuvadis_ALL_Exp<-fread("" %&% dir %&% "PBMC_ALL_R2_-1_model_geuvadis_ALL_expression_spearman_correlation.txt",header=T)
PBMC_ALL_model_geuvadis_CEU_Exp<-fread("" %&% dir %&% "PBMC_ALL_R2_-1_model_geuvadis_CEU_expression_spearman_correlation.txt",header=T)
PBMC_ALL_model_geuvadis_FIN_Exp<-fread("" %&% dir %&% "PBMC_ALL_R2_-1_model_geuvadis_FIN_expression_spearman_correlation.txt",header=T)
PBMC_ALL_model_geuvadis_GBR_Exp<-fread("" %&% dir %&% "PBMC_ALL_R2_-1_model_geuvadis_GBR_expression_spearman_correlation.txt",header=T)
PBMC_ALL_model_geuvadis_TSI_Exp<-fread("" %&% dir %&% "PBMC_ALL_R2_-1_model_geuvadis_TSI_expression_spearman_correlation.txt",header=T)
PBMC_ALL_model_geuvadis_YRI_Exp<-fread("" %&% dir %&% "PBMC_ALL_R2_-1_model_geuvadis_YRI_expression_spearman_correlation.txt",header=T)

PBMC_CAU_model_geuvadis_ALL_Exp<-fread("" %&% dir %&% "PBMC_CAU_R2_-1_model_geuvadis_ALL_expression_spearman_correlation.txt",header=T)
PBMC_CAU_model_geuvadis_CEU_Exp<-fread("" %&% dir %&% "PBMC_CAU_R2_-1_model_geuvadis_CEU_expression_spearman_correlation.txt",header=T)
PBMC_CAU_model_geuvadis_FIN_Exp<-fread("" %&% dir %&% "PBMC_CAU_R2_-1_model_geuvadis_FIN_expression_spearman_correlation.txt",header=T)
PBMC_CAU_model_geuvadis_GBR_Exp<-fread("" %&% dir %&% "PBMC_CAU_R2_-1_model_geuvadis_GBR_expression_spearman_correlation.txt",header=T)
PBMC_CAU_model_geuvadis_TSI_Exp<-fread("" %&% dir %&% "PBMC_CAU_R2_-1_model_geuvadis_TSI_expression_spearman_correlation.txt",header=T)
PBMC_CAU_model_geuvadis_YRI_Exp<-fread("" %&% dir %&% "PBMC_CAU_R2_-1_model_geuvadis_YRI_expression_spearman_correlation.txt",header=T)

PBMC_CHN_model_geuvadis_ALL_Exp<-fread("" %&% dir %&% "PBMC_CHN_R2_-1_model_geuvadis_ALL_expression_spearman_correlation.txt",header=T)
PBMC_CHN_model_geuvadis_CEU_Exp<-fread("" %&% dir %&% "PBMC_CHN_R2_-1_model_geuvadis_CEU_expression_spearman_correlation.txt",header=T)
PBMC_CHN_model_geuvadis_FIN_Exp<-fread("" %&% dir %&% "PBMC_CHN_R2_-1_model_geuvadis_FIN_expression_spearman_correlation.txt",header=T)
PBMC_CHN_model_geuvadis_GBR_Exp<-fread("" %&% dir %&% "PBMC_CHN_R2_-1_model_geuvadis_GBR_expression_spearman_correlation.txt",header=T)
PBMC_CHN_model_geuvadis_TSI_Exp<-fread("" %&% dir %&% "PBMC_CHN_R2_-1_model_geuvadis_TSI_expression_spearman_correlation.txt",header=T)
PBMC_CHN_model_geuvadis_YRI_Exp<-fread("" %&% dir %&% "PBMC_CHN_R2_-1_model_geuvadis_YRI_expression_spearman_correlation.txt",header=T)

PBMC_HIS_model_geuvadis_ALL_Exp<-fread("" %&% dir %&% "PBMC_HIS_R2_-1_model_geuvadis_ALL_expression_spearman_correlation.txt",header=T)
PBMC_HIS_model_geuvadis_CEU_Exp<-fread("" %&% dir %&% "PBMC_HIS_R2_-1_model_geuvadis_CEU_expression_spearman_correlation.txt",header=T)
PBMC_HIS_model_geuvadis_FIN_Exp<-fread("" %&% dir %&% "PBMC_HIS_R2_-1_model_geuvadis_FIN_expression_spearman_correlation.txt",header=T)
PBMC_HIS_model_geuvadis_GBR_Exp<-fread("" %&% dir %&% "PBMC_HIS_R2_-1_model_geuvadis_GBR_expression_spearman_correlation.txt",header=T)
PBMC_HIS_model_geuvadis_TSI_Exp<-fread("" %&% dir %&% "PBMC_HIS_R2_-1_model_geuvadis_TSI_expression_spearman_correlation.txt",header=T)
PBMC_HIS_model_geuvadis_YRI_Exp<-fread("" %&% dir %&% "PBMC_HIS_R2_-1_model_geuvadis_YRI_expression_spearman_correlation.txt",header=T)

PBMC_AFA_model_geuvadis_ALL_Exp<-PBMC_AFA_model_geuvadis_ALL_Exp %>% mutate(model_pop="PBMC AFA",ID="ALL")
PBMC_AFA_model_geuvadis_CEU_Exp<-PBMC_AFA_model_geuvadis_CEU_Exp %>% mutate(model_pop="PBMC AFA",ID="CEU")
PBMC_AFA_model_geuvadis_FIN_Exp<-PBMC_AFA_model_geuvadis_FIN_Exp %>% mutate(model_pop="PBMC AFA",ID="FIN")
PBMC_AFA_model_geuvadis_GBR_Exp<-PBMC_AFA_model_geuvadis_GBR_Exp %>% mutate(model_pop="PBMC AFA",ID="GBR")
PBMC_AFA_model_geuvadis_TSI_Exp<-PBMC_AFA_model_geuvadis_TSI_Exp %>% mutate(model_pop="PBMC AFA",ID="TSI")
PBMC_AFA_model_geuvadis_YRI_Exp<-PBMC_AFA_model_geuvadis_YRI_Exp %>% mutate(model_pop="PBMC AFA",ID="YRI")

PBMC_ALL_model_geuvadis_ALL_Exp<-PBMC_ALL_model_geuvadis_ALL_Exp %>% mutate(model_pop="PBMC ALL",ID="ALL")
PBMC_ALL_model_geuvadis_CEU_Exp<-PBMC_ALL_model_geuvadis_CEU_Exp %>% mutate(model_pop="PBMC ALL",ID="CEU")
PBMC_ALL_model_geuvadis_FIN_Exp<-PBMC_ALL_model_geuvadis_FIN_Exp %>% mutate(model_pop="PBMC ALL",ID="FIN")
PBMC_ALL_model_geuvadis_GBR_Exp<-PBMC_ALL_model_geuvadis_GBR_Exp %>% mutate(model_pop="PBMC ALL",ID="GBR")
PBMC_ALL_model_geuvadis_TSI_Exp<-PBMC_ALL_model_geuvadis_TSI_Exp %>% mutate(model_pop="PBMC ALL",ID="TSI")
PBMC_ALL_model_geuvadis_YRI_Exp<-PBMC_ALL_model_geuvadis_YRI_Exp %>% mutate(model_pop="PBMC ALL",ID="YRI")

PBMC_CAU_model_geuvadis_ALL_Exp<-PBMC_CAU_model_geuvadis_ALL_Exp %>% mutate(model_pop="PBMC CAU",ID="ALL")
PBMC_CAU_model_geuvadis_CEU_Exp<-PBMC_CAU_model_geuvadis_CEU_Exp %>% mutate(model_pop="PBMC CAU",ID="CEU")
PBMC_CAU_model_geuvadis_FIN_Exp<-PBMC_CAU_model_geuvadis_FIN_Exp %>% mutate(model_pop="PBMC CAU",ID="FIN")
PBMC_CAU_model_geuvadis_GBR_Exp<-PBMC_CAU_model_geuvadis_GBR_Exp %>% mutate(model_pop="PBMC CAU",ID="GBR")
PBMC_CAU_model_geuvadis_TSI_Exp<-PBMC_CAU_model_geuvadis_TSI_Exp %>% mutate(model_pop="PBMC CAU",ID="TSI")
PBMC_CAU_model_geuvadis_YRI_Exp<-PBMC_CAU_model_geuvadis_YRI_Exp %>% mutate(model_pop="PBMC CAU",ID="YRI")

PBMC_CHN_model_geuvadis_ALL_Exp<-PBMC_CHN_model_geuvadis_ALL_Exp %>% mutate(model_pop="PBMC CHN",ID="ALL")
PBMC_CHN_model_geuvadis_CEU_Exp<-PBMC_CHN_model_geuvadis_CEU_Exp %>% mutate(model_pop="PBMC CHN",ID="CEU")
PBMC_CHN_model_geuvadis_FIN_Exp<-PBMC_CHN_model_geuvadis_FIN_Exp %>% mutate(model_pop="PBMC CHN",ID="FIN")
PBMC_CHN_model_geuvadis_GBR_Exp<-PBMC_CHN_model_geuvadis_GBR_Exp %>% mutate(model_pop="PBMC CHN",ID="GBR")
PBMC_CHN_model_geuvadis_TSI_Exp<-PBMC_CHN_model_geuvadis_TSI_Exp %>% mutate(model_pop="PBMC CHN",ID="TSI")
PBMC_CHN_model_geuvadis_YRI_Exp<-PBMC_CHN_model_geuvadis_YRI_Exp %>% mutate(model_pop="PBMC CHN",ID="YRI")

PBMC_HIS_model_geuvadis_ALL_Exp<-PBMC_HIS_model_geuvadis_ALL_Exp %>% mutate(model_pop="PBMC HIS",ID="ALL")
PBMC_HIS_model_geuvadis_CEU_Exp<-PBMC_HIS_model_geuvadis_CEU_Exp %>% mutate(model_pop="PBMC HIS",ID="CEU")
PBMC_HIS_model_geuvadis_FIN_Exp<-PBMC_HIS_model_geuvadis_FIN_Exp %>% mutate(model_pop="PBMC HIS",ID="FIN")
PBMC_HIS_model_geuvadis_GBR_Exp<-PBMC_HIS_model_geuvadis_GBR_Exp %>% mutate(model_pop="PBMC HIS",ID="GBR")
PBMC_HIS_model_geuvadis_TSI_Exp<-PBMC_HIS_model_geuvadis_TSI_Exp %>% mutate(model_pop="PBMC HIS",ID="TSI")
PBMC_HIS_model_geuvadis_YRI_Exp<-PBMC_HIS_model_geuvadis_YRI_Exp %>% mutate(model_pop="PBMC HIS",ID="YRI")

AFA_Model_Exp_Correlations<-bind_rows(PBMC_AFA_model_geuvadis_ALL_Exp,
				      PBMC_AFA_model_geuvadis_CEU_Exp,
                                      PBMC_AFA_model_geuvadis_FIN_Exp,
                                      PBMC_AFA_model_geuvadis_GBR_Exp,
                                      PBMC_AFA_model_geuvadis_TSI_Exp,
                                      PBMC_AFA_model_geuvadis_YRI_Exp)
ALL_Model_Exp_Correlations<-bind_rows(PBMC_ALL_model_geuvadis_ALL_Exp,
                                      PBMC_ALL_model_geuvadis_CEU_Exp,
                                      PBMC_ALL_model_geuvadis_FIN_Exp,
                                      PBMC_ALL_model_geuvadis_GBR_Exp,
                                      PBMC_ALL_model_geuvadis_TSI_Exp,
                                      PBMC_ALL_model_geuvadis_YRI_Exp)
CAU_Model_Exp_Correlations<-bind_rows(PBMC_CAU_model_geuvadis_ALL_Exp,
				      PBMC_CAU_model_geuvadis_CEU_Exp,
                                      PBMC_CAU_model_geuvadis_FIN_Exp,
                                      PBMC_CAU_model_geuvadis_GBR_Exp,
                                      PBMC_CAU_model_geuvadis_TSI_Exp,
                                      PBMC_CAU_model_geuvadis_YRI_Exp)
CHN_Model_Exp_Correlations<-bind_rows(PBMC_CHN_model_geuvadis_ALL_Exp,
				      PBMC_CHN_model_geuvadis_CEU_Exp,
                                      PBMC_CHN_model_geuvadis_FIN_Exp,
                                      PBMC_CHN_model_geuvadis_GBR_Exp,
                                      PBMC_CHN_model_geuvadis_TSI_Exp,
                                      PBMC_CHN_model_geuvadis_YRI_Exp)
HIS_Model_Exp_Correlations<-bind_rows(PBMC_HIS_model_geuvadis_ALL_Exp,
				      PBMC_HIS_model_geuvadis_CEU_Exp,
                                      PBMC_HIS_model_geuvadis_FIN_Exp,
                                      PBMC_HIS_model_geuvadis_GBR_Exp,
                                      PBMC_HIS_model_geuvadis_TSI_Exp,
                                      PBMC_HIS_model_geuvadis_YRI_Exp)


## Plot the distributions of concat models

ylimit<-c(-1,1)

violin_boxplot_width=0.1
boxplot_width=0.9
boxplot_xlimit<-c(-1,1)
boxplot_ylimit<-c(0,100)
color<-c('#A1D2CE','#78CAD2','#62A8AC','#5490A7','#50858B','#49516F')

### AFA model
ggplot(data=AFA_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_violin() + geom_boxplot(width=violin_boxplot_width)  +
  scale_y_continuous(limits=ylimit) + ggtitle("TOPMed PBMC AFA Baseline Model in Geuvadis (Model Gene Overlap)") + 
  guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
  scale_fill_manual(values = color)
ggsave("model_gene_overlap_PBMC_AFA_model_geuvadis_estimated_spearman_violin.png",width=12,height=7)

# ggplot(data=AFA_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_boxplot(width=boxplot_width)  +
#   scale_y_continuous(limits=ylimit) + ggtitle("TOPMed PBMC AFA Baseline Model in Geuvadis") + 
#   guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
#   scale_fill_manual(values = color)
# ggsave("PBMC_AFA_model_geuvadis_estimated_spearman_boxplot.png",width=12,height=7)
# 
# 
# ggplot(data=AFA_Model_Exp_Correlations, aes(x=estimate,color=ID)) + geom_density()  +
#   scale_x_continuous(limits=boxplot_xlimit) + scale_y_discrete(limits=boxplot_ylimit) + ggtitle("TOPMed PBMC AFA Baseline Model in Geuvadis") + 
#   guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Estimated Spearman Correlation Coefficient") + ylab("Frequency") +
#   scale_color_manual(values = color)
# ggsave("PBMC_AFA_model_geuvadis_estimated_spearman_density.png",width=12,height=7)
# 
# write.table(summary(AFA_Model_Exp_Correlations),file = "PBMC_AFA_model_geuvadis_estimated_spearman_summary.txt", append = FALSE, sep = "\t", quote = FALSE, dec = ".",
#            row.names = FALSE, col.names = TRUE)


### ALL model
ggplot(data=ALL_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_violin() + geom_boxplot(width=violin_boxplot_width)  +
  scale_y_continuous(limits=ylimit) + ggtitle("TOPMed PBMC ALL Baseline Model in Geuvadis (Model Gene Overlap)") + 
  guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
  scale_fill_manual(values = color)
ggsave("model_gene_overlap_PBMC_ALL_model_geuvadis_estimated_spearman_violin.png",width=12,height=7)

# ggplot(data=ALL_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_boxplot(width=boxplot_width)  +
#   scale_y_continuous(limits=ylimit) + ggtitle("TOPMed PBMC ALL Baseline Model in Geuvadis") + 
#   guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
#   scale_fill_manual(values = color)
# ggsave("PBMC_ALL_model_geuvadis_estimated_spearman_boxplot.png",width=12,height=7)
# 
# ggplot(data=ALL_Model_Exp_Correlations, aes(x=estimate,color=ID)) + geom_density()  +
#   scale_x_continuous(limits=boxplot_xlimit) + scale_y_discrete(limits=boxplot_ylimit) + ggtitle("TOPMed PBMC ALL Baseline Model in Geuvadis") + 
#   guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Estimated Spearman Correlation Coefficient") + ylab("Frequency") +
#   scale_color_manual(values = color)
# ggsave("PBMC_ALL_model_geuvadis_estimated_spearman_density.png",width=12,height=7)
# 
# write.table(summary(ALL_Model_Exp_Correlations),file = "PBMC_ALL_model_geuvadis_estimated_spearman_summary.txt", append = FALSE, sep = "\t", quote = FALSE, dec = ".",
#             row.names = FALSE, col.names = TRUE)

### CAU model
ggplot(data=CAU_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_violin() + geom_boxplot(width=violin_boxplot_width)  +
  scale_y_continuous(limits=ylimit) + ggtitle("TOPMed PBMC CAU Baseline Model in Geuvadis (Model Gene Overlap)") + 
  guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
  scale_fill_manual(values = color)
ggsave("model_gene_overlap_PBMC_CAU_model_geuvadis_estimated_spearman_violin.png",width=12,height=7)

# ggplot(data=CAU_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_boxplot(width=boxplot_width)  +
#   scale_y_continuous(limits=ylimit) + ggtitle("TOPMed PBMC CAU Baseline Model in Geuvadis") + 
#   guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
#   scale_fill_manual(values = color)
# ggsave("PBMC_CAU_model_geuvadis_estimated_spearman_boxplot.png",width=12,height=7)
# 
# ggplot(data=CAU_Model_Exp_Correlations, aes(x=estimate,color=ID)) + geom_density()  +
#   scale_x_continuous(limits=boxplot_xlimit) + scale_y_discrete(limits=boxplot_ylimit) + ggtitle("TOPMed PBMC CAU Baseline Model in Geuvadis") + 
#   guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Estimated Spearman Correlation Coefficient") + ylab("Frequency") +
#   scale_color_manual(values = color)
# ggsave("PBMC_CAU_model_geuvadis_estimated_spearman_density.png",width=12,height=7)
# 
# write.table(summary(CAU_Model_Exp_Correlations),file = "PBMC_CAU_model_geuvadis_estimated_spearman_summary.txt", append = FALSE, sep = "\t", quote = FALSE, dec = ".",
#             row.names = FALSE, col.names = TRUE)

### CHN model
ggplot(data=CHN_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_violin() + geom_boxplot(width=violin_boxplot_width)  +
  scale_y_continuous(limits=ylimit) + ggtitle("TOPMed PBMC CHN Baseline Model in Geuvadis (Model Gene Overlap)") + 
  guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
  scale_fill_manual(values = color)
ggsave("model_gene_overlap_PBMC_CHN_model_geuvadis_estimated_spearman_violin.png",width=12,height=7)

# ggplot(data=CHN_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_boxplot(width=boxplot_width)  +
#   scale_y_continuous(limits=ylimit) + ggtitle("TOPMed PBMC CHN Baseline Model in Geuvadis") + 
#   guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
#   scale_fill_manual(values = color)
# ggsave("PBMC_CHN_model_geuvadis_estimated_spearman_boxplot.png",width=12,height=7)
# 
# ggplot(data=CHN_Model_Exp_Correlations, aes(x=estimate,color=ID)) + geom_density()  +
#   scale_x_continuous(limits=boxplot_xlimit) + scale_y_discrete(limits=boxplot_ylimit) + ggtitle("TOPMed PBMC CHN Baseline Model in Geuvadis") + 
#   guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Estimated Spearman Correlation Coefficient") + ylab("Frequency") +
#   scale_color_manual(values = color)
# ggsave("PBMC_CHN_model_geuvadis_estimated_spearman_denisty.png",width=12,height=7)
# 
# write.table(summary(CHN_Model_Exp_Correlations),file = "PBMC_CHN_model_geuvadis_estimated_spearman_summary.txt", append = FALSE, sep = "\t", quote = FALSE, dec = ".",
#             row.names = FALSE, col.names = TRUE)

### HIS model
ggplot(data=HIS_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_violin() + geom_boxplot(width=violin_boxplot_width)  +
  scale_y_continuous(limits=ylimit) + ggtitle("TOPMed PBMC HIS Baseline Model in Geuvadis (Model Gene Overlap)") + 
  guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
  scale_fill_manual(values = color)
ggsave("model_gene_overlap_PBMC_HIS_model_geuvadis_estimated_spearman_violin.png",width=12,height=7)

# ggplot(data=HIS_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_boxplot(width=boxplot_width)  +
#   scale_y_continuous(limits=ylimit) + ggtitle("TOPMed PBMC HIS Baseline Model in Geuvadis") + 
#   guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
#   scale_fill_manual(values = color)
# ggsave("PBMC_HIS_model_geuvadis_estimated_spearman_boxplot.png",width=12,height=7)
# 
# ggplot(data=HIS_Model_Exp_Correlations, aes(x=estimate,color=ID)) + geom_density()  +
#   scale_x_continuous(limits=boxplot_xlimit) + scale_y_discrete(limits=boxplot_ylimit) + ggtitle("TOPMed PBMC HIS Baseline Model in Geuvadis") + 
#   guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Estimated Spearman Correlation Coefficient") + ylab("Frequency") +
#   scale_color_manual(values = color)
# ggsave("PBMC_HIS_model_geuvadis_estimated_spearman_density.png",width=12,height=7)
# 
# write.table(summary(HIS_Model_Exp_Correlations),file = "PBMC_HIS_model_geuvadis_estimated_spearman_summary.txt", append = FALSE, sep = "\t", quote = FALSE, dec = ".",
#             row.names = FALSE, col.names = TRUE)
