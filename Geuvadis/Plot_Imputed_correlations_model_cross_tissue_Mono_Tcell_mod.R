cat("Compare Imputed Spearman Correlation")

library(data.table)
library(dplyr)
library(ggplot2)
#library(viridis)
library(tidyr)

"%&%" = function(a,b) paste(a,b,sep="")

## Global Variables
tiss="Tcell"
dir="/home/chris/geuvadis/baseline_elasticnet/imputed_results/"

## Read in spearman correlations

#Correlations calculated with `genetic_pi1.R` script

### Tcell AFA models
Tcell_AFA_model_geuvadis_ALL_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_AFA_R2_-1_model_geuvadis_ALL_expression_spearman_correlation.txt",header=T)
Tcell_AFA_model_geuvadis_CEU_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_AFA_R2_-1_model_geuvadis_CEU_expression_spearman_correlation.txt",header=T)
Tcell_AFA_model_geuvadis_FIN_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_AFA_R2_-1_model_geuvadis_FIN_expression_spearman_correlation.txt",header=T)
Tcell_AFA_model_geuvadis_GBR_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_AFA_R2_-1_model_geuvadis_GBR_expression_spearman_correlation.txt",header=T)
Tcell_AFA_model_geuvadis_TSI_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_AFA_R2_-1_model_geuvadis_TSI_expression_spearman_correlation.txt",header=T)
Tcell_AFA_model_geuvadis_YRI_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_AFA_R2_-1_model_geuvadis_YRI_expression_spearman_correlation.txt",header=T)

Tcell_ALL_model_geuvadis_ALL_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_ALL_R2_-1_model_geuvadis_ALL_expression_spearman_correlation.txt",header=T)
Tcell_ALL_model_geuvadis_CEU_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_ALL_R2_-1_model_geuvadis_CEU_expression_spearman_correlation.txt",header=T)
Tcell_ALL_model_geuvadis_FIN_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_ALL_R2_-1_model_geuvadis_FIN_expression_spearman_correlation.txt",header=T)
Tcell_ALL_model_geuvadis_GBR_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_ALL_R2_-1_model_geuvadis_GBR_expression_spearman_correlation.txt",header=T)
Tcell_ALL_model_geuvadis_TSI_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_ALL_R2_-1_model_geuvadis_TSI_expression_spearman_correlation.txt",header=T)
Tcell_ALL_model_geuvadis_YRI_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_ALL_R2_-1_model_geuvadis_YRI_expression_spearman_correlation.txt",header=T)

Tcell_CAU_model_geuvadis_ALL_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_CAU_R2_-1_model_geuvadis_ALL_expression_spearman_correlation.txt",header=T)
Tcell_CAU_model_geuvadis_CEU_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_CAU_R2_-1_model_geuvadis_CEU_expression_spearman_correlation.txt",header=T)
Tcell_CAU_model_geuvadis_FIN_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_CAU_R2_-1_model_geuvadis_FIN_expression_spearman_correlation.txt",header=T)
Tcell_CAU_model_geuvadis_GBR_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_CAU_R2_-1_model_geuvadis_GBR_expression_spearman_correlation.txt",header=T)
Tcell_CAU_model_geuvadis_TSI_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_CAU_R2_-1_model_geuvadis_TSI_expression_spearman_correlation.txt",header=T)
Tcell_CAU_model_geuvadis_YRI_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_CAU_R2_-1_model_geuvadis_YRI_expression_spearman_correlation.txt",header=T)

Tcell_HIS_model_geuvadis_ALL_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_HIS_R2_-1_model_geuvadis_ALL_expression_spearman_correlation.txt",header=T)
Tcell_HIS_model_geuvadis_CEU_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_HIS_R2_-1_model_geuvadis_CEU_expression_spearman_correlation.txt",header=T)
Tcell_HIS_model_geuvadis_FIN_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_HIS_R2_-1_model_geuvadis_FIN_expression_spearman_correlation.txt",header=T)
Tcell_HIS_model_geuvadis_GBR_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_HIS_R2_-1_model_geuvadis_GBR_expression_spearman_correlation.txt",header=T)
Tcell_HIS_model_geuvadis_TSI_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_HIS_R2_-1_model_geuvadis_TSI_expression_spearman_correlation.txt",header=T)
Tcell_HIS_model_geuvadis_YRI_Exp<-fread("" %&% dir %&% "Tcell/spearman/gene_model_overlap/overlap_Tcell_HIS_R2_-1_model_geuvadis_YRI_expression_spearman_correlation.txt",header=T)

Tcell_AFA_model_geuvadis_ALL_Exp<-Tcell_AFA_model_geuvadis_ALL_Exp %>% mutate(model_pop="Tcell AFA",ID="ALL")
Tcell_AFA_model_geuvadis_CEU_Exp<-Tcell_AFA_model_geuvadis_CEU_Exp %>% mutate(model_pop="Tcell AFA",ID="CEU")
Tcell_AFA_model_geuvadis_FIN_Exp<-Tcell_AFA_model_geuvadis_FIN_Exp %>% mutate(model_pop="Tcell AFA",ID="FIN")
Tcell_AFA_model_geuvadis_GBR_Exp<-Tcell_AFA_model_geuvadis_GBR_Exp %>% mutate(model_pop="Tcell AFA",ID="GBR")
Tcell_AFA_model_geuvadis_TSI_Exp<-Tcell_AFA_model_geuvadis_TSI_Exp %>% mutate(model_pop="Tcell AFA",ID="TSI")
Tcell_AFA_model_geuvadis_YRI_Exp<-Tcell_AFA_model_geuvadis_YRI_Exp %>% mutate(model_pop="Tcell AFA",ID="YRI")

Tcell_ALL_model_geuvadis_ALL_Exp<-Tcell_ALL_model_geuvadis_ALL_Exp %>% mutate(model_pop="Tcell ALL",ID="ALL")
Tcell_ALL_model_geuvadis_CEU_Exp<-Tcell_ALL_model_geuvadis_CEU_Exp %>% mutate(model_pop="Tcell ALL",ID="CEU")
Tcell_ALL_model_geuvadis_FIN_Exp<-Tcell_ALL_model_geuvadis_FIN_Exp %>% mutate(model_pop="Tcell ALL",ID="FIN")
Tcell_ALL_model_geuvadis_GBR_Exp<-Tcell_ALL_model_geuvadis_GBR_Exp %>% mutate(model_pop="Tcell ALL",ID="GBR")
Tcell_ALL_model_geuvadis_TSI_Exp<-Tcell_ALL_model_geuvadis_TSI_Exp %>% mutate(model_pop="Tcell ALL",ID="TSI")
Tcell_ALL_model_geuvadis_YRI_Exp<-Tcell_ALL_model_geuvadis_YRI_Exp %>% mutate(model_pop="Tcell ALL",ID="YRI")

Tcell_CAU_model_geuvadis_ALL_Exp<-Tcell_CAU_model_geuvadis_ALL_Exp %>% mutate(model_pop="Tcell CAU",ID="ALL")
Tcell_CAU_model_geuvadis_CEU_Exp<-Tcell_CAU_model_geuvadis_CEU_Exp %>% mutate(model_pop="Tcell CAU",ID="CEU")
Tcell_CAU_model_geuvadis_FIN_Exp<-Tcell_CAU_model_geuvadis_FIN_Exp %>% mutate(model_pop="Tcell CAU",ID="FIN")
Tcell_CAU_model_geuvadis_GBR_Exp<-Tcell_CAU_model_geuvadis_GBR_Exp %>% mutate(model_pop="Tcell CAU",ID="GBR")
Tcell_CAU_model_geuvadis_TSI_Exp<-Tcell_CAU_model_geuvadis_TSI_Exp %>% mutate(model_pop="Tcell CAU",ID="TSI")
Tcell_CAU_model_geuvadis_YRI_Exp<-Tcell_CAU_model_geuvadis_YRI_Exp %>% mutate(model_pop="Tcell CAU",ID="YRI")

Tcell_HIS_model_geuvadis_ALL_Exp<-Tcell_HIS_model_geuvadis_ALL_Exp %>% mutate(model_pop="Tcell HIS",ID="ALL")
Tcell_HIS_model_geuvadis_CEU_Exp<-Tcell_HIS_model_geuvadis_CEU_Exp %>% mutate(model_pop="Tcell HIS",ID="CEU")
Tcell_HIS_model_geuvadis_FIN_Exp<-Tcell_HIS_model_geuvadis_FIN_Exp %>% mutate(model_pop="Tcell HIS",ID="FIN")
Tcell_HIS_model_geuvadis_GBR_Exp<-Tcell_HIS_model_geuvadis_GBR_Exp %>% mutate(model_pop="Tcell HIS",ID="GBR")
Tcell_HIS_model_geuvadis_TSI_Exp<-Tcell_HIS_model_geuvadis_TSI_Exp %>% mutate(model_pop="Tcell HIS",ID="TSI")
Tcell_HIS_model_geuvadis_YRI_Exp<-Tcell_HIS_model_geuvadis_YRI_Exp %>% mutate(model_pop="Tcell HIS",ID="YRI")

AFA_Model_Exp_Correlations<-bind_rows(Tcell_AFA_model_geuvadis_ALL_Exp,
                                      Tcell_AFA_model_geuvadis_CEU_Exp,
                                      Tcell_AFA_model_geuvadis_FIN_Exp,
                                      Tcell_AFA_model_geuvadis_GBR_Exp,
                                      Tcell_AFA_model_geuvadis_TSI_Exp,
                                      Tcell_AFA_model_geuvadis_YRI_Exp)
ALL_Model_Exp_Correlations<-bind_rows(Tcell_ALL_model_geuvadis_ALL_Exp,
                                      Tcell_ALL_model_geuvadis_CEU_Exp,
                                      Tcell_ALL_model_geuvadis_FIN_Exp,
                                      Tcell_ALL_model_geuvadis_GBR_Exp,
                                      Tcell_ALL_model_geuvadis_TSI_Exp,
                                      Tcell_ALL_model_geuvadis_YRI_Exp)
CAU_Model_Exp_Correlations<-bind_rows(Tcell_CAU_model_geuvadis_ALL_Exp,
                                      Tcell_CAU_model_geuvadis_CEU_Exp,
                                      Tcell_CAU_model_geuvadis_FIN_Exp,
                                      Tcell_CAU_model_geuvadis_GBR_Exp,
                                      Tcell_CAU_model_geuvadis_TSI_Exp,
                                      Tcell_CAU_model_geuvadis_YRI_Exp)
HIS_Model_Exp_Correlations<-bind_rows(Tcell_HIS_model_geuvadis_ALL_Exp,
                                      Tcell_HIS_model_geuvadis_CEU_Exp,
                                      Tcell_HIS_model_geuvadis_FIN_Exp,
                                      Tcell_HIS_model_geuvadis_GBR_Exp,
                                      Tcell_HIS_model_geuvadis_TSI_Exp,
                                      Tcell_HIS_model_geuvadis_YRI_Exp)


## Plot the distributions of concat models

ylimit<-c(-1,1)

violin_boxplot_width=0.1
boxplot_width=0.9
boxplot_xlimit<-c(-1,1)
boxplot_ylimit<-c(0,100)
color<-c('#A1D2CE','#78CAD2','#62A8AC','#5490A7','#50858B','#49516F')

### AFA model
ggplot(data=AFA_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_violin() + geom_boxplot(width=violin_boxplot_width)  +
  scale_y_continuous(limits=ylimit) + ggtitle("TOPMed Tcell AFA Baseline Model in Geuvadis (Model Gene Overlap)") + 
  guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
  scale_fill_manual(values = color)
ggsave("model_gene_overlap_Tcell_AFA_model_geuvadis_estimated_spearman_violin.png",width=12,height=7)

ggplot(data=AFA_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_boxplot(width=boxplot_width)  +
  scale_y_continuous(limits=ylimit) + ggtitle("TOPMed Tcell AFA Baseline Model in Geuvadis (Model Gene Overlap)") + 
  guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
  scale_fill_manual(values = color)
ggsave("model_gene_overlap_Tcell_AFA_model_geuvadis_estimated_spearman_boxplot.png",width=12,height=7)

ggplot(data=AFA_Model_Exp_Correlations, aes(x=estimate,color=ID)) + geom_density()  +
  scale_x_continuous(limits=boxplot_xlimit) + scale_y_discrete(limits=boxplot_ylimit) + ggtitle("TOPMed Tcell AFA Baseline Model in Geuvadis (Model Gene Overlap)") + 
  guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Estimated Spearman Correlation Coefficient") + ylab("Frequency") +
  scale_color_manual(values = color)
ggsave("model_gene_overlap_Tcell_AFA_model_geuvadis_estimated_spearman_density.png",width=12,height=7)

write.table(summary(AFA_Model_Exp_Correlations),file = "model_gene_overlap_Tcell_AFA_model_geuvadis_estimated_spearman_summary.txt", append = FALSE, sep = "\t", quote = FALSE, dec = ".",
            row.names = FALSE, col.names = TRUE)

### ALL model
ggplot(data=ALL_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_violin() + geom_boxplot(width=violin_boxplot_width)  +
  scale_y_continuous(limits=ylimit) + ggtitle("TOPMed Tcell ALL Baseline Model in Geuvadis (Model Gene Overlap)") + 
  guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
  scale_fill_manual(values = color)
ggsave("model_gene_overlap_Tcell_ALL_model_geuvadis_estimated_spearman_violin.png",width=12,height=7)

ggplot(data=ALL_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_boxplot(width=boxplot_width)  +
  scale_y_continuous(limits=ylimit) + ggtitle("TOPMed Tcell ALL Baseline Model in Geuvadis (Model Gene Overlap)") + 
  guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
  scale_fill_manual(values = color)
ggsave("model_gene_overlap_Tcell_ALL_model_geuvadis_estimated_spearman_boxplot.png",width=12,height=7)

ggplot(data=ALL_Model_Exp_Correlations, aes(x=estimate,color=ID)) + geom_density()  +
  scale_x_continuous(limits=boxplot_xlimit) + scale_y_discrete(limits=boxplot_ylimit) + ggtitle("TOPMed Tcell ALL Baseline Model in Geuvadis (Model Gene Overlap)") + 
  guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Estimated Spearman Correlation Coefficient") + ylab("Frequency") +
  scale_color_manual(values = color)
ggsave("model_gene_overlap_Tcell_ALL_model_geuvadis_estimated_spearman_density.png",width=12,height=7)

write.table(summary(ALL_Model_Exp_Correlations),file = "model_gene_overlap_Tcell_ALL_model_geuvadis_estimated_spearman_summary.txt", append = FALSE, sep = "\t", quote = FALSE, dec = ".",
            row.names = FALSE, col.names = TRUE)

### CAU model
ggplot(data=CAU_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_violin() + geom_boxplot(width=violin_boxplot_width)  +
  scale_y_continuous(limits=ylimit) + ggtitle("TOPMed Tcell CAU Baseline Model in Geuvadis (Model Gene Overlap)") + 
  guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
  scale_fill_manual(values = color)
ggsave("model_gene_overlap_Tcell_CAU_model_geuvadis_estimated_spearman_violin.png",width=12,height=7)

ggplot(data=CAU_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_boxplot(width=boxplot_width)  +
  scale_y_continuous(limits=ylimit) + ggtitle("TOPMed Tcell CAU Baseline Model in Geuvadis (Model Gene Overlap)") + 
  guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
  scale_fill_manual(values = color)
ggsave("model_gene_overlap_Tcell_CAU_model_geuvadis_estimated_spearman_boxplot.png",width=12,height=7)

ggplot(data=CAU_Model_Exp_Correlations, aes(x=estimate,color=ID)) + geom_density()  +
  scale_x_continuous(limits=boxplot_xlimit) + scale_y_discrete(limits=boxplot_ylimit) + ggtitle("TOPMed Tcell CAU Baseline Model in Geuvadis (Model Gene Overlap)") + 
  guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Estimated Spearman Correlation Coefficient") + ylab("Frequency") +
  scale_color_manual(values = color)
ggsave("model_gene_overlap_Tcell_CAU_model_geuvadis_estimated_spearman_density.png",width=12,height=7)

write.table(summary(CAU_Model_Exp_Correlations),file = "model_gene_overlap_Tcell_CAU_model_geuvadis_estimated_spearman_summary.txt", append = FALSE, sep = "\t", quote = FALSE, dec = ".",
            row.names = FALSE, col.names = TRUE)

### HIS model
ggplot(data=HIS_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_violin() + geom_boxplot(width=violin_boxplot_width)  +
  scale_y_continuous(limits=ylimit) + ggtitle("TOPMed Tcell HIS Baseline Model in Geuvadis (Model Gene Overlap)") + 
  guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
  scale_fill_manual(values = color)
ggsave("model_gene_overlap_Tcell_HIS_model_geuvadis_estimated_spearman_violin.png",width=12,height=7)

ggplot(data=HIS_Model_Exp_Correlations, aes(x=ID,y=estimate,fill=ID)) + geom_boxplot(width=boxplot_width)  +
  scale_y_continuous(limits=ylimit) + ggtitle("TOPMed Tcell HIS Baseline Model in Geuvadis (Model Gene Overlap)") + 
  guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Geuvadis Populations") + ylab("Estimated Spearman Correlation Coefficient") +
  scale_fill_manual(values = color)
ggsave("model_gene_overlap_Tcell_HIS_model_geuvadis_estimated_spearman_boxplot.png",width=12,height=7)

ggplot(data=HIS_Model_Exp_Correlations, aes(x=estimate,color=ID)) + geom_density()  +
  scale_x_continuous(limits=boxplot_xlimit) + scale_y_discrete(limits=boxplot_ylimit) + ggtitle("TOPMed Tcell HIS Baseline Model in Geuvadis (Model Gene Overlap)") + 
  guides(fill = guide_legend(title="Geuvadis Populations")) + xlab("Estimated Spearman Correlation Coefficient") + ylab("Frequency") +
  scale_color_manual(values = color)
ggsave("model_gene_overlap_Tcell_HIS_model_geuvadis_estimated_spearman_density.png",width=12,height=7)

write.table(summary(HIS_Model_Exp_Correlations),file = "model_gene_overlap_Tcell_HIS_model_geuvadis_estimated_spearman_summary.txt", append = FALSE, sep = "\t", quote = FALSE, dec = ".",
            row.names = FALSE, col.names = TRUE)
