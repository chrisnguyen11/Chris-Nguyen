setwd("/home/chris/topmed/scripts")
source("01a_EN_baseline.R")
"%&%" <- function(a,b) paste(a,b, sep='')

#argv <- commandArgs(trailingOnly = TRUE)
#pop <- argv[1]
#chrom <- argv[2]
#pops: AFA,CHN,HIS,CAU,ALL

pop<-"ALL"
chrom<-"22"

snp_annot_file <- "/home/chris/topmed/expression/" %&% pop %&% "/Tcell/pred_db/pred_db_anno_hg38.chr" %&% chrom %&% ".maf0.01.R20.8.dosage.txt"
expression_annot_file <- "/home/chris/gencode/gene_annotation_hg38_protein_coding.txt"
genotype_file <- "/home/chris/topmed/expression/" %&% pop %&% "/Tcell/pred_db/pred_db_geno_hg38.chr" %&% chrom %&% ".maf0.01.R20.8.dosage.txt"
expression_file <- "/home/chris/topmed/expression/" %&% pop %&% "/Tcell/expression/Tcell_ALL_age_sex_exam_10PCS_PF_adj_rinv_TOPMED_expression10_peer_factor_adjusted_sorted.txt"
covariates_file <- ""
prefix <- "baseline_Tcell_ALL"

main(snp_annot_file, expression_annot_file, genotype_file, expression_file, covariates_file, as.numeric(chrom), prefix, null_testing=FALSE)
