# Building TOPMed Models
## Input Files 
SNP annotation files - 
```
chr pos varID refAllele effectAllele rsid
chr1 632373 chr1_632373_G_A_b38 G A chr1:632373
chr1 775010 chr1_775010_A_C_b38 A C chr1:775010
```
expression annotation files - 
```
chr gene_id gene_name start end gene_type
1 ENSG00000186092 OR4F5 65419 71585 protein_coding
1 ENSG00000284733 OR4F29 450703 451697 protein_coding
```
genotype files - 
```
snp_ID 24788 24831 24840 24870 24893 24928 ...
chr1:632373 2 2 2 2 2 2 ...
chr1:782247 0 0 0 0 0.014 0 ...
```
expression files -
```
gene_id 24788   24831   24840   24870   24893   24928 ...
ENSG00000188976 -0.0473655760288239     0.0835636258125305      -0.850793421268463      0.117297530174255       -0.512760400772095      0.000488340854644775 ...
ENSG00000187961 0.0324062705039978      -0.0785085707902908     0.123972654342651       0.80484277009964        0.543906211853027       0.341333776712418 ...
```
## Building the Baseline Elastinet Models
### Output Files
### Plotting the Rho Averages
