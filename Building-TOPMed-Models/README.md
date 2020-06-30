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
snp_ID 24788 24831 24840 ...
chr1:632373 2 2 2  ...
chr1:782247 0 0 0  ...
```
expression files -
```
gene_id 24788   24831   24840   ...
ENSG00000188976 -0.0473655760288239     0.0835636258125305      -0.850793421268463      ...
ENSG00000187961 0.0324062705039978      -0.0785085707902908     0.123972654342651       ...
```
## Building the Baseline Elastic Net Models
### Output Files
covariances files -
```
GENE RSID1 RSID2 VALUE
ENSG00000188976 chr1:1193141 chr1:1193141 0.378501719487062
ENSG00000188976 chr1:1193141 chr1:1194826 0.37932820982495
```
model summaries - 
```
gene_id	gene_name	gene_type	alpha	n_snps_in_window	n_snps_in_model	lambda_min_mse	... 
ENSG00000188976	NOC2L	protein_coding	0.5	245	4	0.0605918843050675	...
ENSG00000187961	KLHL17	protein_coding	0.5	245	5	0.060013935237177	...
```
weights - 
```
gene_id	rsid	varID	ref	alt	beta
ENSG00000188976	chr1:1193141	1_1193141_C_T_b38	C	T	0.00206112053734138
ENSG00000188976	chr1:1194826	1_1194826_G_A_b38	G	A	0.0255065432819892
```
### Plotting the Rho Averages
One method for verifying if the models were built correctly is to plot the rho averages from the models summaries. For each tissue type in population, plot the rho averages and look at the shape of the frequencies of the rho averages. 
