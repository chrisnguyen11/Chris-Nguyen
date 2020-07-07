## Preparing Files 
The Geuvadis genotype files are in the zipped vcf format (vcf.gz), and the expression files are in the zipped gct format (gct.gz). 
1. We need to unzip the genotype files to allow for ease of use.
2. We will only be using the samples with avalible expression. 
- a. The expression files are named by the samples so we pulled all the names of the files and thus the samples with expression. 
```
  ls > temp
  sed 's/.rnaseqccounts.gene_rpkm.gct.gz//g' temp > rna_sample_list.txt 
```
```
  HG00096
  HG00097
  HG00099
  ...
```
- b. Going back to the genotype files, the header was pulled and then numbered.

```
  head -n 20 ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf | tail -n 1 | sed "s/\t/\n/g" > temp
  awk '{print $s " " NR}' temp > header_vcf.txt
```
```
  #CHROM 1
  POS 2
  ID 3
  ...
```
- c. Using the rna_sample_list.txt and header_vcf.txt, we can get the columns of samples in the vcf files with both geotype and expression data.
```
  cat header_vcf.txt rna_sample_list.txt | sort | uniq -d | awk '{print $2}' > col_rna_sample_list.txt
```
```
  10 11 12 ...
```  
- d. Using this list we are able to go back to the genotype files and create new vcf files with samples who also have expression data. Similar to the scripts: 

    - vcf_annotations.txt - writes the header portion of the vcf file to the new vcf file 

    - vcf_pull_RNA_samples.txt - writes the genotype data portion of the vcf file to the new vcf file 

3. To get the samples seperated by populations, Ryan had a file that contained the sample and their population:
```
  HG00096	GBR
  HG00097	GBR
  HG00099	GBR
...
```
- a. The samples seprated by population where then checked aganist the samples with expression data. Steps were similar to steps 2a - 2c: 
  - NOTE: Some samples were lost. 
```
  HG00096	GBR 10
  HG00097	GBR 11
  HG00099	GBR 12
```      
- b. Then the individual population list were used to create new vcf files with samples who also have population data.
    - This was repeated for all the populations: ALL, CEU, FIN, GBR, TIS, YRI
4. With the new vcf files, we need to caululate the minor allele frequency (MAF) and filter by 0.01 MAF.

## Running GWAS-QC
To run GWAS-QC, I used [Ryan's pipline](https://github.com/RyanSchu/gwasqc_pipeline).
The first part of the QC is running 

~/gwasqc/01MissingnessFiltering.txt -b ~/geuvadis/vcf/ALL.chr1 -g 0.01 -h 0.01 -o ~/geuvadis/vcf/QC


~/gwasqc/02RelatednessFiltering.txt -b ~/geuvadis/vcf/QC/missingness_hwe_steps/05filtered_HWE -o ~/geuvadis/vcf/QC/ --rel 0.05


~/gwasqc/03MergeHapMap_TOPMED.txt -b ~/geuvadis/vcf/QC/relatedness_steps/05without_relateds -h ~/HAPMAP3/hg38/hg38_topmed/ -o ~/geuvadis/vcf/QC/ --maf 0.01


## Making Dosages
