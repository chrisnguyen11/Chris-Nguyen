Creating Geuvadis High Coverage vcfs 
- Only using SNPs found in the TOPMed dosage files 
  - Using the script 'make_pos.txt' to get the second column from the TOPMed dosage files (postion column) 
  - Gathers the postion of all SNPs across all tissues and populations 
    - Using the command 'cat *chr22* | sort | uniq > topmed_chr22_pos.txt' to call all the individal tissue population SNP files and only include unique SNPs
    - Using the command 'sed -i 's/chr22://g' topmed_chr22_pos.txt' to remove the 'chr22:' from the postions 
    
