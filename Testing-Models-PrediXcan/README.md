## Impute Expression 
- **impute_baseline_elasticnet.txt** - script that contains the file paths to run PrediXcan
  > * p - refers to the population that is being tested on
  > * p2 - refers to the population the model was built-in
  > * tiss - refers to the tissue type being tested in, refer to cross-tissue when the populations are in different tissue types 
- **qsub_loop.txt** - script that queues multiple jobs of PrediXcan
  > * contains list of populations to be looped through
- **genetic_pi1.R** - script that calculates the Spearmen correlation between predicted and actual expression
  > * NOTE: make the file pi1 - the script will not make the pi1 file if directory does not exist
- **prediXcan_input_from_dosages.txt** - script that takes in dosage file, removes the header, and zips the new prediXcan specific dosage file 
## Refer to the [wiki](https://github.com/chrisnguyen11/TOPMed-Expression/wiki/00-Testing-Models-(Running-PrediXcan)) for more detail and how to use! 
