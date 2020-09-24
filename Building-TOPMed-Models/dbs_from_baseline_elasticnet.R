library(dplyr)
library(RSQLite)
library(qvalue)
"%&%" <- function(a,b) paste(a,b, sep='')

driver <- dbDriver('SQLite')

pop<-"HIS"
tiss<-"PBMC"

# Extra table ----
model_summaries <- read.table('/home/chris/topmed/expression/' %&% pop %&% '/' %&% tiss %&% '/baseline_elasticnet/baseline_'
                              %&% tiss %&% '_' %&% pop %&% '_chr1_model_summaries.txt', header = T, stringsAsFactors = F)
tiss_summary <- read.table('/home/chris/topmed/expression/' %&% pop %&% '/' %&% tiss %&% '/baseline_elasticnet/baseline_'
                           %&% tiss %&% '_' %&% pop %&% '_chr1_tiss_chr_summary.txt', header = T, stringsAsFactors = F)

n_samples <- tiss_summary$n_samples

for (i in 2:22) {
  model_summaries <- rbind(model_summaries,
                           read.table('/home/chris/topmed/expression/' %&% pop %&% '/' %&% tiss %&% '/baseline_elasticnet/baseline_'
                                      %&% tiss %&% '_' %&% pop %&% '_chr' %&% as.character(i) %&% '_model_summaries.txt', header = T, stringsAsFactors = F))
  tiss_summary <- rbind(tiss_summary,
                        read.table('/home/chris/topmed/expression/' %&% pop %&% '/' %&% tiss %&% '/baseline_elasticnet/baseline_'
                                   %&% tiss %&% '_' %&% pop %&% '_chr' %&% as.character(i) %&% '_tiss_chr_summary.txt', header = T, stringsAsFactors = F))
}

pvalues<-model_summaries$zscore_pval
qvalues<-qvalue(pvalues)
model_summaries <- rename(model_summaries, 
                          gene = gene_id, 
                          genename = gene_name,
                          n.snps.in.window = n_snps_in_window,
                          n.snps.in.model = n_snps_in_model,
                          pred.perf.R2 = rho_avg_squared,
                          pred.perf.pval = zscore_pval
) %>% mutate(pred.perf.qval = qvalues$qvalues)


conn <- dbConnect(drv = driver, '/home/chris/topmed/expression/' %&% pop %&% '/' %&% tiss %&% '/dbs/' %&% tiss %&% '_' %&% pop %&% '_WG_unfiltered.db')
dbWriteTable(conn, 'extra', model_summaries, overwrite = TRUE)
dbGetQuery(conn, "CREATE INDEX gene_model_summary ON extra (gene)")

# Weights Table -----
weights <- read.table('/home/chris/topmed/expression/' %&% pop %&% '/' %&% tiss %&% '/baseline_elasticnet/baseline_'
                      %&% tiss %&% '_' %&% pop %&% '_chr1_weights.txt', header = T, stringsAsFactors = F)
for (i in 2:22) {
  weights <- rbind(weights,
                   read.table('/home/chris/topmed/expression/' %&% pop %&% '/' %&% tiss %&% '/baseline_elasticnet/baseline_'
                              %&% tiss %&% '_' %&% pop %&% '_chr' %&% as.character(i) %&% '_weights.txt', header = T, stringsAsFactors = F))
}
weights <- rename(weights,gene = gene_id, weight = beta, ref_allele = ref, eff_allele=alt)
dbWriteTable(conn, 'weights', weights, overwrite = TRUE)
dbGetQuery(conn, "CREATE INDEX weights_rsid ON weights (rsid)")
dbGetQuery(conn, "CREATE INDEX weights_gene ON weights (gene)")
dbGetQuery(conn, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")

# Sample_info Table ----
sample_info <- data.frame(n_samples = n_samples, population = pop, tissue = tiss)
dbWriteTable(conn, 'sample_info', sample_info, overwrite = TRUE)

# Construction Table ----
construction <- tiss_summary %>%
  select(chrom, cv_seed) %>%
  rename(chromosome = chrom)
dbWriteTable(conn, 'construction', construction, overwrite = TRUE)

### write models with no NA

conn2 <- dbConnect(drv = driver, '/home/chris/topmed/expression/' %&% pop %&% '/' %&% tiss %&% '/dbs/' %&% tiss %&% '_' %&% pop %&% '_WG_noNA_unfiltered.db')
complete_model_summaries<-model_summaries[complete.cases(model_summaries),]
dbWriteTable(conn2, 'extra', complete_model_summaries, overwrite = TRUE)
dbGetQuery(conn2, "CREATE INDEX gene_model_summary ON extra (gene)")

# Weights Table -----
weights <- read.table('/home/chris/topmed/expression/' %&% pop %&% '/' %&% tiss %&% '/baseline_elasticnet/baseline_'
                      %&% tiss %&% '_' %&% pop %&% '_chr1_weights.txt', header = T, stringsAsFactors = F)
for (i in 2:22) {
  weights <- rbind(weights,
                   read.table('/home/chris/topmed/expression/' %&% pop %&% '/' %&% tiss %&% '/baseline_elasticnet/baseline_'
                              %&% tiss %&% '_' %&% pop %&% '_chr' %&% as.character(i) %&% '_weights.txt', header = T, stringsAsFactors = F))
}
weights <- rename(weights, gene = gene_id, weight = beta, ref_allele = ref, eff_allele=alt)
dbWriteTable(conn2, 'weights', weights, overwrite = TRUE)
dbGetQuery(conn2, "CREATE INDEX weights_rsid ON weights (rsid)")
dbGetQuery(conn2, "CREATE INDEX weights_gene ON weights (gene)")
dbGetQuery(conn2, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")

# Sample_info Table ----
sample_info <- data.frame(n_samples = n_samples, population = pop, tissue = tiss)
dbWriteTable(conn2, 'sample_info', sample_info, overwrite = TRUE)

# Construction Table ----
construction <- tiss_summary %>%
  select(chrom, cv_seed) %>%
  rename(chromosome = chrom)
dbWriteTable(conn2, 'construction', construction, overwrite = TRUE)

### Write only significant models
conn3 <- dbConnect(drv = driver, '/home/chris/topmed/expression/' %&% pop %&% '/' %&% tiss %&% '/dbs/' %&% tiss %&% '_' %&% pop %&% '_WG_rho0.01_zpval0.05.db')
signif_models<-filter(complete_model_summaries, pred.perf.pval < 0.05 , rho_avg > 0.01)

dbWriteTable(conn3, 'extra', signif_models, overwrite = TRUE)
dbGetQuery(conn3, "CREATE INDEX gene_model_summary ON extra (gene)")

# Weights Table -----
weights <- read.table('/home/chris/topmed/expression/' %&% pop %&% '/' %&% tiss %&% '/baseline_elasticnet/baseline_'
                      %&% tiss %&% '_' %&% pop %&% '_chr1_weights.txt', header = T, stringsAsFactors = F)
for (i in 2:22) {
  weights <- rbind(weights,
                   read.table('/home/chris/topmed/expression/' %&% pop %&% '/' %&% tiss %&% '/baseline_elasticnet/baseline_'
                              %&% tiss %&% '_' %&% pop %&% '_chr' %&% as.character(i) %&% '_weights.txt', header = T, stringsAsFactors = F))
}
weights <- rename(weights, gene = gene_id, weight = beta, ref_allele = ref, eff_allele=alt)
weights<- weights %>% filter ( gene %in% signif_models$gene)
dbWriteTable(conn3, 'weights', weights, overwrite = TRUE)
dbGetQuery(conn3, "CREATE INDEX weights_rsid ON weights (rsid)")
dbGetQuery(conn3, "CREATE INDEX weights_gene ON weights (gene)")
dbGetQuery(conn3, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")


# Sample_info Table ----
sample_info <- data.frame(n_samples = n_samples, population = pop, tissue = tiss)
dbWriteTable(conn3, 'sample_info', sample_info, overwrite = TRUE)

# Construction Table ----
construction <- tiss_summary %>%
  select(chrom, cv_seed) %>%
  rename(chromosome = chrom)
dbWriteTable(conn3, 'construction', construction, overwrite = TRUE)

### write table with no NAs and no empty models

conn4 <- dbConnect(drv = driver, '/home/chris/topmed/expression/' %&% pop %&% '/' %&% tiss %&% '/dbs/' %&% tiss %&% '_' %&% pop %&% '_WG_noNA_noEmptyModels_unfiltered.db')
complete_model_summaries<-model_summaries[complete.cases(model_summaries),]
dbWriteTable(conn4, 'extra', complete_model_summaries, overwrite = TRUE)
dbGetQuery(conn4, "CREATE INDEX gene_model_summary ON extra (gene)")

# Weights Table -----
weights <- read.table('/home/chris/topmed/expression/' %&% pop %&% '/' %&% tiss %&% '/baseline_elasticnet/baseline_'
                      %&% tiss %&% '_' %&% pop %&% '_chr1_weights.txt', header = T, stringsAsFactors = F)
 for (i in 2:22) {
   weights <- rbind(weights,
                    read.table('/home/chris/topmed/expression/' %&% pop %&% '/' %&% tiss %&% '/baseline_elasticnet/baseline_'
                               %&% tiss %&% '_' %&% pop %&% '_chr' %&% as.character(i) %&% '_weights.txt', header = T, stringsAsFactors = F))
}
weights <- rename(weights, gene = gene_id, weight = beta, ref_allele = ref, eff_allele=alt)
dbWriteTable(conn4, 'weights', weights, overwrite = TRUE)
dbGetQuery(conn4, "CREATE INDEX weights_rsid ON weights (rsid)")
dbGetQuery(conn4, "CREATE INDEX weights_gene ON weights (gene)")
dbGetQuery(conn4, "CREATE INDEX weights_rsid_gene ON weights (rsid, gene)")
 
# Sample_info Table ----
sample_info <- data.frame(n_samples = n_samples, population = pop, tissue = tiss)
dbWriteTable(conn4, 'sample_info', sample_info, overwrite = TRUE)

# Construction Table ----
construction <- tiss_summary %>%
   select(chrom, cv_seed) %>%
   rename(chromosome = chrom)
dbWriteTable(conn4, 'construction', construction, overwrite = TRUE)
