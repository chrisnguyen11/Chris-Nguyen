setwd("~/Desktop/Wheeler_Lab/mount/chris/topmed/expression/CAU/Tcell/plots_baseline")
gene_ID_vs_rho_avg <- read.table("gene_ID_vs_rho_avg.txt", header = TRUE, sep = " ")
attach(gene_ID_vs_rho_avg)
names(gene_ID_vs_rho_avg)
rho_avg
hist(rho_avg, breaks = 18, xlim = c(-.3,1),ylim = c(0,5000), main = "CAU Tcell rho_avg")
