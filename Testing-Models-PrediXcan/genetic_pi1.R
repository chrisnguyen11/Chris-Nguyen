##################################################################################################

#SETUP ENVIRONMENT
cat("SETUP ENVIRONMENT\n")
##################################################################################################

library(dplyr)
library(qvalue)
library(data.table)
library(ggplot2) 


##################################################################################################

#DEFINE FUNCTIONS
cat("DEFINE FUNCTIONS\n")
##################################################################################################


"%&%" = function(a,b) paste(a,b,sep="")

match_orders<-function(rowmatrix,colmatrix, idcol=2){
  sampleOrderColMat<-colmatrix[,idcol] %>% unlist() %>% unname()
  col1RowMat<-colnames(rowmatrix)[1]
  rowmatrix<-rowmatrix %>% select(one_of(c(col1RowMat,sampleOrderColMat)))
  return(rowmatrix)
}

get_gene_corr<-function(gene_name,colmatrix,rowmatrix,skipcol=1)
{
  #gene name is the name of the gene
  #colmatrix is a matrix containing genes as columns
  #rowmatrix is a matrix containing genes as rows
  #skipcol = The first n columns of the row matrix do not contain values, skip these for correlation test  

  predicted_exp<-colmatrix %>% 
    select(gene_name) %>% ##select the gene
    unlist() %>% unname()
  
  measured_exp<-rowmatrix %>% 
    filter( gene_id == gene_name) %>% 
    select((skipcol+1):ncol(rowmatrix)) %>%
    unlist() %>% unname()
  
  correlation<-cor.test(measured_exp,predicted_exp, method = "spearman")
  
  expression_corr<-list()
  expression_corr[["gene_id"]]<-gene_name
  expression_corr[["estimate"]]<-correlation$estimate
  expression_corr[["p.value"]]<-correlation$p.value
  
  return(expression_corr)
}


##################################################################################################

#SET GLOBAL VARIABLES
cat("SET GLOBAL VARIABLES\n")
##################################################################################################

tiss="Tcell" 

pop_list<-c("ALL","AFA","CAU","HIS")

#tiss_list<-c("PBMC","Mono","Tcell")

R2_filtering<-c(-1)

columns<-expand.grid(pop_list,R2_filtering) %>% 
  arrange(Var1,Var2) %>%
  mutate(columns=Var1 %&% Var2) %>% rename(pop=Var1,R2=Var2)

col_list<-columns %>% select(columns) %>% unlist() %>% unname()

pi1_matrix<-matrix(NA, nrow = length(pop_list), ncol = length(col_list))
rownames(pi1_matrix)<-pop_list
colnames(pi1_matrix)<-col_list

##################################################################################################

#READ IN & PROCESS DATA
cat("READ IN & PROCESS DATA\n")
##################################################################################################

for (obspop in pop_list)
  {#per each set of observed expression data
  
  observed_expression<-fread("zcat /home/chris/topmed/expression/" %&% obspop %&% "/" %&% tiss %&% "/expression/" %&% tiss %&% "_" %&% obspop %&% "_age_sex_exam_10PCS_PF_adj_rinv_TOPMED_expression10_peer_factor_adjusted.txt.gz", header = T, sep = '\t',stringsAsFactors = F) %>%
  rename_at(vars(1),function(x){return("gene_id")})
  
  for (i in 1:nrow(columns))
    {#How well does each model replicate
    
    predpop<-columns[i,1]
    pi1column<-columns[i,3]
    
    models<-fread("/home/chris/topmed/expression/baseline_elasticnet/imputed_results/" %&% tiss %&% "/" %&% tiss %&% "_" %&% predpop %&% "_unfiltered_model_" %&% obspop %&% "_geno_prediction_summary.txt",stringsAsFactors = F)
    names<-fread("/home/chris/topmed/expression/baseline_elasticnet/imputed_results/geno_prediction_summary_header.txt",header=F,stringsAsFactors = F)
    names<-c(unlist(names))
    colnames(models)<-names
    
    models<-models %>% select(gene_id, test_R2_avg)

    predicted_expression<-fread("/home/chris/topmed/expression/baseline_elasticnet/imputed_results/" %&% tiss %&% "/" %&% tiss %&% "_" %&% predpop %&% "_unfiltered_model_" %&% obspop %&% "_geno_predicted_expression.txt",header=T,stringsAsFactors = F)
    
    obsGenes<-unlist(observed_expression$gene_id)
    predGenes<-colnames(predicted_expression)
    
    gene_list<-data.frame(gene_id=intersect(obsGenes,predGenes),stringsAsFactors = F)

    models<-inner_join(gene_list,models, by = "gene_id")
    cat("Existing Gene Models:\n")
    str(models)

    ##################################################################################################
    
    #ANALYZE DATA
    cat("ANALYZE DATA\n")
    ##################################################################################################
    
    for(j in R2_filtering)
      {#Calculate pi1 at different model R2 thresholds
      
      cat("testing replication rate at R2 threshold of ",j,"\n")
      filtered_gene_R2<-models %>% 
        filter( test_R2_avg > j)
      
      if(dim(filtered_gene_R2)[1] == 0)
      {#check if there are any models that meet R2 threshold
        
        cat("WARNING: ",pi1column," had no models at R2 threshold ", j, " for genes in ",obspop,"\n")
#        cat("WARNING: possible at higher R2 thresholds, but unlikely at low thresholds unless there is a gene name mismatch. Check that your genenames are the same between files")
        cat("WARNING: ASSIGNING 0 as pi1\n")
        pi1_matrix[obspop,pi1column] <- 0
        next
        
      }
      predictive_correlations<-sapply(X=gene_list$gene_id,FUN=get_gene_corr,colmatrix=predicted_expression,rowmatrix=observed_expression,simplify=T,USE.NAMES = T)
      predictive_correlations<-data.frame(gene_id=unlist(predictive_correlations[1,]),estimate=unlist(predictive_correlations[2,]),p.value=unlist(predictive_correlations[3,]))
      
      str(predictive_correlations)
      
      filtered_gene_R2<- filtered_gene_R2 %>% inner_join(predictive_correlations, by = "gene_id")
      
      ##################################################################################################
      
      #WRITE OUT CORRELATION DATA
      cat("WRITE OUT CORRELATION DATA\n")
      ##################################################################################################
      
      fwrite(filtered_gene_R2,"/home/chris/topmed/expression/baseline_elasticnet/imputed_results/" %&% tiss %&% "/" %&% tiss %&% "_" %&% predpop %&% "_R2_" %&% j %&% "_model_" %&% obspop %&% "_expression_spearman_correlation.txt",col.names = T,sep='\t')
      corr_pvals<-filtered_gene_R2$p.value
      qobjCorr <- tryCatch({qvalue(p = corr_pvals)},
                           error=function(cond){
                             cat("Error: ", pi1column, "and ",obspop, " caused error to qvalue, assigning 0 as pi1\n")
                             cond
                           })
      if(inherits(qobjCorr, "error")) {
        pi1_matrix[obspop,pi1column] <- 0
        next
      }
      pi1<- 1 - qobjCorr$pi0
      pi12<- signif(pi1,4)
      pi1_matrix[obspop,pi1column] <- pi12
    }#close
  }
}
##################################################################################################

#WRITE OUT SUMMARY DATA
cat("WRITE OUT SUMMARY DATA")
##################################################################################################

pi1_matrix<-data.table(pi1_matrix)
row.names(pi1_matrix)<-pop_list
fwrite(pi1_matrix, "/home/chris/topmed/expression/baseline_elasticnet/imputed_results/" %&% tiss %&% "/pi1/" %&% tiss %&% "_pi1_R2-1.txt", sep = '\t', col.names=T,row.names = T)
