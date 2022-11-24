#' Associate regenerated genotypes to phenotypes
#'
#'
#'
#' @param clusters Numeric vector.
#' @param regenerate Numeric vector
#' @param sim PhenotypeSimulator object
#' @return dataframe of associated genotypes
#'
#' @examples
#' sim = sim2snp(100,1000,10)
#' pools = phen2pool(sim,3,3)
#' GenotypesCluster = trugen2clust(pools,sim)
#' freq = pool2freq(pools,GenotypesCluster,sim)
#' regen = newgenotypes(GenotypesCluster,pools,sim,freq)
#' newgen_phen_comb = newgen2phen(pools, regen, sim)
#' @export

newgen2phen <- function(clusters,regenerate,sim){
  # clus_pheno below denote a list of K*L  clusters containing phenotype observations
  clus_pheno=list()
  Y=as.data.frame(sim[["phenoComponentsFinal"]][["Y"]])
  for(i in 1:length(clusters)){
    clus_pheno[[i]]=Y %>% filter(rownames(Y) %in% clusters[[i]])
  }
  # 'clusters_genphen_re' below is a list with i th element denoting the a dataframe
  # whose first two columns will be the phenotype observations for all the individuals
  # in i th pool, and rest of the columns denoting regenerated allele freq. of all SNPs
  # for i th pool.
  clusters_genphen_re=list()
  for(i in 1:length(regenerate$freq_re)){
    clusters_genphen_re[[i]]=cbind( clus_pheno[[i]]$Trait_1,
                                    clus_pheno[[i]]$Trait_2,
                                    regenerate$freq_re[[i]][, 2:ncol(regenerate$regen)])

  }
  # the object 'data_combined_re' below will be the regenerated genotype data from the pooled
  # allele frequencies. It is same as 'data_regenerate', except that now we have phenotype observations
  # in the first two columns. Whereas, in 'data_regenerate' there were no phenotypes, and first column there
  # was the cluster assignments.
  # We will use 'data_combined_re' for running LASSO.
  SNPs = sim[["setup"]][["id_snps"]] # denote the IDs of the SNPs.
  data_combined_re=clusters_genphen_re[[1]]
  for(i in 2:length(clusters_genphen_re)){
    data_combined_re=rbind(data_combined_re,clusters_genphen_re[[i]] )
  }
  data_combined_re=as.data.frame(data_combined_re)
  colnames(data_combined_re)=c('y_1','y_2', SNPs)
  return(data_combined_re)
}
