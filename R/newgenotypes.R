#' Regenerate genotype from allele freuencies
#'
#' The goall is to simulate new genotypes from the observed allele frequencies in the pools
#'
#' @param pool_size the desired pool size
#' @param clusters the pools of individuals
#' @param sim a PhenotypeSimulator Object
#' @param allele_freq allele frequencies
#' @return a dataframe of new genotypes
#'
#' @importFrom stats rbinom
#'
#' @examples
#' sim = sim2snp(100,100,10)
#' pools = phen2pool(sim, 3, 3)
#' GenotypesCluster = trugen2clust (pools, sim)
#' alleleFreq = pool2freq (pools, GenotypesCluster, sim)
#' new_genotypes = newgenotypes (GenotypesCluster, pools, sim, alleleFreq)
#' @export
#'

newgenotypes <- function(pool_size, clusters, sim, allele_freq){
  ind_percluster=sapply(pool_size, nrow)
  # 'clusters_genotype_re' below is a list where i th element of the list
  # will be a dataframe with first column all i's (that is, the cluster number)
  # and rest columns will denote regenerated allele frequencies for all the SNPs
  # in cluster i
  clusters_genotype_re=list()
  set.seed(156) # set this seed to have the reproducible regenerated data
  for(i in 1:length(clusters)){
    a=ind_percluster[i]
    y=rep(i,a)
    SNPs = sim[["setup"]][["id_snps"]] # denote the IDs of the SNPs.
    d=matrix(NA, nrow=a, ncol=length(SNPs))
    for(j in 1:length(SNPs)){
      d[,j]= rbinom(a,2, allele_freq[i,j])

    }

    clusters_genotype_re[[i]]=cbind(y,d)

  }
  # the object 'data_regenerate' below is the regenerated genotype data from the pooled
  # allele frequencies. The first column of this data denotes the cluster numbers.
  # So we can use this data for classification purpose, say, random forest.
  # But be sure to shuffle the data before applying to random forest
  data_regenerate=clusters_genotype_re[[1]]
  for(i in 2:length(clusters_genotype_re)){
    data_regenerate=rbind(data_regenerate,clusters_genotype_re[[i]] )
  }

  data_regenerate=as.data.frame(data_regenerate)

  colnames(data_regenerate)[2:ncol(data_regenerate)]=SNPs
  return(list("regen"=data_regenerate,"freq_re"=clusters_genotype_re))
}

