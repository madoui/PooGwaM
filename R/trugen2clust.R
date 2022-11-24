#' Assigning true genotypes to pools
#'
#' This will assign the true genotypes to the cluster based on phenotypes
#'
#' @param clusters cluster of individual
#' @param sim a PhenotypeSimulator object
#' @return genotypes assigned to cluster
#'
#' @examples
#' sim = sim2snp(100,1000,10)
#' pools = phen2pool(sim, 3, 3)
#' GenotypesCluster = trugen2clust (pools, sim)
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @export
trugen2clust <-function(clusters, sim){
  cls_genotype<-function(clusters, genotype){
    k = length(clusters)
    l = list()
    for(i in 1:k){
      #l[[i]] = genotype[ row.names(genotype) == clusters[[i]] ,]
      l[[i]]= genotype %>% filter(row.names(genotype) %in% clusters[[i]])
    }
    return(l)
  }
  genotypes=as.data.frame(sim[["rawComponents"]][["genotypes"]][["genotypes"]])
  clusters_genotype = cls_genotype(clusters, genotypes)
  return(clusters_genotype)
}


