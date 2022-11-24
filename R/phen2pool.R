#' Cluster individuals based on their phenotypes into equal-sized-pool for pool-seq
#'
#' The goal is to cluster the individuals based on their phenotypes,
#' Trait 1 is divided into \code{K} quantiles
#' Trait2 is divided into \code{L} quantiles
#'
#' @param sim a PhenotypeSimulator object.
#' @param K Numeric vector
#' @param L Numeric vector
#' @return list
#'
#' @examples
#' sim = sim2snp(100,100,10)
#' pools = phen2pool(sim, 3, 3)
#' @export
#'
phen2pool <- function(sim, K, L){
  # we want the area between the quantiles to be same (1/K in this case)
  Y = as.data.frame(sim[["phenoComponentsFinal"]][["Y"]])
  probs_1=seq(0,1, by=1/K)
  probs_2=seq(0,1, by=1/L)

  # quantiles of Trait-1 at probability specified in 'probs_1':-
  quantile_1 = stats::quantile(Y$Trait_1, probs=probs_1)

  # Below 'clusters_trait1' is the list of K clusters (almost same sizes)
  # obtained by quantile based clustering of Trait_1:-
  clusters_trait1=list()
  for(i in 1:(length(probs_1)-1)){
    if(i < (length(probs_1)-1)){
      clusters_trait1[[i]]=Y[ Y$Trait_1 >= quantile_1[i] & Y$Trait_1 < quantile_1[i+1] , ]
    }
    else{
      clusters_trait1[[i]]=Y[ Y$Trait_1 >= quantile_1[i] & Y$Trait_1 <= quantile_1[i+1] , ]
    }
  }
  getclusters_trait2=function(df){
    temp_list=list()
    temp_quantile=stats::quantile(df$Trait_2, probs=probs_2)
    for(i in 1:(length(probs_2)-1)){
      if(i< (length(probs_2)-1)){
        temp_list[[i]]=df[ df$Trait_2 >= temp_quantile[i] & df$Trait_2 < temp_quantile[i+1],  ]
      }
      else{
        temp_list[[i]]=df[ df$Trait_2 >= temp_quantile[i] & df$Trait_2 <= temp_quantile[i+1], ]
      }
    }
    return(temp_list)
  }

  clusters_traits=list()
  for(i in 1:length(clusters_trait1)){
    temp=getclusters_trait2(clusters_trait1[[i]])
    clusters_traits=append(clusters_traits, temp)
  }

  # finally all the K*L clusters presenetd in terms of individual ids:-

  clusters=lapply(clusters_traits, rownames) # it will be a list, each denoting the cluster of individuals' ID belonging to a particular cluster
  return(clusters)
}










