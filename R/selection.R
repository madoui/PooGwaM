#' Detection of causal SNPs
#'
#' The goall is to discover causal SNP from regenrated genotypes data and multitrait data
#'
#' @param data phenotypes and genotypes data
#' @return a dataframe of causal SNP
#'
#' @import glmnet
#' @import stabs
#'
#' @examples
#' sim = sim2snp(100,100,10)
#' pools = phen2pool(sim,3,3)
#' GenotypesCluster = trugen2clust(pools,sim)
#' freq = pool2freq(pools,GenotypesCluster,sim)
#' regen = newgenotypes(GenotypesCluster,pools,sim,freq)
#' newgen_phen_comb = newgen2phen(pools, regen, sim)
#' causal_snps = selection (newgen_phen_comb)
#' @export
#'


selection <-function(data){
  n0 = nrow(data) #no. of observations==no. of individuals
  data = data[sample(1:n0, size=n0),] # shuffling the rows
  x = data[,3:ncol(data)] # the feature matrix
  y1 = data$y_1 # Trait 1 observations
  y2 = data$y_2 # Trait 2 observations

  stab.lasso=stabsel(x=x, y=y1, intercept=F, fitfun=glmnet.lasso,
                     cutoff=0.8, PFER=0.5, assumption="none",
                     sampling.type="MB")

  temp1=stab.lasso[["selected"]] # vector of selected variables (or, IDs of selected SNPs)

  # LASSO+ Stab. Selection for Trait 2:-

  stab.lasso=stabsel(x=x, y=y2, intercept=F, fitfun=glmnet.lasso,
                     cutoff=0.8, PFER=0.5, assumption="none",
                     sampling.type="MB")

  temp2=stab.lasso[["selected"]] # vector of selected variables (or, IDs of selected SNPs)

  temp=unique(c(unname(temp1), unname(temp2)) ) # this is the total combined list of discovery by Trait 1 and Trait 2

  return(list("snp_trait1" = temp1, "snp_trait2" = temp2, "snp_unique" = temp))

}
