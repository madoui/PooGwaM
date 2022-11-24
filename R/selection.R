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
#' sim = sim2snp(100,1000,10)
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

  # LASSO + Stability selection for Trait 1:-
  # Below We use 'Original Stability Selection' with tau=0.8 and PFER=0.5.
  # For using CPSS, just change 'sampling.type=SS' and 'assumption=r-concave' in the function below.


  stab.lasso=stabsel(x=x, y=y1, intercept=F, fitfun=glmnet.lasso,
                     cutoff=0.8, PFER=0.5, assumption="none",
                     sampling.type="MB")

  temp1=stab.lasso[["selected"]] # vector of selected variables (or, IDs of selected SNPs)
  length(temp1) # total number of discovery
  #sum(parse_number(causal_snps)%in%temp1) # no. of causal SNPs captured by performing LASSO+ Stab. Selection to Trait 1
  #sum(parse_number(pleitropicsnps) %in% temp1) # no. of 'pleitropic causal SNPs' captured by performing LASSO+ Stab. Selection to Trait 1
  #sum(parse_number(snps_causal_trait1) %in% temp1) # no. of 'causal for Trait 1 only SNPs' captured by performing LASSO+ Stab. Selection to Trait 1

  # LASSO+ Stab. Selection for Trait 2:-

  stab.lasso=stabsel(x=x, y=y2, intercept=F, fitfun=glmnet.lasso,
                     cutoff=0.8, PFER=0.5, assumption="none",
                     sampling.type="MB")

  temp2=stab.lasso[["selected"]] # vector of selected variables (or, IDs of selected SNPs)
  #length(temp2) # total number of discovery
  #sum(parse_number(causal_snps)%in%temp2) # no. of causal SNPs captured by performing LASSO+ Stab. Selection to Trait 2
  #sum(parse_number(pleitropicsnps) %in% temp2) # no. of 'pleitropic causal SNPs' captured by performing LASSO+ Stab. Selection to Trait 2
  #sum(parse_number(snps_causal_trait1) %in% temp2) # no. of 'causal for Trait 2 only SNPs' captured by performing LASSO+ Stab. Selection to Trait 2


  temp=unique(c(unname(temp1), unname(temp2)) ) # this is the total combined list of discovery by Trait 1 and Trait 2


  #causal_snp_trait1 = parse_number(causal_snps) %in% temp1 # total number of causal discovery by Trait 1 and Trait 2 combined
  #causal_snp_trait2 = parse_number(causal_snps) %in% temp2 # total number of causal discovery by Trait 1 and Trait 2 combined

  return(c("snp_trait1" = temp1, "snp_trait2" = temp2))

}
