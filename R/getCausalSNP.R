#' Get causal SNPs from the simulation
#'
#' @param simulation PhenotypeSimulator object
#' @return list
#' @import ggplot2
#'
#' @examples
#' sim = sim2snp(100,100,10)
#' causalSNP = getCausalSNP(sim)
#' @export
#'
getCausalSNP<-function(simulation){
  fixedGen=simulation[["phenoComponentsIntermediate"]][["genFixed"]]
  snpscombined=as.data.frame(t(fixedGen$cov_effect))
  pleitropicsnps=rownames(snpscombined[!snpscombined$Trait_1==0 & !snpscombined$Trait_2==0,])
  snps_causal_trait1=rownames(snpscombined[!snpscombined$Trait_1==0 & snpscombined$Trait_2==0,])
  snps_causal_trait2=rownames(snpscombined[snpscombined$Trait_1==0 & !snpscombined$Trait_2==0,])
  getshared=function(char){
    return( substring(char, 28, nchar(char))       )
  }
  # 'pleitropicsnps' below is a vector of the two pleitropic SNPs
  pleitropicsnps=sapply(pleitropicsnps, getshared)
  getindependent=function(char){
    return( substring(char, 38, nchar(char))       )
  }
  # SNPs which are causal only for Trait 1:-
  snps_causal_trait1 = unname(sapply(snps_causal_trait1, getindependent)) # these are SNPs only causal for Trait 1

  # SNPs which are only causal for Trait 2:-
  snps_causal_trait2 = unname(sapply(snps_causal_trait2, getindependent)) # these are SNPs only causal for Trait 2

  # The totality of 10 causal SNPs:-
  causal_snps=colnames(fixedGen[["cov"]])

  # the N by 2 phenotype matrix(each column of the matrix is obs. on one trait):-

  Y=as.data.frame(simulation[["phenoComponentsFinal"]][["Y"]])

  # plot scatter plot of the two traits:-

  trait_1 = Y$Trait_1
  trait_2 = Y$Trait_2
  g = ggplot2::ggplot(data=Y, ggplot2::aes(x=trait_1, y=trait_2)) + ggplot2::geom_point()
  plot(g)
  causals= list("SNP_causal_trait_1_only" = snps_causal_trait1,
                "SNP_causal_trait_2_only" = snps_causal_trait2,
                "Total" = causal_snps)
  return(causals)
}
