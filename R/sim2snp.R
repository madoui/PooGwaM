#' Simulate phenotypes and genotypes
#'
#' PhenotypeSimulator wrapper to generates 2-dimensional phenotypes  (so, P=2)
#' with 3 components: genetic effects, polygenic effects, and environmental noise.
#' There is no covariaye effect (hence 'delta'=0) and no correlated noise (hence 'rho'=0).
#' Total number of individuals \code{n}, default=2000,
#' total number of SNPs \code{tNrSNP}, default=10000,
#' Number of causal SNPs \code{cNrSNP}, default=10.
#' Total genetic heritability \code{genVar}, default=60%,
#' Note that SNP heritability is controlled by '\code{genVar} * \code{h2s}' and \code{h2s}=0.5/0.6 by default
#' We randomly select 2 SNPs out of 10  to be causal for both the traits
#' (i.e. pleitropic SNPs or SNPs having shared effect);
#' and this quantity is controlled by '1-pIndependerGenetic'.
#' So we set '1-pIndependerGenetic'=2/10=0.2 or 'pIndependerGenetic'=0.8.
#' 'pTraitIndependentGenetic = 0.5' means that the remaining 8 SNPs are randomly selected to be causal for only one trait (i.e 1 out of 2 traits. Hence 1/2=0.5 appears)
#''SNPfrequencies' denote the set of MAFs. MAF of a particular SNP is randomly selected from this set of MAFs.
#'
#' @param n Numeric vector.
#' @param tNrSNP Numeric vector
#' @param cNrSNP Numeric vector
#' @return a PhenotypeSimulator object
#'
#' @examples
#' sim = sim2snp(100, 100, 10)
#' @export
sim2snp<-function(n, tNrSNP,cNrSNP){
  #generating phenotypes
  phenotype=PhenotypeSimulator::runSimulation(N=n, P=2, genVar=0.6,
                                              h2s=0.5/0.6, rho=0, delta=0, tNrSNP=tNrSNP,
                                              cNrSNP=cNrSNP, SNPfrequencies=c(0.05, 0.1,0.2),
                                              pIndependentGenetic = 0.8,
                                              pTraitIndependentGenetic = 0.5, seed=16 )
  #identify causal SNP (trait-specific and peliotrpic)
  fixedGen=phenotype[["phenoComponentsIntermediate"]][["genFixed"]]
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
  snps_causal_trait1=unname(sapply(snps_causal_trait1, getindependent)) # these are SNPs only causal for Trait 1

  # SNPs which are only causal for Trait 2:-
  snps_causal_trait2=unname(sapply(snps_causal_trait2, getindependent)) # these are SNPs only causal for Trait 2

  # The totality of 10 causal SNPs:-
  causal_snps=colnames(fixedGen[["cov"]])

  # the N by 2 phenotype matrix(each column of the matrix is obs. on one trait):-

  Y=as.data.frame(phenotype[["phenoComponentsFinal"]][["Y"]])

  # plot scatter plot of the two traits:-

  trait_1 = Y$Trait_1
  trait_2 = Y$Trait_2
  g = ggplot2::ggplot(data=Y, ggplot2::aes(x=trait_1, y=trait_2)) + ggplot2::geom_point()
  plot(g)
  return(phenotype)
}








