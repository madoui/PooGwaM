#' Simulate phenotypes and genotypes
#'
#' PhenotypeSimulator wrapper to generates p-dimensional phenotypes
#' with 3 components: genetic effects, polygenic effects, and environmental noise.
#' There is no covariaye effect (hence 'delta'=0) and no correlated noise (hence 'rho'=0).
#' Total number of individuals \code{n}, default=2000,
#' Number of phenotypes \code{p}, default = 2
#' total number of SNPs \code{tNrSNP}, default=10000,
#' Number of causal SNPs \code{cNrSNP}, default=10.
#' Total genetic heritability \code{gvar}, default=60%,
#' SNP heritability * \code{h2s}, default=50%
#' Ratio of Pleitropic SNPs \code{pleio}, default=0.2
#' We randomly select 2 SNPs out of 10  to be causal for both the traits
#' (i.e. pleitropic SNPs or SNPs having shared effect);
#' and this quantity is controlled by '1-pIndependerGenetic'.
#' So we set '1-pIndependerGenetic'=2/10=0.2 or 'pIndependerGenetic'=0.8.
#' 'pTraitIndependentGenetic = 0.5' means that the remaining 8 SNPs are randomly selected to be causal for only one trait (i.e 1 out of 2 traits. Hence 1/2=0.5 appears)
#''SNPfrequencies' denote the set of MAFs. MAF of a particular SNP is randomly selected from this set of MAFs.
#'
#' @param p Numeric vector
#' @param n Numeric vector.
#' @param tNrSNP Numeric vector
#' @param cNrSNP Numeric vector
#' @param gvar Numeric vector
#' @param snph2 Numeric vector
#' @param pleio numeric vector
#' @return a PhenotypeSimulator object
#'
#' @examples
#' sim = sim2snp(100, 100, 10)
#' @export
sim2snp<-function(n, tNrSNP, cNrSNP, p=2, gvar = 0.6, snph2 = 0.5, pleio = 0.2){
  #generating phenotypes
  phenotype=PhenotypeSimulator::runSimulation(N = n,
                                              P = p,
                                              genVar = gvar,
                                              h2s = snph2,
                                              rho=0,
                                              delta=0,
                                              tNrSNP=tNrSNP,
                                              cNrSNP=cNrSNP,
                                              SNPfrequencies=c(0.05, 0.1,0.2),
                                              pIndependentGenetic = 1-pleio,
                                              pTraitIndependentGenetic = 0.5, seed=16 )
  return(phenotype)
}








