#' statgenGWASwrap
#'
#' a wrapper for the statgenGWAS package to run on simulated data
#'
#' @param genotypes data frame of genotypes
#' @param phenotypes data frame of phenotypes

#' @return statgenGWAS p-values
#'
#' @importFrom statgenGWAS createGData runSingleTraitGwas
#'
#' @examples
#' sim = PhenoSim ( 100, 100, 10, 0.6, 2 )
#' pv = statgenGWASwrap ( sim$genotypes, sim$phenotypes )
#'
#' @export

statgenGWASwrap <- function (genotypes, phenotypes){
  map <- data.frame( pos = seq(1, SNP), chr = rep(1,SNP) )
  rownames(map) <- colnames(genotypes)
  pheno4statgenGWAS = data.frame(genotype = rownames(genotypes), trait = as.numeric(phenotypes))
  gDataDrops <- createGData(geno = data.frame(genotypes),
                            map = map,
                            pheno = pheno4statgenGWAS)

  # Single traits GWAS analysis
  GWASDrops <- runSingleTraitGwas(gData = gDataDrops, traits ="trait", thrType = "fdr", pThr=0.3)
  return (GWASDrops$GWAResult$pheno4statgenGWAS$pValue)
}
