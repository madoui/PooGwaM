#' poogwam
#'
#' GWAS on pool-seq data using several methods
#'
#' @param Method string
#' @param clusters pool ID
#' @param Freq data frame of allele frequencies
#' @param K interger number of pools
#' @param phenotypes data frame of phenotypes
#' @param eta default=1
#' @param nb.batch float default=1
#' @param max_iter float default=1e3
#' @param tol float default=1e-8
#' @param seed default=1
#' @param trace boolean default=FALSE
#'
#' @return numeric vector of SNP effects
#'
#' @importFrom nnet class.ind
#' @importFrom stats rnorm
#'
#' @examples
#' sim = PhenoSim ( 1000, 100, 10, 0.6, 3 )
#' K = 4
#' clusters = quantile_clustering ( data.frame ( sim$phenotypes ), K )
#' Freq = compute_group_MAFs( sim$genotypes, as.factor(clusters) )
#' pv = poogwam ("poogwam", Freq, sim$phenotypes[,1], clusters )
#'
#' @export
#'
poogwam <- function(Method = "poogwam",
                    Freq,
                    phenotypes,
                    clusters,
                    K=0,
                    eta = 1,
                    max_iter = 1000,
                    tol = 1e-8,
                    seed = 1,
                    nb.batch = 1,
                    trace = FALSE
                    ){
  if (Method == "poogwam"){
    res = adaptRidge( Freq, phenotypes, clusters, eta = eta,
                      max_iter = max_iter,
                      tol = tol,
                      seed = seed,
                      nb.batch = nb.batch,
                      trace = trace )
    beta = res$beta
    pvalues<- pvalues.from.test.stat(beta)
    return (pvalues)
  }

  if (Method == "GWALPHA"){
    return (GWalpha(phenotypes, clusters, Freq, K))
  }

  #if (Method == "Chi2"){
   # return (Chisquare (sgenotypes, clusters))
  #}
}
