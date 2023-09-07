#' Chisquare
#'
#' apply classical chi-square
#'
#' @param genotypes data frame of genotypes
#' @param clusters numeric vector of the clusters
#' @return list
#' @importFrom stats chisq.test
#' @examples
#' sim = PhenoSim (1000, 100, 10, 0.6, 5)
#' clusters = quantile_clustering (data.frame(sim$phenotypes), 7)
#' pv <- Chisquare (sim$genotypes, clusters)
#'
#'
#' @export

Chisquare <- function (genotypes, clusters){
  pvalues<-apply(genotypes,2,function(x){chisq.test(table(clusters,as.factor(x)))$p.value})
  return (pvalues)
}
