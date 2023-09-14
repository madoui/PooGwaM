#' Chisquare
#'
#' apply classical chi-square
#'
#' @param Freq data frame of genotypes
#' @return list
#' @importFrom stats chisq.test
#' @examples
#' sim = PhenoSim (500, 100, 10, 0.7, 2)
#' clusters = quantile_clustering (data.frame(sim$phenotypes), 4)
#' Freq = compute_group_MAFs( sim$genotypes, as.factor(clusters) )
#' pv <- Chisquare (Freq)
#'
#'
#' @export

Chisquare <- function (Freq){
  pvalues<-apply(Freq, 2, function(x){chisq.test(Freq, as.factor(rownames(Freq)))$p.value})
  return (pvalues)
}
