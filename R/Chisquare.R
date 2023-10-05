#' Chisquare
#'
#' apply classical chi-square
#'
#' @param Freq data frame of genotypes
#' @param clusters integers
#' @return list
#' @importFrom stats chisq.test
#' @examples
#' sim = PhenoSim (500, 100, 10, 0.7, 2)
#' clusters = quantile_clustering (data.frame(sim$phenotypes), 4)
#' Freq = compute_group_MAFs( sim$genotypes, as.factor(clusters) )
#' pv <- Chisquare (Freq, clusters)
#'
#'
#' @export

Chisquare <- function (Freq, clusters){
  cluster_size = table(clusters)
  freq2alleleCount <- function(SNPfreq, clusters_size){
    # regenerate allele counts from frequencies assuming Hardy Weinberg equilibrium
    allele_counts = matrix(rep(0, 3*length(cluster_size)), nrow=3)
    allele_counts[1,] = c(round(2*clusters_size*SNPfreq^2))
    allele_counts[2,] = c(round(2*clusters_size*SNPfreq*(1-SNPfreq)))
    allele_counts[3,] = c(round(2*clusters_size*(1-SNPfreq)^2))
    return(chisq.test(allele_counts)$p.value)
  }
  pvalues = apply( Freq, 2, FUN = function(x) freq2alleleCount (x, cluster_size))
  return (p.adjust(pvalues, "BH"))
}
