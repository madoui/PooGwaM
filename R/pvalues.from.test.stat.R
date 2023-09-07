#' pvalues.from.test.stat
#'
#' get p-values method from Fournier-Level et 2017
#'
#' @param test.stat numeric vector of SNP effects
#' @param alpha numeric vector of th error
#'
#' @return numeric vector
#'
#'
#' @examples
#' sim = PhenoSim (1000, 100, 10, 0.6, 5)
#' K = 7
#' clusters = quantile_clustering (data.frame(sim$phenotypes), K)
#' Freq= compute_group_MAFs(sim$genotypes,as.factor(clusters))
#' test.stat<-GWalpha(sim$phenotypes, clusters, Freq, K)
#' pv = pvalues.from.test.stat (test.stat)
#'
#' @importFrom stats pnorm
#'
#' @export

pvalues.from.test.stat <- function( test.stat, alpha = 0.01 ){
  shave <- function ( x, alpha = 0.0 ){
    limits <- quantile ( x, c(alpha/2,1-alpha/2) )
    return(x[(x > limits[1]) & (x <limits[2])])
  }
  # Estimate pvalues assuming that H0 test.stat  are gaussian distributed and
  #  that H1 test.stat are so few that they do not mess up the estimation of the gaussian parameters
  x <- shave ( test.stat, alpha )
  pnorm.eval <- pnorm( test.stat, mean(x), sd(x))
  pvalues = ifelse ( pnorm.eval > 0.5, 2*(1-pnorm.eval), 2*pnorm.eval )
}
