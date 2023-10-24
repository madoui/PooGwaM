#' singValRidge
#'
#' apply singular value decomposition for ridge analysis
#'
#' @param phenotypes data frame of phenotypes
#' @param Freq numeric vector of allele frequencies
#' @param clusters vector of cluster ID
#' @param K number of clusters
#'
#' @return numeric vector
#'
#' @examples
#' sim <- PhenoSim (1000, 100, 10, 0.6, 5)
#' clusters <- quantile_clustering (data.frame(sim$phenotypes))
#' Freq <- compute_group_MAFs(sim$genotypes, as.factor(clusters))
#' criterion <- singValRidge(sim$phenotypes, Freq, clusters, K = 4)
#'
#' @importFrom stats optim pbeta sd
#' @importFrom RSpectra svds
#'
#' @export

singValRidge <- function(phenotypes, Freq, clusters, K){
  y<-scale(phenotypes,center=TRUE)

  n<-nrow(phenotypes)

  x = nnet::class.ind(clusters)%*%as.matrix(Freq)

  alpha<-((2*colMeans(x)*(1-colMeans(x))))

  UDV<-svds(x, K )

  alpha.mean<-mean(alpha)

  S<-UDV$d/(UDV$d^2+n*alpha.mean)

  beta.ridge.svd<-UDV$v%*%(S*t(UDV$u))%*%y
  criterion.ridge.svd<- mean((x%*%beta.ridge.svd - y)^2)+  sum(alpha.mean*beta.ridge.svd^2)

  return(list("beta" = beta.ridge.svd, "crit" = criterion.ridge.svd) )
}

