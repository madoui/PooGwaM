#' GWAlpha
#'
#' apply the GWAlpha algorithm from Fournier-Level et al. 2017
#'
#' @param phenotypes data frame of phenotypes
#' @param d phenotypes index to be analyzed
#' @param clusters numeric vector of the clusters
#' @param Freq numeric vector of allele frequencies
#' @param penalized boolean
#'
#' @return numeric vector
#'
#' @examples
#' sim = PhenoSim (1000, 100, 10, 0.6, 5)
#' clusters = quantile_clustering (data.frame(sim$phenotypes))
#' Freq= compute_group_MAFs(sim$genotypes,as.factor(clusters))
#' test.stat<-GWalpha(sim$phenotypes,d=1, clusters, Freq)
#'
#' @importFrom stats optim pbeta sd
#'
#' @export

GWalpha<-function(phenotypes = phenotypes,
                  d=1,
                  #Yprime = Yprime,
                  clusters = clusters,
                  Freq = Freq,
                  penalized = TRUE
                  ){
  # Attributing each ind to new univariate clusters based on clustering phenotype j
  new.clusters<-factor(sapply(clusters,function(x) substring(x,(d-1)*5+1,(d-1)*5+4)))
  K<-length(levels(new.clusters)) # find the new reduced number of clusters along the phenotype j
  Yprime = seq(0, 1, 1/K)

  # the correspondance between the two clustering is given by the matrix
  correspondance<-dplyr::distinct(data.frame(clusters,new.clusters))
  clusters.of.old.clusters<-correspondance[,2]


  # Attribute new allele frequencies corresponding to the new clusters to each individual
  clusters.size<-table(clusters)
  new.clusters.size<-table(new.clusters)

  recompute.freqj<-function(freqj){
     tapply(freqj*clusters.size ,clusters.of.old.clusters,sum)/ tapply(clusters.size ,clusters.of.old.clusters,sum)}
  Freq<-apply(Freq,2,recompute.freqj) # Passing from Freq corresponding to old clusters to Fred corresponding to new metacluster



  LS<-function(parameters, Yprime, freqj){
    sum((freqj - pbeta(Yprime[-1],parameters[1],parameters[2]) + pbeta(Yprime[-(K+1)],parameters[1],parameters[2]))^2)+
      sum((1-freqj - pbeta(Yprime[-1],parameters[3],parameters[4]) + pbeta(Yprime[-(K+1)],parameters[3],parameters[4]))^2)
  }

  miny<-min(phenotypes);
  maxy<-max(phenotypes)
  Alpha<-rep(1,dim(Freq)[2])

  for (j in 1:dim(Freq)[2]){
    freqj<-Freq[,j]
    solution <- optim( c(1,1,1,1), LS, gr=NULL, Yprime, freqj, control=list(abstol=1e-8) )$par
    muA_hat = miny+(maxy-miny)*solution[1]/(solution[1]+solution[2])
    muB_hat = miny+(maxy-miny)*solution[3]/(solution[3]+solution[4])
    Alpha[j] = ifelse(penalized,2*sqrt(freqj*(1-freqj))*(muA_hat-muB_hat)/sd(phenotypes),(muA_hat-muB_hat)/sd(phenotypes))
  }
  return(Alpha)
}


