#' GWAlpha
#'
#' apply GWAlpha method from Fournier-Level et 2017
#'
#' @param phenotypes data frame of phenotypes
#' @param clusters numeric vector of the clusters
#' @param Freq numeric vector of allele frequencies
#' @param K numeric vector of number of clusters
#' @param penalized boolean
#'
#' @return numeric vector
#'
#' @examples
#' sim = PhenoSim (1000, 100, 10, 0.6, 5)
#' K = 7
#' clusters = quantile_clustering (data.frame(sim$phenotypes), K)
#' Freq= compute_group_MAFs(sim$genotypes,as.factor(clusters))
#' test.stat<-GWalpha(sim$phenotypes[,1], clusters, Freq, K)
#'
#' @importFrom stats optim pbeta sd
#'
#' @export

# Problem with the order of the frequencies
GWalpha<-function(phenotypes = phenotypes,
                  #Yprime = Yprime,
                  clusters = clusters,
                  Freq = Freq,
                  K = K,
                  #Pheno ? 1 2 ...
                  penalized = TRUE
                  ){

  Yprime = seq(0, 1, 1/K)

  LS<-function(parameters, Yprime, freqj){
    sum((freqj - pbeta(Yprime[-1],parameters[1],parameters[2]) + pbeta(Yprime[-(K+1)],parameters[1],parameters[2]))^2)+
      sum((1-freqj - pbeta(Yprime[-1],parameters[3],parameters[4]) + pbeta(Yprime[-(K+1)],parameters[3],parameters[4]))^2)
  }

  miny<-min(phenotypes);
  maxy<-max(phenotypes)
  Alpha<-rep(1,dim(Freq)[2])

  #Y = phenotypes[,1]

  for (j in 1:dim(Freq)[2]){
    freqj<-Freq[,j]
    solution <- optim( c(1,1,1,1), LS, gr=NULL, Yprime, freqj, control=list(abstol=1e-8) )$par
    muA_hat = miny+(maxy-miny)*solution[1]/(solution[1]+solution[2])
    muB_hat = miny+(maxy-miny)*solution[3]/(solution[3]+solution[4])
    Alpha[j] = ifelse(penalized,2*sqrt(freqj*(1-freqj))*(muA_hat-muB_hat)/sd(phenotypes),(muA_hat-muB_hat)/sd(phenotypes))
  }
  return(Alpha)
}


