#' Simulate phenotypes and genotypes using PhenotypeSimulator
#'
#' PhenotypeSimulator wrapper to generates p-dimensional phenotypes
#' Total number of individuals \code{n}, default=2000
#' Total number of SNPs \code{SNP}, default=100
#' Number of phenotypes \code{p}, default=5
#' SNP heritability \code{h2s}, default=50%
#' Number of causal SNPs \code{n_causal}, default=10
#'
#' @param n Numeric vector
#' @param SNP Numeric vector.
#' @param h2 Numeric vector
#' @param n_causal Numeric vector
#' @param p Numeric vector

#' @return a list
#'
#' @importFrom stats na.omit
#'
#' @examples
#' sim = PhenoSim (1000, 100, 10, 0.6, 5)
#'
#' @export
PhenoSim<-function(n = 2000, # number of individuals
                   SNP = 100, # number of SNPs
                   n_causal = 10, # number of causal SNPs
                   h2 = 0.5, # heritability
                   p = 5 #number of traits
                   ){

  # run PhenotypeSimulator
  simulation = PhenotypeSimulator::runSimulation(N = n, # individuals
                             P = p, # traits
                             tNrSNP = SNP, # total SNPs
                             cNrSNP = n_causal, # causal SNPs
                             noiseVar =0.01,
                             h2s=h2,
                             rho=0,
                             delta=0,
                             SNPfrequencies=c(0.5),
                             pIndependentGenetic = 1, seed=16 )

  # store simulation
  phenotypes = simulation[["phenoComponentsFinal"]][["Y"]]
  genotypes = simulation[["rawComponents"]][["genotypes"]][["genotypes"]]
  causal.index<-as.vector(na.omit(as.integer(
  unlist(strsplit(colnames(
  simulation$phenoComponentsIntermediate$genFixed$cov),"SNP_")))))
  causal.true<-rep(0,SNP)
  causal.true[causal.index]<-1

  return(list(h2.estimation=h2,
              causal.true=causal.true,
              phenotypes=data.frame(phenotypes),
              genotypes=data.frame(genotypes)))
}

#' compute_group_MAFs
#'
#' Compute per group major allele frequencies
#'
#' @param df data frame of phenotypes
#' @param groups numeric vector
#' @return data frame
#' @importFrom dplyr %>% group_by select summarise_all
#' @examples
#' sim = PhenoSim (1000, 100, 10, 0.6, 5)
#' Y<-data.frame(sim$phenotypes)
#' clusters <- quantile_clustering(Y, 5)
#'
#' @export
#'
compute_group_MAFs <- function(df, groups) {
  cbind(groups,df) %>%
    group_by(groups) %>%
    summarise_all(mean) %>%
    select(-1)/2 # Careful about the mean ! We have to divide by 2 since there are 2n DNA strands
}



#' make_pooled_data_design_matrix
#'
#' create the design matrix of the pooled data
#'
#' @param genotypes data frame
#' @param phenotypes data frame
#' @param k integer
#' @return list
#' @importFrom stats quantile
#' @examples
#' sim = PhenoSim (1000, 100, 10, 0.6, 5)
#' Y<-data.frame(sim$phenotypes)
#' X<-sim$genotypes
#' X<-make_pooled_data_design_matrix(X, Y, 7)
#'
#' @export
#'
make_pooled_data_design_matrix<-function(genotypes,phenotypes, k = 2){
  clusters <- quantile_clustering(phenotypes, k)
  Freq= compute_group_MAFs(data.frame(genotypes),as.factor(clusters))
  C<-nnet::class.ind(clusters)# Change coding of classing from integer to binary
  return(X = (C%*%as.matrix(Freq)))
}








