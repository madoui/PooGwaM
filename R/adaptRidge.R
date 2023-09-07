#' adaptRidge
#'
#' adaptative Ridge method with descending gradient
#'
#' @param clusters pool ID
#' @param Freq data frame of allele frequencies
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
#' res = adaptRidge( Freq, sim$phenotypes[,1], clusters )
#' plot ( res$criterion )
#' beta = res$beta
#' pvalues<- pvalues.from.test.stat(beta)
#'
#' @export
#'
#'
adaptRidge <- function( Freq = Freq,
                        phenotypes = phenotypes,
                        clusters = clusters,
                        eta = 1,
                        max_iter = 1000,
                        tol = 1e-8,
                        seed = 1,
                        nb.batch = 1,
                        trace = FALSE) {

  y <- scale(phenotypes, center = TRUE)
  Z<-nnet::class.ind(clusters)
  X = Z%*%as.matrix(Freq)
  alpha<-((2*colMeans(X)*(1-colMeans(X))))

  # Set random seed for reproducibility
  if(!is.null(seed)){
    set.seed(seed)
  }

  # Initialize coefficients and intercept
  n <- nrow(X)
  p <- ncol(X)
  batch.size<-round(n/nb.batch)
  meany<-mean(y)
  y<-matrix(y-meany,n,1)
  beta <- t(X)%*%y/(1.5*n)
  grad<-matrix(1,p,1)

  # Initialize variables for adaptive learning rate
  k<-1
  criterion<-rep( 0, max_iter*nb.batch )
  # SGD loop
  for (iter in 1:max_iter) {
    # Randomly shuffle the data
    idx <- sample(n, n)
    X <- X[idx, ]
    y <- y[idx]
    # SGD update
    for (b in 1:nb.batch) {
      criterion[k] <- mean((X%*%beta - y)^2)+ sum(alpha*beta^2)
      batchi<-sample(n,batch.size,replace=FALSE)
      xi <- matrix(X[batchi, ],batch.size,p) # get a matrix even if batch.size==1
      yi <- cbind(y[batchi])
      y_pred <- xi %*% beta
      grad <- (-2/batch.size * t(xi) %*% (yi - y_pred) + 2*alpha * beta )
      beta.canditate <- beta - eta * grad
      while ( (mean((X%*%beta.canditate - y)^2)+ sum(alpha*beta.canditate^2)) < criterion[k]){
        beta<-beta.canditate
        grad <- (-2/batch.size * t(xi) %*% (yi - y_pred) + 2*alpha * beta )
        eta<-eta*1.1
        beta.canditate <- beta - eta* grad
      }
      eta<-eta*0.9
      k<-k+1
    }

    if (trace==TRUE) print(paste("Iteration: ",k," gradient norm : ", (sqrt(sum(grad^2)))))
    # Check for convergence
    if (sqrt(mean(grad^2)) < tol) {
      break
    }
  }

  # Return coefficients
  return( list( beta = beta , criterion = criterion ) )
}


