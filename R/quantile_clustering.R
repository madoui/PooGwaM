#' Cluster individuals based on their phenotypes into equal-sized pool for pool-seq
#'
#' The goal is to cluster the individuals into pool based on their phenotypes,
#'
#' @param data data frame of phenotypes
#' @param n_quantiles numeric vector for the number of quantiles for each phenotype
#' @return list
#'
#' @examples
#' clusters <- quantile_clustering(y<-data.frame(matrix(rnorm(500),500,1)))
#' plot(data.frame(1:500,y),col=clusters,main="1D phenotypes",xlab="Index",ylab="Trait 1")
#'
#' @export
#'
quantile_clustering <- function(data, n_quantiles = 4) {
  if (!is.data.frame(data)) {
    stop("Input data must be a data frame.")
  }
  if (any(!sapply(data, is.numeric))) {
    stop("All columns in the data frame must be numeric.")
  }
  # Compute quantile boundaries for each column
  quantile_boundaries <- lapply(data, function(column) {
    quantile(column, probs = seq(0, 1, length.out = n_quantiles + 1), na.rm = TRUE)
  })
  # Function to assign a vector to a cluster based on its quantiles
  assign_cluster <- function(vector) {
    cluster <- sapply(seq_along(vector), function(i) {
      # Find the quantile index of the current variable
      quantile_index <- sum(vector[i] >= quantile_boundaries[[i]]) - 1
      # Convert the quantile index to a character and concatenate it
      paste0("Q", quantile_index)
    })
    paste(cluster, collapse = "-")
  }

  # Assign each row to a cluster
  clusters <- apply(data, 1, assign_cluster)
  return(clusters=as.integer(factor(clusters)))
}









