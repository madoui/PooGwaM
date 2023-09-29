#' quantile_clustering
#'
#' Cluster individuals based on their phenotypes into equal-sized pool for pool-seq
#'
#' @param data data frame of phenotypes
#' @param n_quantiles numeric vector for the number of quantiles for each phenotype
#' @return factor
#' @importFrom stats quantile
#' @examples
#' par(mfrow=c(1,2))
#' clusters <- quantile_clustering(y<-data.frame(matrix(rnorm(500),500,1)))
#' plot(data.frame(1:500,y),col=clusters,main="1D phenotypes",xlab="Index",ylab="Trait 1")
#' clusters <- quantile_clustering(y<-data.frame(matrix(rnorm(1000),500,2)))
#' plot(y,col=clusters,main="2D phenotypes",xlab="Trait 1",ylab="Trait 2")
#'
#' @export
#'
quantile_clustering <- function(data, n_quantiles = 4) {
  cluster_labels<-unlist(data.frame(sapply(letters[1:2],function(x) paste(x,letters,sep=""))))

  if (!is.data.frame(data)) {
    stop("Input data must be a data frame.")
  }
  if (any(!sapply(data, is.numeric))) {
    stop("All columns in the data frame must be numeric.")
  }
  # Compute quantile boundaries for each column
  quantile_boundaries <- lapply(data, function(column) {
    quantile(column,
             probs = seq(0, 1, length.out = n_quantiles+1),
             na.rm = TRUE)
  })
  # Function to assign a vector to a cluster based on its quantiles
  assign_cluster <- function(vector) {
    cluster <- sapply(seq_along(vector), function(i) {
      # Find the quantile index of the current variable
      quantile_index <- sum(vector[i] >= quantile_boundaries[[i]][-length(quantile_boundaries[[i]])])
      # Convert the quantile index to a character and concatenate it
      paste0("Q.", cluster_labels[quantile_index])
    })
    paste(cluster, collapse = "-")
  }

  # Assign each row to a cluster
  clusters <- apply(data, 1, assign_cluster)
  return(clusters=factor(clusters))
}







