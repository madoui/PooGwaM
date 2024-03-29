% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PhenoSim.R
\name{compute_group_MAFs}
\alias{compute_group_MAFs}
\title{compute_group_MAFs}
\usage{
compute_group_MAFs(df, groups)
}
\arguments{
\item{df}{data frame of phenotypes}

\item{groups}{numeric vector}
}
\value{
data frame
}
\description{
Compute per group major allele frequencies
}
\examples{
sim = PhenoSim (1000, 100, 10, 0.6, 5)
Y<-data.frame(sim$phenotypes)
clusters <- quantile_clustering(Y, 5)

}
