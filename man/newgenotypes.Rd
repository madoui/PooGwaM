% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/newgenotypes.R
\name{newgenotypes}
\alias{newgenotypes}
\title{Regenerate genotype from allele freuencies}
\usage{
newgenotypes(pool_size, clusters, sim, allele_freq)
}
\arguments{
\item{pool_size}{the desired pool size}

\item{clusters}{the pools of individuals}

\item{sim}{a PhenotypeSimulator Object}

\item{allele_freq}{allele frequencies}
}
\value{
a dataframe of new genotypes
}
\description{
The goall is to simulate new genotypes from the observed allele frequencies in the pools
}
\examples{
sim = sim2snp(100,100,10)
pools = phen2pool(sim, 3, 3)
GenotypesCluster = trugen2clust (pools, sim)
alleleFreq = pool2freq (pools, GenotypesCluster, sim)
new_genotypes = newgenotypes (GenotypesCluster, pools, sim, alleleFreq)
}
