\name{design_pools}
\alias{design_pools}
\alias{quantile_clustering}
\alias{print.pool_design}
\title{Deterministic phenotype-defined pooling}
\description{Constructs a Cartesian grid of marginal phenotype quantiles or a recursively balanced partition with a fixed total pool count.}
\usage{
design_pools(phenotypes, k = 4, method = c("quantile", "balanced"), min_size = 2L)
quantile_clustering(data, n_quantiles = 4, min_size = 1L)
\method{print}{pool_design}(x, ...)
}
\arguments{
\item{phenotypes, data}{Finite numeric individual-by-trait matrix or data frame with unique row and trait IDs.}
\item{k}{For quantiles, one count or a count per trait. For balanced partitioning, one total number of pools.}
\item{method}{Quantile grid or deterministic recursive balanced partition.}
\item{min_size}{Minimum individual count per occupied pool. Infeasible designs are rejected.}
\item{n_quantiles}{One quantile count or one per trait for the compatibility helper.}
\item{x}{A pooling design.}
\item{...}{Ignored print arguments.}
}
\details{The quantile grid uses observed order statistics. Equal values are never split along a cut. Thus tied traits may yield fewer marginal bins. Grid pools need not have equal sizes; a violated minimum size raises an error. The balanced method standardizes the traits, considers near-half splits at distinct values in each trait and chooses the split with largest between-child sum-of-squares reduction. Row ID ordering makes both designs invariant to input row order. Neither procedure evaluates genotypes. Balanced multivariate pools may overlap along a single trait and must not automatically be exported to GWAlpha.}
\value{A \code{pool_design} list containing named memberships, observed bounds, sizes, settings and, for quantiles, marginal bin indices and thresholds. The compatibility function returns only named factor memberships. Print returns its input invisibly.}
\examples{
Y <- cbind(height = 1:40, mass = sin(1:40))
rownames(Y) <- paste0("id",1:40)
design <- design_pools(Y, k = 4, method = "balanced", min_size = 5)
print(design)
}
