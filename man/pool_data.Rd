\name{as_pool_gwas_data}
\alias{as_pool_gwas_data}
\alias{validate_pool_data}
\alias{pool_from_individuals}
\alias{compute_group_MAFs}
\alias{print.pool_gwas_data}
\title{Construct and validate explicitly identified pooled data}
\description{Aligns all inputs by unique identifiers and validates frequencies, diploid sizes and typed counts. Frequencies use one consistently coded allele across pools.}
\usage{
as_pool_gwas_data(frequencies = NULL, phenotypes, pool_sizes,
  counts = NULL, covariates = NULL, variants = NULL,
  missing = c("error", "drop_variants"), min_maf = 0, min_depth = 0,
  design = list(), provenance = list())
validate_pool_data(data)
pool_from_individuals(genotypes, phenotypes, clusters, covariates = NULL,
  retain_individuals = FALSE, ...)
compute_group_MAFs(df, groups)
\method{print}{pool_gwas_data}(x, ...)
}
\arguments{
\item{frequencies}{Numeric pools-by-variants matrix in [0,1], with row and column IDs. May be omitted if counts are supplied.}
\item{phenotypes}{Numeric matrix of pool means, or individual phenotypes for \code{pool_from_individuals}. Row and trait IDs are required.}
\item{pool_sizes}{Positive integer diploid individual counts, named by pool ID.}
\item{counts}{List with identified integer matrices \code{alt}, \code{depth}, and explicit \code{type}: \code{"reads"} or \code{"chromosomes"}. Chromosome depth must equal twice pool size. Counts and frequencies must agree.}
\item{covariates}{Optional numeric matrix with matching row IDs, aggregated to pool means in the individual interface. No intercept column should be supplied.}
\item{variants}{Optional metadata data frame with a unique \code{variant_id} column.}
\item{missing}{Reject missing measurements, or drop every variant with a missing measurement or zero depth in any pool. No denominator-changing imputation occurs.}
\item{min_maf}{Minimum globally weighted minor allele frequency. Default zero retains monomorphic variants for null calibration.}
\item{min_depth}{Minimum count depth in every pool. Requires counts when positive.}
\item{design, provenance}{Lists documenting the pooling design and data source.}
\item{data, x}{A \code{pool_gwas_data} object.}
\item{genotypes, df}{Individual-by-variant matrix or numeric data frame of observed diploid 0/1/2 dosages with IDs.}
\item{clusters}{Named individual-to-pool vector, or a \code{pool_design} object. Unused factor levels are ignored.}
\item{groups}{Named individual-to-pool vector.}
\item{retain_individuals}{Whether to retain individual inputs in the result. Default avoids unnecessary retention.}
\item{...}{Additional constructor arguments in \code{pool_from_individuals}; ignored by print.}
}
\details{The constructor requires at least two occupied pools. Constant traits and variants are retained and flagged in QC. Metadata for excluded variants records explicit reasons. The historical name \code{compute_group_MAFs} is retained, but returns the coded-allele frequency, not an independently flipped minor allele frequency per pool. Samples with missing genotypes or phenotypes are rejected; resolve missingness before pooling with an explicit documented policy. This data object is not evidence of a valid association model.}
\value{The constructor and individual aggregator return \code{pool_gwas_data} with frequencies, phenotype means, pool sizes, typed counts, covariates, metadata and QC. \code{validate_pool_data} returns a list with \code{valid}, \code{data} and \code{qc}. \code{compute_group_MAFs} returns an identified frequency matrix. Print returns its input invisibly.}
\examples{
F <- matrix(c(.1,.3,.7,.9,.2,.4,.5,.8), 4, 2,
  dimnames = list(paste0("pool", 1:4), c("snp1", "snp2")))
Y <- matrix(c(-2,-1,1,2), 4, 1,
  dimnames = list(rownames(F), "trait"))
dat <- as_pool_gwas_data(F, Y, setNames(rep(20,4), rownames(F)))
validate_pool_data(dat)$valid
}
