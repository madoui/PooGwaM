\name{as_pool_gwas_data}
\alias{as_pool_gwas_data}
\alias{validate_pool_data}
\alias{pool_from_individuals}
\alias{compute_group_MAFs}
\alias{print.pool_gwas_data}
\title{Construct and validate explicitly identified pooled data}
\description{Aligns all inputs by unique identifiers and validates frequencies, diploid sizes and typed counts. Frequencies use one consistently coded allele across pools. Fixed cofactors can be supplied as a numeric matrix or encoded from a formula and identified metadata.}
\usage{
as_pool_gwas_data(frequencies = NULL, phenotypes, pool_sizes,
  counts = NULL, covariates = NULL, variants = NULL,
  missing = c("error", "drop_variants"), min_maf = 0, min_depth = 0,
  design = list(), provenance = list(), metadata = NULL)
validate_pool_data(data)
pool_from_individuals(genotypes, phenotypes, clusters, covariates = NULL,
  retain_individuals = FALSE, ..., metadata = NULL)
compute_group_MAFs(df, groups)
\method{print}{pool_gwas_data}(x, ...)
}
\arguments{
\item{frequencies}{Numeric pools-by-variants matrix in [0,1], with row and column IDs. May be omitted if counts are supplied.}
\item{phenotypes}{Numeric matrix of pool means, or individual phenotypes for \code{pool_from_individuals}. Row and trait IDs are required.}
\item{pool_sizes}{Positive integer diploid individual counts, named by pool ID.}
\item{counts}{List with identified integer matrices \code{alt}, \code{depth}, and explicit \code{type}: \code{"reads"} or \code{"chromosomes"}. Chromosome depth must equal twice pool size. Counts and frequencies must agree.}
\item{covariates}{Optional numeric matrix with matching row IDs, or a one-sided
  formula evaluated in \code{metadata}, for example
  \code{~ sex + age + batch + site + treatment}. Numeric matrices are treated
  as already encoded and must not contain an intercept column. A formula must
  include its usual intercept; the resulting intercept column is removed
  because the solver adds it separately. \code{~ 1} represents an intercept-only
  covariate model. In the individual interface, encoding precedes averaging
  within each pool. The solver treats these covariates as unpenalized fixed
  effects.}
\item{metadata}{Optional data frame used with a formula in \code{covariates}.
  Explicit unique row names must identify exactly the pools for
  \code{as_pool_gwas_data}, or exactly the individuals for
  \code{pool_from_individuals}. Rows are aligned by ID. Required columns can
  be numeric, logical, factors or characters; character variables are encoded
  as factors. Missing values in the variables used by the formula are rejected,
  without dropping rows. Logical variables retain levels \code{FALSE} and
  \code{TRUE}; a factor declared with fewer than two levels is rejected.
  Supply this argument by name, especially after the
  dots in \code{pool_from_individuals}.}
\item{variants}{Optional metadata data frame with a unique \code{variant_id} column.}
\item{missing}{Reject missing measurements, or drop every variant with a missing measurement or zero depth in any pool. No denominator-changing imputation occurs.}
\item{min_maf}{Minimum globally weighted minor allele frequency. Default zero retains monomorphic variants for null calibration.}
\item{min_depth}{Minimum count depth in every pool. Requires counts when positive.}
\item{design, provenance}{Lists documenting the pooling design and data source.}
\item{data, x}{A \code{pool_gwas_data} object.}
\item{genotypes, df}{Individual-by-variant matrix or numeric data frame of observed diploid 0/1/2 dosages with IDs.}
\item{clusters}{Named individual-to-pool vector, or a \code{pool_design} object. Unused factor levels are ignored.}
\item{groups}{Named individual-to-pool vector.}
\item{retain_individuals}{Whether to retain individual genotypes and phenotypes
  in the result's \code{individuals} component. Formula metadata and cluster
  memberships are retained for encoding validation and prediction even when
  this argument is \code{FALSE}.}
\item{...}{Additional constructor arguments in \code{pool_from_individuals}; ignored by print.}
}
\details{
The constructor requires at least two occupied pools. Constant traits and
variants are retained and flagged in QC. Metadata for excluded variants records
explicit reasons. The historical name \code{compute_group_MAFs} is retained,
but returns the coded-allele frequency, not an independently flipped minor
allele frequency per pool. Samples with missing genotypes or phenotypes are
rejected; resolve missingness before pooling with an explicit documented policy.
This data object is not evidence of a valid association model.

With an individual-level formula, factor indicators, supported transformations
and interactions are constructed for each individual before pooling. An
indicator's pool mean is the proportion of individuals in that category.
For \code{I(age^2)} the covariate is the mean of individual squared ages, not
the square of mean age. Similarly, a pooled interaction is the mean of the
individual interaction. In contrast, a formula in \code{as_pool_gwas_data}
is evaluated directly on pool-level metadata; individual transformations
cannot generally be recovered from pooled means alone.

Factor levels, contrasts and the training encoding are stored for prediction.
Set factor levels explicitly before construction to choose references under
treatment contrasts. Prediction reuses the training encoding rather than
estimating a new coding scheme. Unknown categories and unsupported formula
functions are errors. List variables explicitly: \code{.}, a response on the
left-hand side and removal of the intercept are not supported.

Allowed transformations are arithmetic with \code{I}, \code{log},
\code{log1p}, \code{sqrt}, \code{exp}, \code{abs}, \code{scale} and
\code{poly}, alongside formula interactions. Direct \code{scale(age)} and
\code{poly(age, 2)} terms retain their learned centering, scaling and polynomial
basis for prediction. Stateful transformations nested inside other transforms,
such as \code{I(scale(age))} or \code{poly(scale(age), 2)}, are rejected instead
of recomputing an inconsistent prediction encoding. Custom functions and
variables looked up outside \code{metadata} are not supported.

The full pooled matrix formed by the intercept and encoded covariates must
have full column rank and fewer columns than occupied pools when fitting.
Pooling can create collinearity even when the individual design has full rank.
Neither construction nor fitting silently removes redundant fixed-effect terms.
Inspect \code{fit$diagnostics$covariate_design} after a successful fit.

Formula-based covariates are supported for effect estimation and prediction.
The permutation, count-bootstrap and knockoff APIs do not provide adjusted
inference for this extension. Their existing restrictions still apply.
}
\value{
The constructor and individual aggregator return \code{pool_gwas_data} with
frequencies, phenotype means, pool sizes, typed counts, numeric pooled
covariates, variant metadata and QC. Formula inputs additionally retain
\code{covariate_spec}, \code{covariate_metadata} and, for individual inputs,
\code{covariate_clusters}, to reproduce the encoding at prediction.
\code{validate_pool_data} returns a list with \code{valid}, \code{data} and
\code{qc}. \code{compute_group_MAFs} returns an identified frequency matrix.
Print returns its input invisibly.
}
\examples{
F <- matrix(c(.1,.3,.7,.9,.2,.4,.5,.8), 4, 2,
  dimnames = list(paste0("pool", 1:4), c("snp1", "snp2")))
Y <- matrix(c(-2,-1,1,2), 4, 1,
  dimnames = list(rownames(F), "trait"))
dat <- as_pool_gwas_data(F, Y, setNames(rep(20,4), rownames(F)))
validate_pool_data(dat)$valid

# Encode individual cofactors before taking pool means.
ids <- paste0("person", 1:12)
G <- matrix(rep(c(0, 1, 2), 8), 12, 2,
  dimnames = list(ids, c("snp1", "snp2")))
individual_y <- matrix(150 + 1:12, 12, 1,
  dimnames = list(ids, "height"))
groups <- setNames(rep(paste0("pool", 1:4), each = 3), ids)
metadata <- data.frame(
  sex = factor(c("F", "F", "F", "F", "F", "M",
                 "F", "M", "M", "M", "F", "F"), levels = c("F", "M")),
  age = c(21, 22, 23, 30, 32, 34, 40, 42, 44, 50, 52, 54),
  row.names = ids)
adjusted_data <- pool_from_individuals(G, individual_y, groups,
  covariates = ~ sex + age, metadata = metadata)
adjusted_data$covariates
}
