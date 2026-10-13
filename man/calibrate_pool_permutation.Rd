\name{calibrate_pool_permutation}
\alias{calibrate_pool_permutation}
\title{Individual permutation calibration under global independence}
\usage{
calibrate_pool_permutation(genotypes, phenotypes, pool_builder,
  B = 999, seed = 1, strata = NULL, fit_args = list())
}
\arguments{
\item{genotypes}{Finite individual-by-variant matrix of observed diploid 0/1/2 dosages, with unique row and column identifiers.}
\item{phenotypes}{Finite individual-by-trait numeric matrix with unique identifiers. Individuals are aligned by ID.}
\item{pool_builder}{Function taking the complete phenotype matrix and returning named individual pool labels or a \code{pool_design}. Its design must depend on phenotypes only.}
\item{B}{Positive integer number of random permutations.}
\item{seed}{Nonnegative integer seed. The caller's RNG kind and state are restored, including when an error occurs.}
\item{strata}{Optional vector of exchangeability blocks named by every individual ID. Permutations remain within blocks. Blocks require a scientific exchangeability justification.}
\item{fit_args}{Uniquely named list of fixed arguments for \code{fit_pool_ridge}, excluding \code{data}. No tuning or model selection is performed.}
}
\description{Permutes complete genotype rows relative to fixed multivariate phenotype vectors and reconstructs the pool frequencies at each iteration. Returns Monte Carlo calibration of absolute ridge effects and standardized multivariate effect norms.}
\details{
The phenotype-only pooling rule is evaluated once with a controlled RNG. Conditional on the phenotypes and realized design, labels remain fixed while genotype rows are permuted and pool allele frequencies are rebuilt. Thus all phenotype correlations and all within-row genomic linkage patterns are retained. Permuting whole phenotype vectors and rebuilding their design is an equivalent representation when the same phenotype-only rule and its randomness are coupled appropriately.

The statistic for a trait is the absolute original-unit coefficient. The omnibus statistic is the Euclidean norm of coefficients divided by weighted pool-phenotype standard deviations. It is a score, not a chi-squared statistic. Null exceedances include ties, with a relative tolerance of 64 machine epsilons times the observed score, preserving invariance to trait measurement units. Every Monte Carlo p-value is \code{(1 + exceedances)/(B + 1)}. Minimum attainable p-value is \code{1/(B+1)}.

Validity assumes exchangeability of complete genotype rows with respect to phenotype vectors, globally or within the supplied blocks. This tests the global null of independence. Per-variant scores do not automatically test individual causal nulls when other loci are associated or in linkage disequilibrium. The maximum omnibus statistic supplies \code{p_maxT_omnibus} under the same global null; no strong familywise error or partial-null FDR guarantee is claimed.

Families, ancestry, batches and continuous covariates cannot be corrected by declaring arbitrary blocks. Covariate adjustment and data-dependent tuning are unsupported here. Settings and pooling rules must be specified independently of genotype-phenotype association. This exact-dosage procedure does not simulate sequencing reads; use of read-based observations would require incorporating their generation in the entire null experiment.
}
\value{A \code{poogwam_calibration} list containing the observed \code{fit}, \code{p_raw} (variant-by-trait matrix), named \code{p_omnibus} and \code{p_maxT_omnibus} vectors, \code{B}, \code{resolution}, RNG seed, assumptions, limitations, and streaming \code{null_summary}. The inference scope is always \code{global_null_only}.}
\seealso{\code{\link{association_table}}, \code{\link{calibrate_pool_counts}}, \code{\link{fit_pool_ridge}}}
\examples{
ids <- paste0("i", 1:12)
g <- matrix(rep(c(0, 1, 2), 8), 12, 2,
            dimnames = list(ids, c("s1", "s2")))
y <- cbind(height = seq_len(12), mass = rep(1:6, 2))
rownames(y) <- ids
cal <- calibrate_pool_permutation(g, y,
  pool_builder = function(y) design_pools(y, k = 3, method = "balanced"),
  B = 9, seed = 7)
association_table(cal$fit, cal, adjustment = "holm")
}
