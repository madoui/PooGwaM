\name{calibrate_pool_permutation}
\alias{calibrate_pool_permutation}
\title{Individual permutation calibration under global independence}
\usage{
calibrate_pool_permutation(genotypes, phenotypes, pool_builder,
  B = 999, seed = 1, strata = NULL, fit_args = list(),
  maxT = c("raw", "studentized"), B_pilot = 199,
  scale_floor = 1e-8, rng_kind = NULL)
}
\arguments{
\item{genotypes}{Finite individual-by-variant matrix of observed diploid 0/1/2 dosages, with unique row and column identifiers.}
\item{phenotypes}{Finite individual-by-trait numeric matrix with unique identifiers. Individuals are aligned by ID.}
\item{pool_builder}{Function taking the complete phenotype matrix and returning named individual pool labels or a \code{pool_design}. Its design must depend on phenotypes only.}
\item{B}{Positive integer number of random permutations.}
\item{seed}{Nonnegative integer seed. The caller's RNG kind and state are restored, including when an error occurs.}
\item{strata}{Optional vector of exchangeability blocks named by every individual ID. Permutations remain within blocks. Blocks require a scientific exchangeability justification.}
\item{fit_args}{Uniquely named list of fixed arguments for \code{fit_pool_ridge}, excluding \code{data}. No tuning or model selection is performed.}
\item{maxT}{Maximum omnibus statistic: \code{"raw"} preserves the V2 procedure; \code{"studentized"} normalizes each variant's omnibus score using a separate permutation pilot. Raw remains the default pending confirmatory evaluation of the studentized option.}
\item{B_pilot}{Integer number of pilot permutations, at least two in studentized mode. The default 199 is a provisional computational setting, not a validated power recommendation. Ignored in raw mode.}
\item{scale_floor}{Relative lower bound on pilot standard deviations, in \eqn{(0,1]}. The effective floor is this value times the largest pilot mean or standard deviation, bounded below by \code{.Machine$double.xmin}. Ignored in raw mode.}
\item{rng_kind}{Optional character vector of the three values returned by \code{RNGkind()}. By default the caller's RNG algorithms are inherited. Pass a previous result's \code{rng$kind}, together with its seed, inputs and settings, to replay that configuration. The caller's kind and state are restored on return or error. Box-Muller normal generation and user-supplied RNGs are rejected before changing the RNG: their hidden/external state cannot be restored from \code{.Random.seed}. Select a built-in generator with \code{normal.kind="Inversion"} before calling if necessary.}
}
\description{Permutes complete genotype rows relative to fixed multivariate phenotype vectors and reconstructs the pool frequencies at each iteration. Returns Monte Carlo calibration of absolute ridge effects and standardized multivariate effect norms.}
\details{
The phenotype-only pooling rule is evaluated once with a controlled RNG. Conditional on the phenotypes and realized design, labels remain fixed while genotype rows are permuted and pool allele frequencies are rebuilt. Thus all phenotype correlations and all within-row genomic linkage patterns are retained. Permuting whole phenotype vectors and rebuilding their design is an equivalent representation when the same phenotype-only rule and its randomness are coupled appropriately.

The statistic for a trait is the absolute original-unit coefficient. The omnibus statistic is the Euclidean norm of coefficients divided by weighted pool-phenotype standard deviations. It is a score, not a chi-squared statistic. Null exceedances include ties, with a relative tolerance of 64 machine epsilons times the observed score, preserving invariance to trait measurement units. Every Monte Carlo p-value is \code{(1 + exceedances)/(B + 1)}. Minimum attainable p-value is \code{1/(B+1)}.

In studentized mode the pilot estimates a mean and standard deviation for each variant's omnibus score. These define a fixed increasing transformation \eqn{t_j(s)=(s-\widehat\mu_j)/\mathrm{scale}_j}. The main series compares each observed transformed score to the genome-wide maximum of transformed null scores. The transform is not an assumed normal or Student t distribution. Pilot and main series use successive, non-overlapping segments of the inherited pseudorandom generator, sampling uniformly with replacement from the permitted group, including the identity; duplicates are not removed. The pilot is never replayed as the main series and does not enter its exceedance counts. There are \code{B_pilot + B} null fits and one observed fit.

Conditional on the permutation orbit, independent uniform pilot draws are independent of the observed data and main draws. Conditioning on that pilot fixes a common transformation for their rank comparison. This argument requires the stated exchangeability and fixed settings; it does not extend the inference scope beyond the global null. A two-pass normalization fitted only to the same main null replicates would not supply this argument.

The scale is the maximum of the pilot standard deviation and its positive floor. A zero pilot standard deviation does not establish a constant null statistic. If every pilot score is zero, the procedure warns and uses a fixed unit scale for every variant, reporting \code{status="all_zero_pilot_unit_scale"}; it does not force p-values to one. Overflow or a failed fit stops the calculation rather than dropping a draw. Centered scores may be negative. Their tie tolerance is 64 machine epsilons times the maximum of the absolute centered score, observed score divided by scale, and pilot mean divided by scale. Inspect the reported scales and sensitivity to \code{B_pilot} and \code{scale_floor}.

Validity assumes exchangeability of complete genotype rows with respect to phenotype vectors, globally or within the supplied blocks. This tests the global null of independence. Per-variant scores do not automatically test individual causal nulls when other loci are associated or in linkage disequilibrium. The maximum omnibus statistic supplies \code{p_maxT_omnibus} under the same global null; no strong familywise error or partial-null FDR guarantee is claimed.

Families, ancestry, batches and continuous covariates cannot be corrected by declaring arbitrary blocks. Covariate adjustment and data-dependent tuning are unsupported here. Settings and pooling rules must be specified independently of genotype-phenotype association. This exact-dosage procedure does not simulate sequencing reads; use of read-based observations would require incorporating their generation in the entire null experiment.
}
\value{A \code{poogwam_calibration} list containing the observed \code{fit}, \code{p_raw} (variant-by-trait matrix), named \code{p_omnibus} and \code{p_maxT_omnibus} vectors, \code{B}, \code{resolution}, \code{maxT_method}, RNG seed, assumptions, limitations, and streaming \code{null_summary}. The latter always summarizes the untransformed main-series scores, including their raw maxima. In studentized mode, \code{studentization} stores pilot size \code{B}, pilot moments, effective scales and floor, zero/floored variant IDs, status, observed transformed scores, transformed maxima, and comparison tolerances; otherwise it is \code{NULL}.

\code{maxT_diagnostics} contains the first variant attaining each raw/transformed maximum and the number of exact ties. First-winner IDs depend on column order when maxima are tied, so inspect tie counts when attributing maxima to rare or common variants.

\code{rng} records the complete RNG kind, seed, R/package versions, initial state and states at pooling, pilot and main-series boundaries. Package version is \code{NA} with \code{source_mode="sourced"} when functions are sourced outside a package namespace. \code{provenance} records MD5 fingerprints of R version-2 serializations of the inputs and fixed settings, the realized pool design, and the builder's formals/body (not its enclosing environment). Seed-based replay requires the same builder behavior, including its random-number consumption. Substituting a saved fixed design for a random builder alone does not restore the later RNG stream; the recorded \code{after_pooling} state identifies that boundary for explicit lower-level replay. These metadata do not constitute an automatic checkpoint/resume interface. Checksums identify inputs, not statistical validity; identical numerical replay also requires compatible R/package versions and numerical libraries. The inference scope is always \code{global_null_only}.}
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
studentized <- calibrate_pool_permutation(g, y,
  pool_builder = function(y) design_pools(y, k = 3, method = "balanced"),
  B = 9, seed = 7, maxT = "studentized", B_pilot = 9)
studentized$maxT_method
}
