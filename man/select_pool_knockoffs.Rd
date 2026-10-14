\name{select_pool_knockoffs}
\alias{select_pool_knockoffs}
\title{Experimental Knockoff+ Selection After Phenotype-Defined Pooling}
\description{
Constructs symmetric importance contrasts from a supplied pair of individual
genotype matrices, pools both with the same phenotype-only plan, and applies
the knockoff+ selection rule. The procedure is experimental and does not
certify a supplied knockoff distribution or generate p-values.
}
\usage{
select_pool_knockoffs(genotypes, knockoffs, phenotypes, pool_builder,
  q = 0.1, method = c("ridge", "multitrait"), fit_args = list(),
  generator_provenance = NULL, allele_metadata = NULL, seed = 1)
}
\arguments{
\item{genotypes, knockoffs}{Numeric matrices of exact, observed diploid
  dosages 0, 1 or 2. Rows identify individuals and columns identify variants.
  Both matrices require unique nonempty row and column names, identical sets
  of IDs, and the same alternate-allele orientation for each paired column.
  Missing values, fractional dosages and pooled read frequencies are rejected.
  Inputs are aligned by their IDs, not their positions.}
\item{phenotypes}{Finite numeric matrix with named individuals matching the
  genotype matrices, and unique named trait columns. Multiple traits define
  a single vector-valued response for inference.}
\item{pool_builder}{Function taking only the phenotype matrix and returning
  a \code{pool_design} or a complete named individual-to-pool vector accepted
  by \code{\link{pool_from_individuals}}. The function is called once, with
  individuals in canonical ID order. Its implementation must not access
  genotypes or asymmetric tuning information through its enclosing environment.
  Any randomness must be independent of the original/knockoff labels.}
\item{q}{Target FDR for the knockoff+ rule, strictly greater than zero and
  at most one. The target has its usual theoretical meaning only under the
  Model-X and sign-flip assumptions described below.}
\item{method}{Joint solver for the augmented pooled matrix. \code{"ridge"}
  uses \code{\link{fit_pool_ridge}}; \code{"multitrait"} uses
  \code{\link{fit_pool_multitrait}}, including sparse-group penalties.}
\item{fit_args}{Uniquely named list of fixed arguments accepted by the selected
  solver, excluding \code{data}. Parameters must be chosen without asymmetric
  use of original and knockoff labels. A nonconverged solver causes an error
  instead of returning selections.}
\item{generator_provenance}{Optional named list recording the source of
  knockoffs, genotype model, parameters and estimation data. If omitted, a
  \code{knockoff_provenance} attribute on \code{knockoffs} is retained if
  available. Such declarations are metadata, not evidence certifying a model.}
\item{allele_metadata}{Optional list with exactly two named entries,
  \code{genotypes} and \code{knockoffs}. Each is a data frame with columns
  \code{variant_id}, \code{ref} and \code{alt}. Every variant must occur once;
  reference and alternate labels must be nonempty and distinct. The tables
  are aligned by ID and must agree exactly. No strand or allele flip is inferred.
  Without these tables, common allele orientation is an explicit unchecked
  assumption. Matching supplied labels do not independently verify the inputs.}
\item{seed}{Nonnegative integer seed for any randomness in
  \code{pool_builder}, using the caller's current RNG kinds. The generator's
  own seed is separate. The caller's RNG kinds and state, including absence
  of an initial random state, are restored on return or error. Box-Muller
  normal generation and user-supplied RNGs are rejected before mutation
  because their complete state cannot be restored; use \code{Inversion}
  normal generation instead.}
}
\details{
The initial domain is independent individuals, exact individual dosages and
a known genotype model. No covariate adjustment, ancestry model, kinship
model, LD estimation or sequencing observation model is provided here.

The target null for variant \eqn{j} is that the complete individual phenotype
vector \eqn{Y} is independent of \eqn{G_j} conditional on the remaining
genotypes \eqn{G_{-j}}. It is not a causal conclusion, nor a null for a separate
variant-trait pair. Multitrait importance is the Euclidean norm of the
trait-standardized coefficient vector. The contrast is original importance
minus knockoff importance. Identical genotype copies have contrast exactly zero.

All complete pairs are retained; no outcome screening or MAF filter is applied.
The two matrices share one pool plan, equal DNA contributions within each pool,
and exactly the same preprocessing and solver rules. Internally the individuals
and variant pairs are ordered by ID. Within each unordered pair, the first
distinct individual dosage determines which column is fitted first. Synthetic
column IDs and a returned mapping avoid collisions with user-supplied variant
names. Consequently the numerical matrix passed to the solver is identical
after any original/knockoff swaps, while the score mapping changes orientation.
This preserves the sign-flip property even when summation order or iterative
convergence otherwise depends on column order. Identical pairs are explicitly
assigned zero contrast.

For finite contrasts \eqn{W}, knockoff+ chooses the smallest positive value
among \eqn{|W_j|} such that
\deqn{\frac{1 + \#\{j:W_j\leq -t\}}{\max(1,\#\{j:W_j\geq t\})}\leq q.}
Variants with \eqn{W_j\geq t} are selected. If no threshold qualifies, the
threshold is infinity and the selected set is empty. Tied magnitudes are
evaluated together. At \eqn{q=0.05}, at least 20 positive contrasts must pass
the threshold even when no negative contrast passes its counterpart.

The Model-X assumptions require exchangeability of original/knockoff pairs
under every subset of pair swaps, and knockoffs independent of phenotypes
conditional on the originals. These imply joint invariance with the outcome
under swaps of null variants. A phenotype-only pooling map applied identically
to both matrices commutes with these swaps. Composing this map with a symmetric
fit and antisymmetric contrasts gives the sign-flip property used by knockoff+.
The argument concerns the original individual hypotheses; it does not assume
that the aggregated copies satisfy a new conditional independence property
given the aggregated originals.

These mathematical conditions cannot be certified from a finite supplied pair.
An estimated or misspecified LD model need not yield exact Model-X knockoffs.
Pooling noisy sequencing reads on only the original side also falls outside
this construction. The return status therefore remains
\code{"experimental_knockoff"} even for outputs of the reference generator.
The existing permutation calibrations retain their separate global-null scope.
}
\value{
A \code{poogwam_knockoff} list containing:
\item{W}{Named contrasts in the original input variant order.}
\item{threshold, selected, q, filter}{Knockoff+ threshold, selected variant
  IDs, target level and filter name.}
\item{importance, beta_standardized, knockoff_beta_standardized}{Paired
  importance scores and trait-standardized coefficients, mapped back to
  original variant IDs.}
\item{fit, pair_map, pool_design}{The augmented fit, mapping from external
  variant IDs to internal original and knockoff columns, and the realized
  pool plan. The internal \code{fit} is exploratory and is not a calibration.}
\item{inference_status, inference_scope, inference_unit}{Explicit experimental
  status and the conditional, variant-level target for the full trait vector.}
\item{generator_provenance, model_validation, allele_metadata, allele_validation}{
  Declared provenance and allele checks. \code{model_validation} is always
  \code{"not_certified"}.}
\item{diagnostics, rng, R_version, package_version, source_mode, assumptions, limitations}{Pair and solver
  diagnostics, RNG kinds and states before and after pooling, R version,
  package version, whether code is loaded from its namespace or sourced,
  and the conditions required for interpretation. The package version is
  missing for directly sourced development code.}
\item{fingerprint_method, input_fingerprint, settings_fingerprint}{Provenance
  checksums of the aligned individual inputs and realized plan/settings,
  using MD5 of uncompressed R serialization version 2. These are checksums,
  not certificates of the model or of input authenticity.}
No p-values or unconditional guarantee for externally supplied pairs are returned.
}
\references{
Candes, E., Fan, Y., Janson, L. and Lv, J. (2018).
Panning for gold: Model-X knockoffs for high dimensional controlled variable
selection. \emph{Journal of the Royal Statistical Society B}, 80, 551--577.
\doi{10.1111/rssb.12265}.

Barber, R. F., Candes, E. J. and Samworth, R. J. (2020).
Robust inference with knockoffs. \emph{Annals of Statistics}, 48, 1409--1431.
\doi{10.1214/19-AOS1852}.
}
\seealso{\code{\link{generate_independent_knockoffs}},
  \code{\link{fit_pool_ridge}}, \code{\link{fit_pool_multitrait}},
  \code{\link{calibrate_pool_permutation}}}
\examples{
# A small demonstration with a known product-binomial genotype model.
ids <- paste0("i", 1:40)
variants <- paste0("v", 1:12)
f <- setNames(rep(0.3, 12), variants)
set.seed(21)
g <- matrix(rbinom(40 * 12, 2, 0.3), 40, 12,
            dimnames = list(ids, variants))
y <- matrix(g[, 1] + rnorm(40), 40, 1, dimnames = list(ids, "trait"))
copies <- generate_independent_knockoffs(g, f, seed = 22)
builder <- function(y) design_pools(y, k = 4, min_size = 2)
selection <- select_pool_knockoffs(g, copies, y, builder, q = 0.1)
selection$selected
selection$inference_status
}
