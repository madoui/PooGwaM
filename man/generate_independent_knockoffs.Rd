\name{generate_independent_knockoffs}
\alias{generate_independent_knockoffs}
\title{Reference Knockoffs for a Known Independent Binomial Genotype Model}
\description{
Generates independent diploid dosage copies under a product-binomial model
with externally supplied, known allele probabilities. This reference
generator is intended for methodological examples and controlled simulations,
not for correlated genetic data with an estimated LD distribution.
}
\usage{
generate_independent_knockoffs(genotypes, frequencies, seed = 1)
}
\arguments{
\item{genotypes}{Numeric matrix of observed diploid 0/1/2 dosages, with at
  least two individuals and unique nonempty individual and variant names.
  Its dimensions and IDs determine the output. Observed dosages are validated
  but never used to estimate the supplied allele probabilities.}
\item{frequencies}{Numeric probabilities in [0,1], uniquely named by every
  genotype variant ID. These must be externally known population alternate
  allele probabilities under the specified model, not frequencies estimated
  from \code{genotypes}. The function cannot verify their external origin.
  At probabilities exactly zero or one, incompatible observed dosages are
  rejected as outside the model's support.}
\item{seed}{Nonnegative integer seed. The caller's current RNG kinds are
  used and recorded. RNG kinds and state are restored on return or error.
  Box-Muller normal generation and user-supplied RNGs are rejected before
  mutation because their complete state cannot be restored.}
}
\details{
The mathematical model assumes independent individuals and mutually
independent variants with \eqn{G_{ij}\sim\mathrm{Binomial}(2,f_j)}.
Independent copies from the same known product distribution then satisfy
pairwise Model-X exchangeability. Generation takes no phenotype argument.
If the actual genotypes exhibit LD, population structure, relatedness or a
different distribution, these copies do not automatically satisfy the model.
Supplying estimated probabilities does not establish exact finite-sample FDR
control, even if they were estimated without phenotypes.

Sampling follows canonical individual and variant ID order before returning
to the input order. Consequently reordering named inputs does not change the
paired sampled values at a fixed seed and RNG configuration. Frequencies,
model assumptions, seed, RNG kinds and states, R and package versions, and
input/settings/output serialization checksums are retained
in the matrix attribute \code{knockoff_provenance}. This attribute is a record
of how copies were generated, not a certificate for the original data.
}
\value{
A numeric matrix of diploid knockoff dosages with exactly the input dimensions
and dimnames, and a \code{knockoff_provenance} attribute. It can be supplied
directly as \code{knockoffs} to \code{\link{select_pool_knockoffs}}.
}
\seealso{\code{\link{select_pool_knockoffs}}}
\examples{
g <- matrix(c(0, 1, 2, 0, 1, 1, 0, 2), 4, 2,
            dimnames = list(letters[1:4], c("snp1", "snp2")))
copies <- generate_independent_knockoffs(g,
  frequencies = c(snp1 = 0.25, snp2 = 0.5), seed = 10)
copies
}
