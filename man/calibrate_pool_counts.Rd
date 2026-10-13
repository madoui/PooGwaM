\name{calibrate_pool_counts}
\alias{calibrate_pool_counts}
\title{Experimental homogeneous pool-count bootstrap}
\usage{
calibrate_pool_counts(data, B = 999, seed = 1, fit_args = list())
}
\arguments{
\item{data}{A \code{pool_gwas_data} object with aligned observed \code{counts$alt}, positive \code{counts$depth}, and explicit \code{counts$type}, either \code{"reads"} or \code{"chromosomes"}. Frequencies must equal \code{alt/depth}.}
\item{B}{Positive integer number of bootstrap replicates.}
\item{seed}{Nonnegative integer seed; caller RNG kind and state are restored.}
\item{fit_args}{Uniquely named list of fixed arguments for \code{fit_pool_ridge}, excluding \code{data}.}
}
\description{Simulates an explicitly restricted global-null model with biological chromosome sampling and, for sequencing reads, a second observation stage. This is an experimental parametric bootstrap, not a validated general Pool-seq association test.}
\details{
For variant j, estimate a common null frequency as the pool-size-weighted mean of observed frequencies. Holding phenotypes, pool sizes and depths fixed, generate \eqn{A_{kj} \sim Binomial(2 n_k, p_j)} independently across pools and variants. For \code{type="reads"}, then generate \eqn{R_{kj} \sim Binomial(C_{kj}, A_{kj}/(2 n_k))}; the null frequencies are \eqn{R_{kj}/C_{kj}}. For \code{type="chromosomes"}, use \eqn{A_{kj}/(2 n_k)} directly and require observed depths exactly equal to twice each pool size. Biological sampling is not counted twice.

The model assumes unrelated diploid individuals, homogeneous allele frequencies, independent chromosomes, independent variants, equal DNA contributions, and no sequencing error or extra dispersion. It does not preserve LD. Covariates and structure are unsupported and covariate-bearing inputs are rejected. The model conditions on the observed phenotype-based pool design under independence. Positive integer pool sizes and depths are required.

The null allele frequencies are plug-in estimates; the pool-size-weighted estimator is not claimed to be a likelihood optimum for read data. The Monte Carlo formula \code{(1+exceedances)/(B+1)} prevents zero p-values but does not establish finite-sample validity. Model checking, LD/structure modelling and calibration on independent null data are required before confirmatory use. No FDR, strong FWER or partial-null guarantee is claimed. Maximum-statistic values share these restrictions.
}
\value{A \code{poogwam_calibration}, as described for \code{calibrate_pool_permutation}, with \code{inference_status="experimental_parametric_global_null"}, null frequency estimates, count type, assumptions and limitations.}
\seealso{\code{\link{calibrate_pool_permutation}}, \code{\link{association_table}}}
