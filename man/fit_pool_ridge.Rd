\name{fit_pool_ridge}
\alias{fit_pool_ridge}
\title{Exact pooled multi-response ridge regression}
\usage{
fit_pool_ridge(data, lambda = 1, penalty = c("within_pool", "identity"),
  penalty_floor = 1e-8, block_size = 1000L, standardize_traits = FALSE)
}
\arguments{
\item{data}{A \code{pool_gwas_data} object.}
\item{lambda}{Strictly positive ridge regularization.}
\item{penalty}{Diagonal ridge penalty: weighted mean of \eqn{2F(1-F)} or identity.}
\item{penalty_floor}{Strictly positive lower bound on every diagonal penalty.}
\item{block_size}{Positive integer number of variants per computational block.}
\item{standardize_traits}{Divide responses by their weighted population standard deviations before optimization. Reported effects are converted back to original units.}
}
\description{Fits all phenotype columns simultaneously using a deterministic dual system with one row per pool. Includes an unpenalized intercept and optional full-rank covariates.}
\details{Let \eqn{X=2F} denote alternate-allele dosage, \eqn{w_k=n_k/\sum n_k}, and \eqn{Z} contain the intercept and covariates. The criterion is the weighted residual sum of squares plus \eqn{\lambda\sum_{jt}D_j B_{jt}^2}. Nuisance effects are projected out using weighted QR; variants are processed by blocks and the dual system is \eqn{K\times K}. Variants constant across pools receive exactly zero coefficients. Covariates must be linearly independent of the intercept and leave at least one residual pool degree of freedom. The within-pool penalty is a regularization heuristic, not an assertion of Hardy-Weinberg equilibrium or a sequencing-error model. The method fits several responses but ridge itself does not share selection across traits. No p-values are generated.}
\value{A \code{poogwam_fit} object containing \code{beta} (variant by trait effects in original phenotype units per alternate allele), \code{beta_standardized} (effects divided by weighted phenotype standard deviations), \code{beta_working} (optimization units), \code{intercept}, \code{covariate_effects}, \code{fitted}, \code{residuals}, the optimization \code{objective}, \code{penalty}, \code{weights}, \code{trait_scales}, \code{diagnostics}, \code{settings}, and original \code{data}. Diagnostics include gradient, condition number, rank, constant traits and penalty flooring. Constant traits use scale one. The inference status is exploratory.}
\seealso{\code{\link{fit_pool_multitrait}}, \code{\link{predict.poogwam_fit}}}
