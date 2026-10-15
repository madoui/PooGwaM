\name{fit_pool_multitrait}
\alias{fit_pool_multitrait}
\title{Pooled sparse-group multitrait regression}
\usage{
fit_pool_multitrait(data, lambda_group = 0.1, lambda_l1 = 0,
  lambda_ridge = 1, penalty = c("within_pool", "identity"),
  penalty_floor = 1e-8, standardize_traits = TRUE,
  max_iter = 10000L, tolerance = 1e-8, block_size = 1000L)
}
\arguments{
\item{data}{A \code{pool_gwas_data} object, optionally constructed with an
  encoded numeric covariate matrix or a formula and identified metadata;
  see \code{\link{as_pool_gwas_data}}. The fixed-effect design must have full
  rank after pooling and leave a positive residual pool dimension.}
\item{lambda_group}{Nonnegative penalty on the Euclidean norm of each variant's vector of effects.}
\item{lambda_l1}{Nonnegative penalty on the absolute value of each individual effect.}
\item{lambda_ridge}{Strictly positive quadratic regularization; guarantees uniqueness of penalized effects.}
\item{penalty}{Diagonal quadratic penalty, as in \code{fit_pool_ridge}.}
\item{penalty_floor}{Strictly positive minimum diagonal penalty.}
\item{standardize_traits}{Standardize traits for shared selection; defaults to \code{TRUE}. Effects are returned in original units.}
\item{max_iter}{Positive integer iteration limit.}
\item{tolerance}{Strictly positive relative proximal KKT tolerance.}
\item{block_size}{Positive integer number of variants per block.}
}
\description{Fits variant-wise shared selection with optional trait-specific sparsity, an unpenalized intercept and covariates.}
\details{The criterion is one half of the weighted residual sum of squares plus one half of \code{lambda_ridge} times the diagonal quadratic penalty, plus \code{lambda_group} times the sum of row Euclidean norms, plus \code{lambda_l1} times the sum of absolute coefficients. The solver uses accelerated proximal gradient with backtracking and a monotone restart. The soft-thresholding step is followed by row group shrinkage. Convergence is assessed using the unit-step proximal-gradient mapping, normalized by the maximum of one and the largest absolute working coefficient. When both sparse penalties are zero the exact ridge solution is used. If the iteration limit is reached a warning is issued and \code{diagnostics$converged} is false; such fits should not be used for calibrated testing. Penalty choices should be fixed independently or selected inside a valid resampling scheme. No p-values are generated.}
\value{A \code{poogwam_fit} object with the same coefficient and prediction fields as \code{fit_pool_ridge}, including \code{covariate_effects}, \code{covariate_spec} for formula inputs and \code{diagnostics$covariate_design}. Additional diagnostics include convergence, KKT residual, objective history, line-search steps, and selected variants. The objective and sparse penalties operate on standardized phenotype units by default. The intercept and fixed covariates are never subject to the sparse or genetic ridge penalties; their reported effects are in original phenotype units. Nonzero genetic coefficients do not by themselves imply statistically significant associations.}
\seealso{\code{\link{fit_pool_ridge}}, \code{\link{predict.poogwam_fit}}}
