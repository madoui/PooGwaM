\name{poogwam_fit-methods}
\alias{print.poogwam_fit}
\alias{predict.poogwam_fit}
\title{Print and predict pooled regression fits}
\usage{
\method{print}{poogwam_fit}(x, ...)
\method{predict}{poogwam_fit}(object, newdata = NULL, covariates = NULL, ...,
  metadata = NULL, clusters = NULL)
}
\arguments{
\item{x, object}{A fit returned by \code{fit_pool_ridge} or \code{fit_pool_multitrait}.}
\item{newdata}{A pool-by-variant frequency matrix with unique row and column IDs, or a \code{pool_gwas_data} object. \code{NULL} returns training predictions. The complete training variant set is required.}
\item{covariates}{Optional already-encoded numeric prediction covariates with
  matching pool IDs and exactly the training covariate column names, without
  an intercept column. This explicit matrix interface is used with a frequency
  matrix in \code{newdata} and preserves compatibility with numeric training
  covariates. For a \code{pool_gwas_data} object, covariate information is taken
  from \code{newdata}. Do not combine an explicit numeric matrix with
  formula-prediction metadata.}
\item{metadata}{Named argument providing a data frame for prediction from a
  formula learned during training. Its explicit row names identify pools
  for a pool-level training formula, or individuals for an individual-level
  training formula. All variables required by the learned formula must be
  present. Encoding reuses the training terms, levels and contrasts; missing
  values and new factor levels are rejected. For a \code{pool_gwas_data} object,
  stored metadata are used unless this argument supplies replacement metadata.}
\item{clusters}{Named argument mapping every individual metadata row to its
  prediction pool, or a \code{pool_design}. Required with an individual-level
  training formula when \code{newdata} is a frequency matrix. With a
  \code{pool_gwas_data} object, stored clusters are used unless overridden.
  Individual encoding is followed by pooling; the resulting pool IDs must
  match \code{newdata}. Pool membership counts must also match pool sizes
  when \code{newdata} supplies them. Not used for pool-level metadata.}
\item{...}{Additional arguments; unused prediction arguments are rejected.}
}
\description{Display an exploratory fit or predict pool phenotype means in their original units. Matrices are aligned by identifiers. Variant allele orientation must match the training data.}
\details{
The fitted formula and encoding are stored in \code{object$covariate_spec}.
For an individual-level formula, nonlinear transformations and interactions
are evaluated before averaging individuals in the prediction pools, using
the learned preprocessing. A formula fitted to pool-level metadata instead
operates directly on pool rows. The two data levels are not interchangeable.

Fixed-effect coefficients are in \code{object$covariate_effects}, one row per
encoded covariate and one column per trait. The intercept is returned
separately in \code{object$intercept}; neither is genetically penalized.
These coefficients describe adjustment in the pooled model and do not by
themselves establish individual causal effects or significance.

\code{object$diagnostics$covariate_design} records the pooled design columns,
rank, dimensions, constant or aliased columns, input level, factor levels and
contrasts. Its \code{residual_df} is the number of pools minus the rank of the
intercept-plus-covariate design. It is not the residual degrees of freedom of
an ordinary least-squares model for deriving tests on penalized genetic effects.
Rank-deficient fixed effects or absence of residual pool dimension are rejected
during fitting, without silently selecting a subset of cofactor columns.

A model fitted from numeric covariates does not acquire a formula encoding
from a new formula-based data object. In that case, supply the frequency
matrix and an explicitly encoded numeric covariate matrix instead. Rank
errors during fitting have class \code{poogwam_covariate_rank_error} and
carry the design diagnostics in \code{covariate_diagnostics}.

\code{metadata} and \code{clusters} follow the dots and must be supplied by
name. Passing \code{newdata = NULL} returns training fitted values; it is not
a request to refit the model with alternative metadata. Covariate-adjusted
calibration and knockoff inference are outside this prediction extension.
}
\value{Printing returns the fit invisibly. Prediction returns a pool-by-trait numeric matrix.}
\examples{
pool_ids <- paste0("p", 1:6)
F <- matrix(c(.1, .3, .2, .6, .8, .7), 6, 1,
  dimnames = list(pool_ids, "snp1"))
Y <- matrix(c(10, 12, 11, 16, 18, 17), 6, 1,
  dimnames = list(pool_ids, "trait"))
metadata <- data.frame(
  site = factor(c("A", "A", "A", "B", "B", "B")),
  age = c(20, 30, 25, 22, 33, 28), row.names = pool_ids)
dat <- as_pool_gwas_data(F, Y, setNames(rep(10, 6), pool_ids),
  covariates = ~ site + age, metadata = metadata)
fit <- fit_pool_ridge(dat)
fit$covariate_effects
fit$diagnostics$covariate_design
predict(fit, newdata = F, metadata = metadata)
}
