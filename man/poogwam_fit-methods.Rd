\name{poogwam_fit-methods}
\alias{print.poogwam_fit}
\alias{predict.poogwam_fit}
\title{Print and predict pooled regression fits}
\usage{
\method{print}{poogwam_fit}(x, ...)
\method{predict}{poogwam_fit}(object, newdata = NULL, covariates = NULL, ...)
}
\arguments{
\item{x, object}{A fit returned by \code{fit_pool_ridge} or \code{fit_pool_multitrait}.}
\item{newdata}{A pool-by-variant frequency matrix with unique row and column IDs, or a \code{pool_gwas_data} object. \code{NULL} returns training predictions. The complete training variant set is required.}
\item{covariates}{Prediction covariates with matching pool IDs and the training covariate names. For a \code{pool_gwas_data} object these are taken from \code{newdata}.}
\item{...}{Additional arguments; unused prediction arguments are rejected.}
}
\description{Display an exploratory fit or predict pool phenotype means in their original units. Matrices are aligned by identifiers. Variant allele orientation must match the training data.}
\value{Printing returns the fit invisibly. Prediction returns a pool-by-trait numeric matrix.}
