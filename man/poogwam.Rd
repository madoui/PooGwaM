\name{poogwam}
\alias{poogwam}
\title{Fit the PooGwaM V2 pooled model}
\description{Dispatches to deterministic ridge or coupled sparse-group multitrait estimation. Returns an exploratory fit with explicit diagnostics.}
\usage{poogwam(data, method = c("ridge", "multitrait"), ...)}
\arguments{
\item{data}{Validated \code{pool_gwas_data}.}
\item{method}{Either separable multi-response ridge or sparse-group multitrait estimation.}
\item{...}{Arguments to \code{fit_pool_ridge} or \code{fit_pool_multitrait}.}
}
\details{V2 no longer converts the distribution of coefficients across SNPs into Gaussian p-values. See \code{association_table} for output and the explicitly bounded calibration interfaces. GWAlpha, Chi-square and individual GWAS are external comparators and are not dispatched by this function. The V1 positional API has been replaced intentionally; see the migration guide.}
\value{A \code{poogwam_fit} containing variant-by-trait effects in original phenotype units, fitted pool means, diagnostics, settings and input data.}
\seealso{\code{\link{as_pool_gwas_data}}, \code{\link{fit_pool_ridge}}, \code{\link{fit_pool_multitrait}}, \code{\link{association_table}}}
