\name{association_table}
\alias{association_table}
\title{Effects, scores and separately adjusted calibration results}
\usage{
association_table(fit, calibration = NULL,
  adjustment = c("BH", "BY", "holm", "none"))
}
\arguments{
\item{fit}{A \code{poogwam_fit} from \code{fit_pool_ridge}.}
\item{calibration}{Optional matching \code{poogwam_calibration}. Data, settings, coefficients and standardized coefficients must match exactly. Using \code{calibration$fit} is the simplest way to ensure this.}
\item{adjustment}{Multiple-testing arithmetic applied once to raw p-values. Available methods are BH, BY, Holm or none. BH is not a claim that its dependence assumptions hold.}
}
\description{Builds a long table retaining effect estimates, trait and omnibus scores, raw p-values and corrected p-values in distinct columns. Without calibration, every p-value is missing and status is exploratory.}
\details{
Trait scores are absolute ridge effects. Omnibus scores are norms of standardized effects, and omnibus rows have no scalar effect. There are two separate declared testing families: all variant-by-trait pairs, and all variant omnibus scores. Each family is adjusted once from its raw p-values. Jointly selecting discoveries across these two families does not inherit either family's error control.

The \code{p_maxT} column is supplied only on omnibus rows and is copied from the null maximum distribution. It is not readjusted. Global-null calibration does not establish individual causal-null validity, FDR under partial alternatives or strong FWER. BY or Holm arithmetic cannot repair invalid raw p-values. Inspect \code{inference_status}, \code{inference_scope}, Monte Carlo resolution and calibration assumptions before interpretation.
}
\value{A data frame with \code{variant_id}, \code{trait}, \code{test}, \code{effect}, \code{effect_standardized}, \code{score}, \code{p_raw}, \code{p_adjusted}, \code{p_maxT}, \code{test_family}, \code{adjustment}, \code{inference_status}, \code{inference_scope} and \code{p_resolution}. An \code{inference_note} attribute explains restrictions.}
\seealso{\code{\link{calibrate_pool_permutation}}, \code{\link{calibrate_pool_counts}}}
