\name{evaluate_associations}
\alias{evaluate_associations}
\title{Evaluate associations against trait-specific truth}
\usage{
evaluate_associations(values, truth, input_type, alpha = 0.05,
  adjustment = NULL, family = c("trait", "global"),
  higher_better = TRUE, score_threshold = NULL, top_n = NULL)
}
\arguments{
\item{values}{Numeric matrix with unique variant row IDs and trait column IDs, or named vector.}
\item{truth}{Binary or logical matrix/vector with exactly matching IDs, in any order.}
\item{input_type}{Explicitly \code{"score"}, \code{"p_raw"}, or \code{"p_adjusted"}.}
\item{alpha}{Threshold for p-values.}
\item{adjustment}{For raw p-values only, a \code{p.adjust} method. NULL defaults to BH.}
\item{family}{Adjustment within traits or across all tests.}
\item{higher_better}{Whether larger scores rank first. Ignored for p-values.}
\item{score_threshold}{Optional threshold for exploratory score selection.}
\item{top_n}{Optional leading scores per trait; includes all ties at the boundary.}
}
\description{Compute ranking and selection metrics without confusing scores with p-values or adjusting p-values twice. The truth can describe the scientific target of interest (for example, fixed effects) but must not be interpreted as strict null hypotheses unless the data-generating model justifies that interpretation.}
\details{Average precision is the sum of recall increments times precision at distinct score thresholds. AUROC uses average ranks for ties. AUROC is NA when either class is absent; recall and average precision are NA without positive truth. FDP is zero when there are no selections. A single dataset's FDP does not estimate FDR across repeated experiments. Score selections are exploratory and never receive p-values.}
\value{List containing \code{summary} (per-trait metrics), \code{results} (variant-trait table), \code{curves} (distinct thresholds), and \code{metadata}.}
\examples{
scores <- setNames(c(4, 3, 1, 2), c("s1", "s2", "s3", "s4"))
truth <- setNames(c(TRUE, FALSE, FALSE, TRUE), names(scores))
evaluate_associations(scores, truth, "score", top_n = 2)$summary
}
