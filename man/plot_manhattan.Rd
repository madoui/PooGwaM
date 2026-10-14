\name{plot_manhattan}
\alias{plot_manhattan}
\alias{manhattan_data}
\title{Prepare and display Manhattan plots with qqman or CMplot}
\description{Aligns genomic coordinates by variant ID and displays one explicitly chosen test family. Plotting uses optional packages; preparing data requires no plotting dependency. No p-values are calculated, adjusted or inferred by these functions.}
\usage{
manhattan_data(x, variants = NULL,
  value = c("p_maxT", "p_raw", "p_adjusted", "score", "W"),
  test = c("omnibus", "trait"), trait = NULL,
  chromosome = "chromosome", position = "position",
  chromosome_order = NULL, na_action = c("error", "omit"))

plot_manhattan(x, variants = NULL,
  value = c("p_maxT", "p_raw", "p_adjusted", "score", "W"),
  test = c("omnibus", "trait"), trait = NULL,
  backend = c("qqman", "CMplot"),
  chromosome = "chromosome", position = "position",
  chromosome_order = NULL, na_action = c("error", "omit"),
  threshold = NULL, highlight = NULL, ...)
}
\arguments{
\item{x}{A calibration, fit, knockoff result, or data frame containing \code{variant_id} and the requested value column. Association tables can contain \code{test} and \code{trait}.}
\item{variants}{Metadata with unique \code{variant_id}, chromosome and position columns. Supersets are allowed, but every selected result must match. If absent, uses fit metadata or columns in the supplied data frame.}
\item{value}{Value to display. P-values use minus log10; \code{score} and signed \code{W} use their original scale. Knockoff objects require \code{value="W"}. Uncalibrated fits require \code{value="score"}.}
\item{test}{Select omnibus results or results for one trait. The two families are never combined.}
\item{trait}{Trait name. Required when multiple traits are available and \code{test="trait"}.}
\item{backend}{Optional plotting package to use. Install it separately with \code{install.packages}. Use \code{qqman} for a single chromosome; this wrapper restricts CMplot to at least two chromosomes because its single-chromosome axis can fail on short intervals or discard labels.}
\item{chromosome, position}{Metadata column names. Positions must be numeric, finite, nonnegative integers in base pairs. Chromosome labels may be arbitrary strings.}
\item{chromosome_order}{Optional character vector listing each plotted chromosome exactly once. By default, numeric labels (optionally prefixed by chr) are ordered numerically, followed by other labels lexically. No species-specific chromosome numbering is assumed.}
\item{na_action}{Reject missing values, or explicitly omit them with a warning and recorded IDs. Missing identifier matches are always errors.}
\item{threshold}{Optional single red horizontal line in the original value units: a p-value in (0,1], or a finite score, including zero. No conventional GWAS threshold is added by default. It does not automatically amplify points; use \code{highlight} for that.}
\item{highlight}{Variant IDs to highlight, all of which must be plotted.}
\item{...}{Named backend styling arguments, such as \code{col}, \code{main} and \code{cex}. Arguments changing input, scales, output files or testing semantics must use the explicit interface.}
}
\details{Zero p-values are rejected rather than silently clipped. The returned coordinates use consecutive integer chromosome codes for the plotting backend and retain the original chromosome labels separately. A plot does not strengthen the inference status of its input. In particular, maxT global-null calibration and experimental knockoffs retain their stated limitations.

Plots are drawn on the current graphics device. CMplot is called with \code{file.output=FALSE}. To save a figure, open a device such as \code{pdf()} or \code{png()}, call this function, then close it with \code{dev.off()}.}
\value{\code{manhattan_data} returns a \code{poogwam_manhattan_data} list containing the sorted backend data (SNP, CHR, BP, P), chromosome labels, scale, chosen family, omitted IDs and available inference metadata. \code{plot_manhattan} invisibly returns the same list with backend and threshold. It draws a base graphics plot, not a ggplot object.}
\examples{
tab <- data.frame(variant_id = c("a", "b", "c"),
  chromosome = c("1", "1", "X"), position = c(100, 200, 100),
  p_maxT = c(.2, .01, .8))
manhattan_data(tab)
if (requireNamespace("qqman", quietly = TRUE)) {
  plot_manhattan(tab, threshold = .05)
}
}
\seealso{\code{\link{association_table}}, \code{\link{select_pool_knockoffs}}}
