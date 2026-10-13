#' Evaluate associations against trait-specific truth
#'
#' Inputs are named variant-by-trait matrices (or named vectors for one trait).
#' The type of the supplied values must be stated explicitly. Scores remain
#' scores: only raw p-values are adjusted. Selection of scores is optional and
#' exploratory. At a top-n boundary, every tied score is retained.
#'
#' @param values Numeric named matrix or vector of scores or p-values.
#' @param truth Logical or binary matrix/vector with the same variant and trait IDs.
#' @param input_type Required: "score", "p_raw", or "p_adjusted".
#' @param alpha Decision threshold for p-values in [0,1].
#' @param adjustment Adjustment for raw p-values; NULL selects BH. Must be NULL
#'   for scores and already adjusted p-values.
#' @param family Adjust across each trait or across all variant-trait tests.
#' @param higher_better For scores, whether larger values rank first.
#' @param score_threshold Optional score selection threshold.
#' @param top_n Optional number of leading scores per trait, including ties.
#' @return A list with a summary per trait, a long results table, threshold
#'   curves, and explicit evaluation metadata. AUROC is NA with either class
#'   absent; average precision and recall are NA without positives. FDP is zero
#'   with zero selections. These are descriptive metrics for this dataset,
#'   not a repeated-experiment FDR estimate.
#' @export
evaluate_associations <- function(values, truth, input_type, alpha = 0.05,
                                  adjustment = NULL,
                                  family = c("trait", "global"),
                                  higher_better = TRUE,
                                  score_threshold = NULL, top_n = NULL) {
  if (missing(input_type)) stop("input_type must explicitly identify score, p_raw, or p_adjusted")
  input_type <- match.arg(input_type, c("score", "p_raw", "p_adjusted"))
  family <- match.arg(family)
  to_matrix <- function(x, label) {
    if (is.null(dim(x))) {
      x <- matrix(x, ncol = 1L, dimnames = list(names(x), "trait"))
    }
    if (!is.matrix(x) || length(dim(x)) != 2L || any(dim(x) == 0L))
      stop(label, " must be a nonempty matrix or a named vector")
    for (ids in dimnames(x)) {
      if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
        stop(label, " requires unique nonempty variant and trait IDs")
    }
    if (is.null(rownames(x)) || is.null(colnames(x)))
      stop(label, " requires unique nonempty variant and trait IDs")
    x
  }
  values <- to_matrix(values, "values")
  truth <- to_matrix(truth, "truth")
  if (!is.numeric(values) || any(!is.finite(values)))
    stop("values must be finite numeric values")
  if (!(is.logical(truth) || is.numeric(truth)) || anyNA(truth) ||
      any(!truth %in% c(0, 1))) stop("truth must contain only logical or binary values")
  if (!setequal(rownames(values), rownames(truth)) ||
      !setequal(colnames(values), colnames(truth)))
    stop("values and truth must have exactly matching variant and trait IDs")
  truth <- truth[rownames(values), colnames(values), drop = FALSE] != 0
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1)
    stop("alpha must be a single number in [0,1]")
  if (!is.logical(higher_better) || length(higher_better) != 1L || is.na(higher_better))
    stop("higher_better must be TRUE or FALSE")
  if (input_type != "score") {
    if (any(values < 0 | values > 1)) stop("p-values must lie in [0,1]")
    if (!is.null(score_threshold) || !is.null(top_n))
      stop("score_threshold and top_n are only available for scores")
  }
  if (input_type != "p_raw" && !is.null(adjustment))
    stop("adjustment is only allowed for raw p-values; adjusted p-values are never readjusted")
  if (!is.null(score_threshold) && (!is.numeric(score_threshold) ||
      length(score_threshold) != 1L || !is.finite(score_threshold)))
    stop("score_threshold must be a finite scalar")
  if (!is.null(top_n) && (!is.numeric(top_n) || length(top_n) != 1L ||
      !is.finite(top_n) || top_n < 0 || top_n > nrow(values) || top_n != floor(top_n)))
    stop("top_n must be an integer between zero and the number of variants")
  if (!is.null(top_n) && !is.null(score_threshold))
    stop("choose either top_n or score_threshold")
  p_adjusted <- matrix(NA_real_, nrow(values), ncol(values), dimnames = dimnames(values))
  if (input_type == "p_raw") {
    if (is.null(adjustment)) adjustment <- "BH"
    adjustment <- match.arg(adjustment, stats::p.adjust.methods)
    if (family == "global") {
      p_adjusted[] <- stats::p.adjust(as.vector(values), method = adjustment)
    } else {
      for (j in seq_len(ncol(values)))
        p_adjusted[, j] <- stats::p.adjust(values[, j], method = adjustment)
    }
  } else if (input_type == "p_adjusted") p_adjusted <- values

  summaries <- rows <- curves <- vector("list", ncol(values))
  names(curves) <- colnames(values)
  divide <- function(x, n) if (n == 0) NA_real_ else x / n
  for (j in seq_len(ncol(values))) {
    observed <- values[, j]
    positive <- truth[, j]
    score <- if (input_type == "score" && higher_better) observed else -observed
    npos <- sum(positive)
    nneg <- length(positive) - npos
    # Evaluate distinct score thresholds, making both AUC and AP invariant to
    # arbitrary row ordering within ties.
    ord <- order(score, decreasing = TRUE)
    ends <- c(which(diff(score[ord]) != 0), length(score))
    tp_curve <- cumsum(positive[ord])[ends]
    fp_curve <- ends - tp_curve
    precision <- tp_curve / ends
    recall <- if (npos == 0) rep(NA_real_, length(ends)) else tp_curve / npos
    ap <- if (npos == 0) NA_real_ else sum(diff(c(0, recall)) * precision)
    auroc <- if (npos == 0 || nneg == 0) NA_real_ else
      (sum(rank(score, ties.method = "average")[positive]) - npos * (npos + 1) / 2) / (npos * nneg)
    selected <- rep(NA, length(score))
    if (input_type != "score") selected <- p_adjusted[, j] <= alpha
    if (!is.null(score_threshold)) selected <- if (higher_better)
      observed >= score_threshold else observed <= score_threshold
    if (!is.null(top_n)) selected <- if (top_n == 0) rep(FALSE, length(score)) else
      score >= sort(score, decreasing = TRUE)[top_n]
    if (anyNA(selected)) {
      tp <- fp <- tn <- fn <- nselected <- fdp <- sensitivity <- fpr <- precision_selected <- NA_real_
    } else {
      tp <- sum(selected & positive)
      fp <- sum(selected & !positive)
      tn <- sum(!selected & !positive)
      fn <- sum(!selected & positive)
      nselected <- tp + fp
      fdp <- if (nselected == 0) 0 else fp / nselected
      sensitivity <- divide(tp, npos)
      fpr <- divide(fp, nneg)
      precision_selected <- divide(tp, nselected)
    }
    summaries[[j]] <- data.frame(trait = colnames(values)[j], n_variants = length(score),
      n_positive = npos, n_negative = nneg, n_selected = nselected,
      TP = tp, FP = fp, TN = tn, FN = fn, FDP = fdp, recall = sensitivity,
      precision = precision_selected, false_positive_rate = fpr,
      AUROC = auroc, average_precision = ap, stringsAsFactors = FALSE)
    rows[[j]] <- data.frame(variant_id = rownames(values), trait = colnames(values)[j],
      truth = positive, value = observed, rank = rank(-score, ties.method = "min"),
      p_raw = if (input_type == "p_raw") observed else NA_real_,
      p_adjusted = p_adjusted[, j], selected = selected, stringsAsFactors = FALSE)
    curves[[j]] <- data.frame(threshold = observed[ord][ends], n_selected = ends,
      TP = tp_curve, FP = fp_curve, precision = precision, recall = recall,
      false_positive_rate = if (nneg == 0) NA_real_ else fp_curve / nneg)
  }
  list(summary = do.call(rbind, summaries), results = do.call(rbind, rows), curves = curves,
    metadata = list(input_type = input_type, adjustment = adjustment, family = family,
      alpha = if (input_type == "score") NULL else alpha,
      higher_better = if (input_type == "score") higher_better else FALSE,
      score_threshold = score_threshold, top_n = top_n,
      score_selection = if (input_type == "score") "exploratory; all boundary ties retained" else NULL,
      FDP_definition = "false selections / max(1, all selections); one dataset, not FDR",
      average_precision_definition = "sum of recall increments times precision at distinct thresholds"))
}
