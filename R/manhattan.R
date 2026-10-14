# Plotting is optional. Preparing a plot never creates or recalibrates p-values.
manhattan_data <- function(x, variants = NULL,
                           value = c("p_maxT", "p_raw", "p_adjusted", "score", "W"),
                           test = c("omnibus", "trait"), trait = NULL,
                           chromosome = "chromosome", position = "position",
                           chromosome_order = NULL, na_action = c("error", "omit")) {
  value <- match.arg(value)
  test <- match.arg(test)
  na_action <- match.arg(na_action)
  scalar_string <- function(z) is.character(z) && length(z) == 1L &&
    !is.na(z) && nzchar(z)
  if (!scalar_string(chromosome) || !scalar_string(position) || chromosome == position)
    stop("chromosome and position must name distinct metadata columns.", call. = FALSE)
  if (!is.null(trait) && !scalar_string(trait))
    stop("trait must be one nonempty trait name.", call. = FALSE)
  if (inherits(x, "poogwam_calibration")) {
    if (is.null(variants)) variants <- x$fit$data$variants
    tab <- association_table(x$fit, x)
  } else if (inherits(x, "poogwam_fit")) {
    if (is.null(variants)) variants <- x$data$variants
    tab <- association_table(x)
  } else if (inherits(x, "poogwam_knockoff")) {
    if (value != "W" || test != "omnibus" || !is.null(trait))
      stop("Knockoff plots require value='W', test='omnibus' and no trait.", call. = FALSE)
    tab <- data.frame(variant_id = names(x$W), W = as.numeric(x$W),
                      inference_status = "experimental_knockoff")
  } else {
    if (!is.data.frame(x)) stop("x must be a fit, calibration, knockoff result or data frame.", call. = FALSE)
    tab <- x
  }
  if (!all(c("variant_id", value) %in% names(tab)))
    stop("The results require variant_id and the requested value column.", call. = FALSE)
  if ("test" %in% names(tab)) tab <- tab[!is.na(tab$test) & tab$test == test, , drop = FALSE]
  if (test == "omnibus" && !is.null(trait))
    stop("trait is only used with test='trait'.", call. = FALSE)
  if (test == "trait") {
    if (!"trait" %in% names(tab)) stop("Trait plots require a trait column.", call. = FALSE)
    choices <- unique(as.character(tab$trait[!is.na(tab$trait)]))
    if (is.null(trait)) {
      if (length(choices) != 1L) stop("Choose one trait explicitly for a trait plot.", call. = FALSE)
      trait <- choices
    }
    tab <- tab[!is.na(tab$trait) & tab$trait == trait, , drop = FALSE]
  }
  if (!nrow(tab)) stop("No results match the requested test and trait.", call. = FALSE)
  ids <- as.character(tab$variant_id)
  .pg_ids(ids, nrow(tab), "Result variants")
  if (is.null(variants)) variants <- tab
  if (!is.data.frame(variants) ||
      !all(c("variant_id", chromosome, position) %in% names(variants)))
    stop("variants requires variant_id, chromosome and position columns (or their supplied names).",
         call. = FALSE)
  map_ids <- as.character(variants$variant_id)
  .pg_ids(map_ids, nrow(variants), "Metadata variants")
  if (!all(ids %in% map_ids)) stop("Missing coordinates for some result variant IDs.", call. = FALSE)
  map <- variants[match(ids, map_ids), , drop = FALSE]
  chr <- as.character(map[[chromosome]])
  bp <- map[[position]]
  val <- tab[[value]]
  if (!is.numeric(bp) || !is.numeric(val))
    stop("Positions and plot values must be numeric.", call. = FALSE)
  missing <- is.na(chr) | !nzchar(trimws(chr)) | is.na(bp) | is.na(val)
  omitted <- ids[missing]
  if (any(missing) && na_action == "error")
    stop("Missing coordinates or values; use na_action='omit' to exclude explicitly. Uncalibrated fits have no p-values.",
         call. = FALSE)
  if (any(missing)) {
    warning(length(omitted), " variants omitted because coordinates or values are missing.", call. = FALSE)
    tab <- tab[!missing, , drop = FALSE]
    ids <- ids[!missing]; chr <- chr[!missing]; bp <- bp[!missing]; val <- val[!missing]
  }
  if (!length(ids)) stop("No complete variants remain for plotting.", call. = FALSE)
  if (any(!is.finite(bp) | bp < 0 | bp != floor(bp)))
    stop("Positions must be finite, nonnegative base-pair integers.", call. = FALSE)
  if (any(!is.finite(val))) stop("Plot values must be finite.", call. = FALSE)
  logp <- value %in% c("p_maxT", "p_raw", "p_adjusted")
  if (logp && any(val <= 0 | val > 1))
    stop("P-values must be in (0, 1]; zero values are not silently clipped.", call. = FALSE)
  labels <- unique(chr)
  if (is.null(chromosome_order)) {
    bare <- sub("^chr", "", labels, ignore.case = TRUE)
    numeric_chr <- grepl("^[0-9]+$", bare)
    number <- rep(Inf, length(bare))
    number[numeric_chr] <- as.numeric(bare[numeric_chr])
    labels <- labels[order(!numeric_chr, number, labels, method = "radix")]
  } else {
    if (!is.character(chromosome_order) || anyNA(chromosome_order) ||
        anyDuplicated(chromosome_order) || !setequal(chromosome_order, labels))
      stop("chromosome_order must list each plotted chromosome label exactly once.", call. = FALSE)
    labels <- chromosome_order
  }
  plot_data <- data.frame(SNP = ids, CHR = match(chr, labels), BP = bp, P = val,
                          stringsAsFactors = FALSE)
  plot_data <- plot_data[order(plot_data$CHR, plot_data$BP, plot_data$SNP, method = "radix"), ]
  rownames(plot_data) <- NULL
  field <- function(name) if (name %in% names(tab)) unique(tab[[name]]) else NULL
  structure(list(data = plot_data, chromosome_labels = labels, value = value,
    logp = logp, test = test, trait = trait, omitted_variant_ids = omitted,
    inference_status = field("inference_status"), inference_scope = field("inference_scope"),
    maxT_method = field("maxT_method")), class = "poogwam_manhattan_data")
}

plot_manhattan <- function(x, variants = NULL,
                           value = c("p_maxT", "p_raw", "p_adjusted", "score", "W"),
                           test = c("omnibus", "trait"), trait = NULL,
                           backend = c("qqman", "CMplot"),
                           chromosome = "chromosome", position = "position",
                           chromosome_order = NULL, na_action = c("error", "omit"),
                           threshold = NULL, highlight = NULL, ...) {
  backend <- match.arg(backend)
  prepared <- manhattan_data(x, variants, value, test, trait, chromosome, position,
                             chromosome_order, na_action)
  if (!is.null(threshold) && (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || (prepared$logp && (threshold <= 0 || threshold > 1))))
    stop("threshold must be finite; for p-values it must be in (0,1].", call. = FALSE)
  if (!is.null(highlight) && (!is.character(highlight) || anyNA(highlight) ||
      any(!highlight %in% prepared$data$SNP)))
    stop("highlight must contain only plotted variant IDs.", call. = FALSE)
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || anyNA(names(dots)) ||
      any(!nzchar(names(dots))) || anyDuplicated(names(dots))))
    stop("Additional plot arguments must have unique names.", call. = FALSE)
  reserved <- c("x", "Pmap", "chr", "bp", "p", "snp", "logp", "LOG10",
                "chrlabs", "chr.labels", "suggestiveline", "genomewideline",
                "threshold", "highlight", "file.output", "plot.type", "multraits", "multracks")
  if (any(names(dots) %in% reserved))
    stop("Use the explicit plotting arguments for data, scales, thresholds and chromosome labels.", call. = FALSE)
  if (backend == "CMplot" && length(prepared$chromosome_labels) < 2L)
    stop("Use backend='qqman' for a single chromosome. CMplot's single-chromosome axis can lose labels or fail on short intervals.",
         call. = FALSE)
  if (!requireNamespace(backend, quietly = TRUE))
    stop("Install the optional plotting package with install.packages('", backend, "').", call. = FALSE)
  y <- if (prepared$logp) -log10(prepared$data$P) else prepared$data$P
  yrange <- range(c(0, y, if (!is.null(threshold))
    if (prepared$logp) -log10(threshold) else threshold))
  if (diff(yrange) == 0) yrange <- yrange + c(0, 1)
  else {
    padding <- .05 * diff(yrange)
    yrange[2L] <- yrange[2L] + padding
    if (yrange[1L] < 0) yrange[1L] <- yrange[1L] - padding
  }
  default_label <- if (prepared$logp) paste0("-log10(", prepared$value, ")") else
    if (prepared$value == "W") "Knockoff W (experimental)" else "Exploratory score"
  label <- if (prepared$test == "trait") prepared$trait else "Omnibus"
  defaults <- list(main = paste(label, prepared$value), ylab = default_label, ylim = yrange)
  if (backend == "qqman") {
    # qqman 0.1.9 does not define its labels object for a single chromosome.
    defaults$xlab <- if (length(prepared$chromosome_labels) == 1L)
      paste("Chromosome", prepared$chromosome_labels, "position") else "Chromosome"
    args <- list(x = prepared$data, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
      logp = prepared$logp, chrlabs = if (length(prepared$chromosome_labels) > 1L)
        prepared$chromosome_labels else NULL,
      suggestiveline = FALSE, genomewideline = FALSE, highlight = highlight)
    do.call(qqman::manhattan, c(args, dots, defaults[!names(defaults) %in% names(dots)]))
  } else {
    args <- list(Pmap = prepared$data, plot.type = "m", LOG10 = prepared$logp,
      chr.labels = prepared$chromosome_labels, file.output = FALSE,
      threshold = NULL, highlight = highlight, multracks = FALSE, multraits = FALSE)
    defaults$mar <- c(4, 6, 4, 2)
    defaults$lab.cex <- 1
    defaults$lab.font <- 1
    defaults$main.cex <- 1.2
    defaults$verbose <- FALSE
    args["chr.den.col"] <- list(NULL)
    if ("chr.den.col" %in% names(dots)) args$chr.den.col <- NULL
    do.call(CMplot::CMplot, c(args, dots, defaults[!names(defaults) %in% names(dots)]))
  }
  # Both backends interpret a threshold of zero as "no line". Draw the explicit
  # user threshold ourselves, also allowing signed W thresholds at zero.
  if (!is.null(threshold))
    graphics::abline(h = if (prepared$logp) -log10(threshold) else threshold, col = "red")
  prepared$backend <- backend
  prepared$threshold <- threshold
  invisible(prepared)
}
