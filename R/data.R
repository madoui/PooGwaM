# All data alignment is by ID; implicit positional matching is deliberately absent.
.pg_ids <- function(ids, n, what) {
  if (is.null(ids) || length(ids) != n || anyNA(ids) ||
      any(!nzchar(ids)) || anyDuplicated(ids))
    stop(what, " must have unique, nonempty IDs", call. = FALSE)
  as.character(ids)
}

.pg_matrix <- function(x, what, allow_na = FALSE) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1))))
      stop(what, " must contain only numeric columns", call. = FALSE)
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x) || !length(x) ||
      nrow(x) < 1L || ncol(x) < 1L)
    stop(what, " must be a nonempty numeric matrix with row and column IDs", call. = FALSE)
  .pg_ids(rownames(x), nrow(x), paste(what, "rows"))
  .pg_ids(colnames(x), ncol(x), paste(what, "columns"))
  if (any(is.infinite(x)) || (!allow_na && anyNA(x)))
    stop(what, " contains missing or nonfinite values", call. = FALSE)
  storage.mode(x) <- "double"
  x
}

.pg_align <- function(x, rows, cols = NULL, what) {
  if (!setequal(rownames(x), rows) || (!is.null(cols) && !setequal(colnames(x), cols)))
    stop(what, " IDs do not match", call. = FALSE)
  if (is.null(cols)) x[rows, , drop = FALSE] else x[rows, cols, drop = FALSE]
}

.pg_scalar <- function(x, what, lower = -Inf, upper = Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > upper || (integer && x != floor(x)))
    stop(what, " must be a finite ", if (integer) "integer " else "number ",
         "in [", lower, ", ", upper, "]", call. = FALSE)
  invisible(x)
}

as_pool_gwas_data <- function(frequencies = NULL, phenotypes, pool_sizes,
                              counts = NULL, covariates = NULL, variants = NULL,
                              missing = c("error", "drop_variants"),
                              min_maf = 0, min_depth = 0,
                              design = list(), provenance = list()) {
  missing <- match.arg(missing)
  .pg_scalar(min_maf, "min_maf", 0, .5)
  .pg_scalar(min_depth, "min_depth", 0)
  if (!is.list(design) || !is.list(provenance))
    stop("design and provenance must be lists", call. = FALSE)
  if (!is.null(counts)) {
    if (!is.list(counts) || !all(c("alt", "depth", "type") %in% names(counts)) ||
        !is.character(counts$type) || length(counts$type) != 1L || is.na(counts$type) ||
        !counts$type %in% c("reads", "chromosomes"))
      stop("counts requires alt, depth and explicit type 'reads' or 'chromosomes'", call. = FALSE)
    counts$alt <- .pg_matrix(counts$alt, "counts$alt", allow_na = TRUE)
    counts$depth <- .pg_matrix(counts$depth, "counts$depth", allow_na = TRUE)
    counts$depth <- .pg_align(counts$depth, rownames(counts$alt), colnames(counts$alt), "counts$depth")
    if (any(counts$alt < 0 | counts$depth < 0 | counts$alt > counts$depth, na.rm = TRUE) ||
        any(counts$alt != floor(counts$alt), na.rm = TRUE) ||
        any(counts$depth != floor(counts$depth), na.rm = TRUE))
      stop("Counts must be nonnegative integers with alt <= depth", call. = FALSE)
    if (is.null(frequencies)) {
      frequencies <- counts$alt / counts$depth
      frequencies[counts$depth == 0 & !is.na(counts$depth)] <- NA_real_
    }
  } else if (is.null(frequencies)) {
    stop("Provide frequencies or explicit counts", call. = FALSE)
  } else if (min_depth > 0) {
    stop("min_depth requires counts", call. = FALSE)
  }
  F <- .pg_matrix(frequencies, "frequencies", allow_na = TRUE)
  if (nrow(F) < 2L) stop("At least two occupied pools are required", call. = FALSE)
  if (any(F < 0 | F > 1, na.rm = TRUE)) stop("Frequencies must lie in [0,1]", call. = FALSE)
  pool_ids <- rownames(F)
  variant_ids <- colnames(F)
  Y <- .pg_align(.pg_matrix(phenotypes, "phenotypes"), pool_ids, what = "phenotypes")
  if (!is.numeric(pool_sizes)) stop("pool_sizes must be numeric", call. = FALSE)
  .pg_ids(names(pool_sizes), length(pool_sizes), "pool_sizes")
  if (!setequal(names(pool_sizes), pool_ids)) stop("pool_sizes IDs do not match", call. = FALSE)
  n <- setNames(as.numeric(pool_sizes[pool_ids]), pool_ids)
  if (any(!is.finite(n)) || any(n <= 0 | n != floor(n)) || !is.finite(sum(n)))
    stop("pool_sizes must contain positive finite integer sizes", call. = FALSE)
  C <- if (is.null(covariates)) NULL else
    .pg_align(.pg_matrix(covariates, "covariates"), pool_ids, what = "covariates")
  if (!is.null(C) && "(Intercept)" %in% colnames(C))
    stop("Do not supply an intercept column; the model includes one", call. = FALSE)
  missing_variant <- colSums(is.na(F)) > 0
  low_depth <- rep(FALSE, ncol(F))
  if (!is.null(counts)) {
    counts$alt <- .pg_align(counts$alt, pool_ids, variant_ids, "counts$alt")
    counts$depth <- .pg_align(counts$depth, pool_ids, variant_ids, "counts$depth")
    missing_variant <- missing_variant | colSums(is.na(counts$alt) | is.na(counts$depth) |
                                                  counts$depth == 0, na.rm = TRUE) > 0
    low_depth <- colSums(counts$depth < min_depth, na.rm = TRUE) > 0
    known <- is.finite(counts$alt) & is.finite(counts$depth) & counts$depth > 0 & is.finite(F)
    if (any(abs(F[known] - counts$alt[known] / counts$depth[known]) > 1e-10))
      stop("frequencies are inconsistent with counts", call. = FALSE)
    if (counts$type == "chromosomes" &&
        any(abs(sweep(counts$depth, 1L, 2 * n, "-")) > 0, na.rm = TRUE))
      stop("Chromosome depth must equal twice the diploid pool size", call. = FALSE)
  }
  if (any(missing_variant) && missing == "error")
    stop("Missing frequency/count or zero depth; use missing='drop_variants' to exclude affected variants", call. = FALSE)
  global_frequency <- colSums(F * (n / sum(n)))
  maf <- pmin(global_frequency, 1 - global_frequency)
  low_maf <- !is.na(maf) & maf < min_maf
  drop <- missing_variant | low_depth | low_maf
  reason <- vapply(seq_along(drop), function(j)
    paste(c(if (missing_variant[j]) "missing_or_zero_depth",
            if (low_depth[j]) "low_depth", if (low_maf[j]) "low_maf"), collapse = ";"), "")
  if (all(drop)) stop("No variants remain after quality control", call. = FALSE)
  if (is.null(variants)) {
    variants <- data.frame(variant_id = variant_ids, stringsAsFactors = FALSE)
  } else {
    if (!is.data.frame(variants) || !"variant_id" %in% names(variants))
      stop("variants must be a data frame with variant_id", call. = FALSE)
    .pg_ids(as.character(variants$variant_id), nrow(variants), "variants")
    if (!setequal(variants$variant_id, variant_ids)) stop("variants IDs do not match", call. = FALSE)
    variants <- variants[match(variant_ids, variants$variant_id), , drop = FALSE]
  }
  excluded <- data.frame(variant_id = variant_ids[drop], reason = reason[drop], stringsAsFactors = FALSE)
  F <- F[, !drop, drop = FALSE]
  variants <- variants[!drop, , drop = FALSE]
  rownames(variants) <- variants$variant_id
  if (!is.null(counts)) {
    counts$alt <- counts$alt[, !drop, drop = FALSE]
    counts$depth <- counts$depth[, !drop, drop = FALSE]
  }
  between <- vapply(seq_len(ncol(F)), function(j) diff(range(F[, j])), 0)
  constant_traits <- colnames(Y)[vapply(seq_len(ncol(Y)), function(j) diff(range(Y[, j])) == 0, FALSE)]
  qc <- list(excluded = excluded, global_frequency = global_frequency[!drop], maf = maf[!drop],
             no_between_pool_variation = colnames(F)[between == 0], constant_traits = constant_traits,
             n_input_variants = length(variant_ids), n_retained_variants = ncol(F),
             missing_policy = missing, min_maf = min_maf, min_depth = min_depth)
  structure(list(frequencies = F, phenotypes = Y, pool_sizes = n, counts = counts,
                 covariates = C, variants = variants, qc = qc, design = design,
                 provenance = provenance, schema_version = "2.0"), class = "pool_gwas_data")
}

validate_pool_data <- function(data) {
  if (!inherits(data, "pool_gwas_data")) stop("Expected a pool_gwas_data object", call. = FALSE)
  validated <- as_pool_gwas_data(data$frequencies, data$phenotypes, data$pool_sizes,
    counts = data$counts, covariates = data$covariates, variants = data$variants,
    design = data$design, provenance = data$provenance)
  # Recompute current measurements but retain the provenance of prior filtering
  # and ancillary individual-level summaries carried by a constructed object.
  for (name in setdiff(names(data), names(validated))) validated[[name]] <- data[[name]]
  if (is.list(data$qc)) {
    if (is.data.frame(data$qc$excluded) &&
        identical(names(data$qc$excluded), c("variant_id", "reason")))
      validated$qc$excluded <- data$qc$excluded
    for (name in c("n_input_variants", "missing_policy", "min_maf", "min_depth"))
      if (!is.null(data$qc[[name]])) validated$qc[[name]] <- data$qc[[name]]
  }
  list(valid = TRUE, data = validated, qc = validated$qc)
}

print.pool_gwas_data <- function(x, ...) {
  cat("PooGwaM data:", nrow(x$frequencies), "pools,", ncol(x$frequencies),
      "variants,", ncol(x$phenotypes), "traits;", sum(x$pool_sizes), "diploid individuals\n")
  cat("Measurement:", if (is.null(x$counts)) "frequencies only" else x$counts$type, "\n")
  cat("Excluded variants:", nrow(x$qc$excluded), "; constant traits:", length(x$qc$constant_traits), "\n")
  invisible(x)
}

pool_from_individuals <- function(genotypes, phenotypes, clusters, covariates = NULL,
                                  retain_individuals = FALSE, ...) {
  G <- .pg_matrix(genotypes, "genotypes")
  if (any(!G %in% 0:2)) stop("genotypes must be observed diploid 0/1/2 dosages", call. = FALSE)
  Y <- .pg_align(.pg_matrix(phenotypes, "phenotypes"), rownames(G), what = "phenotypes")
  design <- if (inherits(clusters, "pool_design")) clusters else list()
  if (inherits(clusters, "pool_design")) clusters <- clusters$clusters
  ids <- rownames(G)
  .pg_ids(names(clusters), length(clusters), "clusters individuals")
  if (!setequal(names(clusters), ids)) stop("clusters individual IDs do not match", call. = FALSE)
  groups <- as.character(clusters[ids])
  if (anyNA(groups) || any(!nzchar(groups))) stop("clusters contains missing pool IDs", call. = FALSE)
  pools <- sort(unique(groups))
  if (length(pools) < 2L) stop("At least two occupied pools are required", call. = FALSE)
  n <- setNames(as.numeric(table(factor(groups, levels = pools))), pools)
  aggregate_sum <- function(x) rowsum(x, group = groups, reorder = TRUE)[pools, , drop = FALSE]
  alt <- aggregate_sum(G)
  depth <- matrix(rep(2 * n, ncol(G)), nrow = length(pools), dimnames = dimnames(alt))
  means <- aggregate_sum(Y) / as.numeric(n)
  C <- if (is.null(covariates)) NULL else
    .pg_align(.pg_matrix(covariates, "covariates"), ids, what = "covariates")
  Cmeans <- if (is.null(C)) NULL else aggregate_sum(C) / as.numeric(n)
  if (!is.logical(retain_individuals) || length(retain_individuals) != 1L || is.na(retain_individuals))
    stop("retain_individuals must be TRUE or FALSE", call. = FALSE)
  out <- as_pool_gwas_data(alt / depth, means, n,
      counts = list(alt = alt, depth = depth, type = "chromosomes"),
      covariates = Cmeans, design = design, ...)
  out$within_pool_ss <- aggregate_sum((Y - means[groups, , drop = FALSE])^2)
  out$individual_trait_sd <- apply(Y, 2L, stats::sd)
  if (retain_individuals)
    out$individuals <- list(genotypes = G, phenotypes = Y,
                            clusters = setNames(groups, ids), covariates = C)
  out
}

compute_group_MAFs <- function(df, groups) {
  G <- .pg_matrix(df, "genotypes")
  if (any(!G %in% 0:2)) stop("genotypes must contain observed 0/1/2 dosages", call. = FALSE)
  .pg_ids(names(groups), length(groups), "groups individuals")
  if (!setequal(names(groups), rownames(G))) stop("groups IDs do not match genotypes", call. = FALSE)
  groups <- as.character(groups[rownames(G)])
  if (anyNA(groups) || any(!nzchar(groups))) stop("Missing pool IDs", call. = FALSE)
  pools <- sort(unique(groups))
  sizes <- as.numeric(table(factor(groups, levels = pools)))
  rowsum(G, groups, reorder = TRUE)[pools, , drop = FALSE] / (2 * sizes)
}
