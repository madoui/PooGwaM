# Experimental Model-X statistics on individual genotypes followed by pooling.
# Pairwise canonicalization makes the fitted input identical after any swaps,
# including the summation order used by the numerical solvers.

.knockoff_genotypes <- function(x, label) {
  x <- .pg_matrix(x, label)
  if (nrow(x) < 2L || any(!x %in% 0:2))
    stop(label, " requires at least two individuals with observed diploid 0/1/2 dosages.",
         call. = FALSE)
  x
}

.knockoff_rng <- function(seed) {
  # Use the same restoration policy as calibration, including rejection of
  # Box-Muller caches and external user-supplied RNG state before any mutation.
  restore <- .calibration_rng(seed)
  recorded <- attr(restore, "provenance")
  list(seed = recorded$seed, kind = recorded$kind,
       start = recorded$initial_state, restore = restore,
       r_version = recorded$r_version, package_version = recorded$package_version,
       source_mode = recorded$source_mode)
}

.knockoff_provenance <- function(x) {
  if (is.null(x)) return(list(source = "external_pair_without_generator_metadata"))
  if (!is.list(x) || (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
      any(!nzchar(names(x))) || anyDuplicated(names(x)))))
    stop("generator_provenance must be a uniquely named list.", call. = FALSE)
  x
}

.knockoff_alleles <- function(x, ids) {
  if (is.null(x)) return(NULL)
  if (!is.list(x) || !identical(sort(names(x)), c("genotypes", "knockoffs")))
    stop("allele_metadata requires genotypes and knockoffs data frames.", call. = FALSE)
  aligned <- lapply(names(x), function(label) {
    tab <- x[[label]]
    if (!is.data.frame(tab) || !all(c("variant_id", "ref", "alt") %in% names(tab)))
      stop("Each allele_metadata table requires variant_id, ref and alt.", call. = FALSE)
    tab <- tab[, c("variant_id", "ref", "alt"), drop = FALSE]
    tab[] <- lapply(tab, as.character)
    .pg_ids(tab$variant_id, nrow(tab), paste("allele_metadata", label))
    if (!setequal(tab$variant_id, ids) || anyNA(tab) ||
        any(!nzchar(tab$ref)) || any(!nzchar(tab$alt)) || any(tab$ref == tab$alt))
      stop("allele_metadata must contain matching variant IDs and distinct nonempty alleles.",
           call. = FALSE)
    tab <- tab[match(ids, tab$variant_id), , drop = FALSE]
    rownames(tab) <- ids
    tab
  })
  if (!identical(aligned[[1L]], aligned[[2L]]))
    stop("Original and knockoff alleles must have identical orientation; no implicit flips are made.",
         call. = FALSE)
  aligned[[1L]]
}

.knockoff_threshold <- function(W, q) {
  .pg_scalar(q, "q", 0, 1)
  if (q == 0) stop("q must be strictly positive.", call. = FALSE)
  if (!is.numeric(W) || !length(W) || any(!is.finite(W)))
    stop("W must contain finite numeric contrasts.", call. = FALSE)
  # Accumulate tied magnitudes together. This is O(p log p), rather than
  # evaluating the complete vector for each of p candidate thresholds.
  magnitudes <- sort(unique(abs(W[W != 0])), decreasing = TRUE)
  if (!length(magnitudes)) return(Inf)
  bins <- match(abs(W), magnitudes)
  positives <- cumsum(tabulate(bins[W > 0], nbins = length(magnitudes)))
  negatives <- cumsum(tabulate(bins[W < 0], nbins = length(magnitudes)))
  eligible <- which((1 + negatives) / pmax(1, positives) <= q)
  if (!length(eligible)) Inf else magnitudes[max(eligible)]
}

.knockoff_fit_args <- function(fit_args, solver) {
  if (!is.list(fit_args) || (length(fit_args) &&
      (is.null(names(fit_args)) || anyNA(names(fit_args)) ||
       any(!nzchar(names(fit_args))) || anyDuplicated(names(fit_args)))))
    stop("fit_args must be a uniquely named list of fixed solver settings.", call. = FALSE)
  if (any(!names(fit_args) %in% setdiff(names(formals(solver)), "data")))
    stop("Unknown or forbidden fit_args for the selected method.", call. = FALSE)
  fit_args
}

#' Generate independent reference knockoffs under a known product-binomial model
#' @export
generate_independent_knockoffs <- function(genotypes, frequencies, seed = 1) {
  G <- .knockoff_genotypes(genotypes, "genotypes")
  if (!is.numeric(frequencies) || !length(frequencies) ||
      any(!is.finite(frequencies)) || any(frequencies < 0 | frequencies > 1))
    stop("frequencies must contain externally known allele probabilities in [0,1].",
         call. = FALSE)
  .pg_ids(names(frequencies), length(frequencies), "frequencies")
  if (!setequal(names(frequencies), colnames(G)))
    stop("frequencies must be named by every genotype variant ID.", call. = FALSE)
  rows <- sort(rownames(G))
  variants <- sort(colnames(G))
  frequencies <- frequencies[variants]
  if (any(G[, variants[frequencies == 0], drop = FALSE] != 0) ||
      any(G[, variants[frequencies == 1], drop = FALSE] != 2))
    stop("Observed genotypes are outside the support of the supplied frequencies.", call. = FALSE)
  rng <- .knockoff_rng(seed)
  on.exit(rng$restore(), add = TRUE)
  copies <- matrix(stats::rbinom(length(G), size = 2,
    prob = rep(frequencies, each = nrow(G))), nrow(G), ncol(G),
    dimnames = list(rows, variants))
  copies <- copies[rownames(G), colnames(G), drop = FALSE]
  attr(copies, "knockoff_provenance") <- list(
    generator = "independent_binomial", frequencies = frequencies,
    frequency_source = "externally_supplied_as_known_not_estimated_by_this_function",
    model = "Independent individuals and variants with G_ij ~ Binomial(2, f_j).",
    certification = "The original genotype distribution and the frequency source are not certified.",
    rng = list(seed = rng$seed, kind = rng$kind, start = rng$start,
               end = get(".Random.seed", envir = .GlobalEnv)),
    R_version = rng$r_version, package_version = rng$package_version,
    source_mode = rng$source_mode,
    fingerprint_method = "MD5 of uncompressed R serialization version 2",
    input_fingerprint = .calibration_fingerprint(G[rows, variants, drop = FALSE]),
    settings_fingerprint = .calibration_fingerprint(list(frequencies = frequencies,
                                                         seed = rng$seed, kind = rng$kind)),
    output_fingerprint = .calibration_fingerprint(copies[rows, variants, drop = FALSE]))
  copies
}

#' Experimental knockoff+ selection after phenotype-defined pooling
#' @export
select_pool_knockoffs <- function(genotypes, knockoffs, phenotypes, pool_builder,
                                  q = 0.1, method = c("ridge", "multitrait"),
                                  fit_args = list(), generator_provenance = NULL,
                                  allele_metadata = NULL, seed = 1) {
  method <- match.arg(method)
  .pg_scalar(q, "q", 0, 1)
  if (q == 0) stop("q must be strictly positive.", call. = FALSE)
  inherited_provenance <- attr(knockoffs, "knockoff_provenance", exact = TRUE)
  provenance <- .knockoff_provenance(if (is.null(generator_provenance))
    inherited_provenance else generator_provenance)
  G <- .knockoff_genotypes(genotypes, "genotypes")
  K <- .knockoff_genotypes(knockoffs, "knockoffs")
  Y <- .pg_matrix(phenotypes, "phenotypes")
  variant_ids <- colnames(G)
  ids <- sort(variant_ids)
  individuals <- sort(rownames(G))
  if (!setequal(rownames(K), individuals) || !setequal(colnames(K), ids) ||
      !setequal(rownames(Y), individuals))
    stop("genotypes, knockoffs and phenotypes must have matching individual and variant IDs.",
         call. = FALSE)
  G <- G[individuals, ids, drop = FALSE]
  K <- K[individuals, ids, drop = FALSE]
  Y <- Y[individuals, , drop = FALSE]
  alleles <- .knockoff_alleles(allele_metadata, ids)
  if (!is.function(pool_builder))
    stop("pool_builder must be a phenotype-only function.", call. = FALSE)
  solver <- if (method == "ridge") fit_pool_ridge else fit_pool_multitrait
  fit_args <- .knockoff_fit_args(fit_args, solver)
  p <- ncol(G)
  # The first distinct dosage determines a stable unordered-pair orientation.
  # The outcome and original/knockoff labels never choose the fitted order.
  original_first <- identical_pair <- logical(p)
  for (j in seq_len(p)) {
    difference <- which(G[, j] != K[, j])
    identical_pair[j] <- !length(difference)
    original_first[j] <- identical_pair[j] || G[difference[1L], j] < K[difference[1L], j]
  }
  first <- 2 * seq_len(p) - 1
  second <- first + 1
  original_slot <- ifelse(original_first, first, second)
  knockoff_slot <- ifelse(original_first, second, first)
  internal_ids <- paste0("pair_column_", seq_len(2 * p))
  augmented <- matrix(0, nrow(G), 2 * p,
                      dimnames = list(individuals, internal_ids))
  augmented[, original_slot] <- G
  augmented[, knockoff_slot] <- K
  rng <- .knockoff_rng(seed)
  on.exit(rng$restore(), add = TRUE)
  design <- pool_builder(Y)
  after_design <- get(".Random.seed", envir = .GlobalEnv)
  data <- pool_from_individuals(augmented, Y, design)
  if (!is.null(data$covariates) || !identical(colnames(data$frequencies), internal_ids))
    stop("Knockoff pooling must preserve every pair and cannot add covariates.", call. = FALSE)
  fit <- do.call(solver, c(list(data = data), fit_args))
  if (!isTRUE(fit$diagnostics$converged))
    stop("Knockoff selection requires a converged fit; no selection was computed.", call. = FALSE)
  B <- fit$beta_standardized
  largest <- apply(abs(B), 1L, max)
  divisor <- largest
  divisor[divisor == 0] <- 1
  importance <- largest * sqrt(rowSums((B / divisor)^2))
  if (any(!is.finite(importance)))
    stop("Knockoff importance scores are not finite; rescale or regularize.", call. = FALSE)
  W <- importance[original_slot] - importance[knockoff_slot]
  # Identical copies are exchangeable even at the same input: their W must be
  # exactly zero, regardless of roundoff across coefficient blocks.
  W[identical_pair] <- 0
  names(W) <- ids
  W <- W[variant_ids]
  threshold <- .knockoff_threshold(W, q)
  order_out <- match(variant_ids, ids)
  original_beta <- B[original_slot[order_out], , drop = FALSE]
  knockoff_beta <- B[knockoff_slot[order_out], , drop = FALSE]
  rownames(original_beta) <- rownames(knockoff_beta) <- variant_ids
  mapping <- data.frame(variant_id = variant_ids,
    original_column = internal_ids[original_slot[order_out]],
    knockoff_column = internal_ids[knockoff_slot[order_out]], stringsAsFactors = FALSE)
  structure(list(
    W = W, threshold = threshold, selected = names(W)[W >= threshold],
    q = q, filter = "knockoff+", method = method,
    importance = data.frame(variant_id = variant_ids,
      original = importance[original_slot[order_out]],
      knockoff = importance[knockoff_slot[order_out]], row.names = NULL),
    beta_standardized = original_beta, knockoff_beta_standardized = knockoff_beta,
    fit = fit, pair_map = mapping, pool_design = design,
    inference_status = "experimental_knockoff",
    inference_scope = "conditional_variant_null",
    inference_unit = "variant_association_with_the_trait_vector",
    generator_provenance = provenance,
    model_validation = "not_certified",
    allele_metadata = if (is.null(alleles)) NULL else alleles[variant_ids, , drop = FALSE],
    allele_validation = if (is.null(alleles)) "shared_orientation_assumed_not_checked" else
      "supplied_pair_labels_match_not_independently_verified",
    diagnostics = list(identical_pairs = ids[identical_pair],
      qc = "Complete exact dosages required; no variant or pair is filtered.",
      pair_order = "Lexicographic individual dosages within each unordered pair.",
      solver = fit$diagnostics),
    rng = list(seed = rng$seed, kind = rng$kind, start = rng$start,
               after_pool_builder = after_design,
               end = get(".Random.seed", envir = .GlobalEnv)),
    R_version = rng$r_version, package_version = rng$package_version,
    source_mode = rng$source_mode,
    fingerprint_method = "MD5 of uncompressed R serialization version 2",
    input_fingerprint = .calibration_fingerprint(list(genotypes = G, knockoffs = K,
                                                      phenotypes = Y, alleles = alleles)),
    settings_fingerprint = .calibration_fingerprint(list(design = design,
      method = method, fit_args = fit$settings, q = q, seed = rng$seed, kind = rng$kind)),
    assumptions = c(
      "Individual original/knockoff pairs must satisfy Model-X exchangeability and knockoffs independent of phenotypes conditional on genotypes.",
      "The initial domain is independent individuals, exact diploid dosages and a known genotype model.",
      "Pooling uses phenotypes only, with independent randomness if any, and identical weights for each pair.",
      "Solver settings are fixed independently of original/knockoff labels; no outcome-based asymmetric tuning is allowed.",
      "Allele orientation must agree across each original/knockoff pair."),
    limitations = c(
      "Matrix validation and declared provenance cannot certify the original genotype model or exchangeability.",
      "Conditional-variant FDR follows only under the stated Model-X and sign-flip assumptions; no unconditional guarantee is claimed for external pairs.",
      "No adjustment for covariates, ancestry, kinship, estimated LD or sequencing reads is supplied.",
      "Selection concerns the full trait vector, not separate variant-trait hypotheses or causal effects.",
      "Few pools and highly correlated copies may greatly reduce power.")),
    class = "poogwam_knockoff")
}
