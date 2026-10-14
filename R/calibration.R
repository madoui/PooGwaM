# Calibration deliberately never estimates a null distribution across SNPs.
# A randomization test here targets global genotype/phenotype independence.

.calibration_integer <- function(x, name, minimum = 1) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max) {
    stop(name, " must be one integer between ", minimum,
         " and .Machine$integer.max.", call. = FALSE)
  }
  as.integer(x)
}

.calibration_rng <- function(seed, rng_kind = NULL) {
  seed <- .calibration_integer(seed, "seed", 0)
  if (!is.null(rng_kind) && (!is.character(rng_kind) ||
      length(rng_kind) != 3L || anyNA(rng_kind) || any(!nzchar(rng_kind))))
    stop("rng_kind must be NULL or the three strings returned by RNGkind().", call. = FALSE)
  kind <- RNGkind()
  if (kind[2L] == "Box-Muller" || (!is.null(rng_kind) && rng_kind[2L] == "Box-Muller"))
    stop("Box-Muller has a normal cache that R cannot restore from .Random.seed; ",
         "select normal.kind = 'Inversion' before calibration.", call. = FALSE)
  if (any(kind == "user-supplied") || any(rng_kind == "user-supplied"))
    stop("user-supplied RNG state cannot be captured and restored by calibration.", call. = FALSE)
  present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  previous <- if (present) get(".Random.seed", envir = .GlobalEnv) else NULL
  restore <- function() {
    do.call(RNGkind, as.list(kind))
    if (present) assign(".Random.seed", previous, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }
  tryCatch({
    if (!is.null(rng_kind)) do.call(RNGkind, as.list(unname(rng_kind)))
    set.seed(seed)
  }, error = function(e) {
    restore()
    stop(e)
  })
  attr(restore, "provenance") <- list(
    seed = seed, kind = setNames(RNGkind(), c("kind", "normal.kind", "sample.kind")),
    initial_state = get(".Random.seed", envir = .GlobalEnv),
    r_version = R.version.string,
    package_version = if (isNamespace(environment(.calibration_rng)))
      as.character(getNamespaceVersion(environment(.calibration_rng))) else NA_character_,
    source_mode = if (isNamespace(environment(.calibration_rng))) "namespace" else "sourced")
  restore
}

.calibration_rng_state <- function() get(".Random.seed", envir = .GlobalEnv)

# Hash a streamed serialization rather than allocating another genome-sized
# raw vector. This is a provenance checksum, not a credential or signature.
.calibration_fingerprint <- function(object) {
  path <- tempfile("poogwam-provenance-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, compress = FALSE, version = 2)
  unname(tools::md5sum(path))
}

.calibration_fit_args <- function(fit_args) {
  if (!is.list(fit_args) || (length(fit_args) &&
      (is.null(names(fit_args)) || anyNA(names(fit_args)) ||
       any(!nzchar(names(fit_args))) || anyDuplicated(names(fit_args))))) {
    stop("fit_args must be a uniquely named list of fixed solver settings.",
         call. = FALSE)
  }
  allowed <- setdiff(names(formals(fit_pool_ridge)), "data")
  if (any(!names(fit_args) %in% allowed)) {
    stop("Unknown or forbidden fit_args: ",
         paste(setdiff(names(fit_args), allowed), collapse = ", "), call. = FALSE)
  }
  fit_args
}

.calibration_matrix <- function(x, label) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x)) ||
      nrow(x) < 2L || ncol(x) < 1L)
    stop(label, " must be a finite numeric matrix with at least two rows.",
         call. = FALSE)
  for (ids in list(rownames(x), colnames(x))) {
    if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
      stop(label, " must have unique, nonempty row and column identifiers.",
           call. = FALSE)
  }
  x
}

.calibration_scores <- function(fit) {
  if (!inherits(fit, "poogwam_fit") || !is.matrix(fit$beta) ||
      !is.matrix(fit$beta_standardized) ||
      !identical(dimnames(fit$beta), dimnames(fit$beta_standardized)) ||
      !identical(dim(fit$beta), dim(fit$beta_standardized)) ||
      any(!is.finite(fit$beta)) || any(!is.finite(fit$beta_standardized))) {
    stop("fit must be a finite poogwam_fit with standardized effects.",
         call. = FALSE)
  }
  absolute <- abs(fit$beta_standardized)
  largest <- apply(absolute, 1L, max)
  divisor <- largest
  divisor[divisor == 0] <- 1
  # Scale before squaring to avoid intermediate overflow in Euclidean norms.
  omnibus <- largest * sqrt(rowSums((absolute / divisor)^2))
  if (any(!is.finite(omnibus)))
    stop("Standardized omnibus scores overflow; rescale or regularize the fit.", call. = FALSE)
  list(trait = abs(fit$beta), omnibus = omnibus)
}

.calibration_pilot <- function(fit, simulate_fit, B_pilot, scale_floor) {
  mean <- m2 <- setNames(numeric(nrow(fit$beta)), rownames(fit$beta))
  amplitude <- mean
  for (b in seq_len(B_pilot)) {
    null_fit <- simulate_fit(b)
    if (!identical(dimnames(null_fit$beta), dimnames(fit$beta)))
      stop("The pilot changed variant or trait identifiers.", call. = FALSE)
    scores <- .calibration_scores(null_fit)$omnibus
    # Welford moments on a dynamically scaled score avoid squaring very
    # large or tiny finite numbers. Scaling is separate for each variant.
    new_amplitude <- pmax(amplitude, scores)
    divisor <- new_amplitude
    divisor[divisor == 0] <- 1
    ratio <- amplitude / divisor
    mean <- mean * ratio
    m2 <- m2 * ratio^2
    amplitude <- new_amplitude
    normalized <- scores / divisor
    delta <- normalized - mean
    mean <- mean + delta / b
    m2 <- m2 + delta * (normalized - mean)
  }
  mean <- mean * amplitude
  sd <- sqrt(pmax(m2, 0) / (B_pilot - 1)) * amplitude
  if (any(!is.finite(mean)) || any(!is.finite(sd)))
    stop("Pilot moments overflow; rescale or regularize the fit.", call. = FALSE)
  reference <- max(c(mean, sd))
  status <- "pilot_scaled"
  if (reference == 0) {
    # A zero pilot can miss rare events. This fixed fallback is a valid
    # increasing transformation, not a claim that the score is constant.
    scale <- setNames(rep(1, length(mean)), names(mean))
    floor <- NA_real_
    status <- "all_zero_pilot_unit_scale"
    warning("All pilot omnibus scores are zero; using a fixed unit scale. ",
            "This does not establish zero null variance; inspect the pilot or increase B_pilot.",
            call. = FALSE)
  } else {
    floor <- max(scale_floor * reference, .Machine$double.xmin)
    scale <- pmax(sd, floor)
  }
  list(B = B_pilot, mean_omnibus = mean, sd_omnibus = sd, scale = scale,
       scale_floor = scale_floor, reference_scale = reference,
       effective_floor = floor, floored_variants = names(sd)[sd < scale],
       zero_sd_variants = names(sd)[sd == 0], status = status)
}

.calibration_studentize <- function(score, studentization) {
  result <- (score - studentization$mean_omnibus) / studentization$scale
  if (any(!is.finite(result)))
    stop("Studentized omnibus scores overflow; inspect pilot scales and regularization.",
         call. = FALSE)
  result
}

# A streaming implementation avoids retaining B full genome-wide fits.
.calibration_run <- function(fit, simulate_fit, B, method, metadata,
                             studentization = NULL) {
  observed <- .calibration_scores(fit)
  exceed <- matrix(0, nrow(fit$beta), ncol(fit$beta), dimnames = dimnames(fit$beta))
  exceed_omnibus <- exceed_max <- setNames(numeric(nrow(fit$beta)), rownames(fit$beta))
  mean_trait <- m2_trait <- exceed
  mean_omnibus <- m2_omnibus <- exceed_omnibus
  maxima <- numeric(B)
  maximum_variants <- character(B)
  maximum_ties <- integer(B)
  # Relative tolerance preserves p-values when a trait's measurement unit
  # changes, even to very small units. A zero observed score always has p = 1.
  tolerance_trait <- 64 * .Machine$double.eps * observed$trait
  tolerance_omnibus <- 64 * .Machine$double.eps * observed$omnibus
  if (!is.null(studentization)) {
    observed_t <- .calibration_studentize(observed$omnibus, studentization)
    # Account for cancellation when centering; unlike raw norms, t can be
    # negative. The scale is fixed using the independent pilot only.
    epsilon <- 64 * .Machine$double.eps
    scaled_roundoff <- function(value) {
      ratio <- value / studentization$scale
      out <- epsilon * ratio
      overflow <- !is.finite(ratio)
      # A centered score may be finite even when either uncentered ratio
      # overflows. Multiply first only in that case, avoiding both overflow
      # for large values and premature underflow for very small values.
      out[overflow] <- (epsilon * value[overflow]) / studentization$scale[overflow]
      out
    }
    tolerance_t <- pmax(epsilon * abs(observed_t),
      scaled_roundoff(observed$omnibus), scaled_roundoff(studentization$mean_omnibus))
    if (any(!is.finite(tolerance_t)))
      stop("Studentized comparison tolerance overflows; increase scale_floor or regularization.",
           call. = FALSE)
    maxima_t <- numeric(B)
    maximum_t_variants <- character(B)
    maximum_t_ties <- integer(B)
  }
  for (b in seq_len(B)) {
    null_fit <- simulate_fit(b)
    null <- .calibration_scores(null_fit)
    if (!identical(dimnames(null_fit$beta), dimnames(fit$beta))) {
      stop("The null procedure changed variant or trait identifiers.", call. = FALSE)
    }
    exceed <- exceed + (null$trait >= observed$trait - tolerance_trait)
    exceed_omnibus <- exceed_omnibus +
      (null$omnibus >= observed$omnibus - tolerance_omnibus)
    maxima[b] <- max(null$omnibus)
    maximum_variants[b] <- names(null$omnibus)[which.max(null$omnibus)]
    maximum_ties[b] <- sum(null$omnibus == maxima[b])
    if (is.null(studentization)) {
      exceed_max <- exceed_max + (maxima[b] >= observed$omnibus - tolerance_omnibus)
    } else {
      null_t <- .calibration_studentize(null$omnibus, studentization)
      maxima_t[b] <- max(null_t)
      maximum_t_variants[b] <- names(null_t)[which.max(null_t)]
      maximum_t_ties[b] <- sum(null_t == maxima_t[b])
      exceed_max <- exceed_max + (maxima_t[b] >= observed_t - tolerance_t)
    }
    delta <- null$trait - mean_trait
    mean_trait <- mean_trait + delta / b
    m2_trait <- m2_trait + delta * (null$trait - mean_trait)
    delta <- null$omnibus - mean_omnibus
    mean_omnibus <- mean_omnibus + delta / b
    m2_omnibus <- m2_omnibus + delta * (null$omnibus - mean_omnibus)
  }
  if (!is.null(studentization)) {
    studentization$observed_omnibus <- observed_t
    studentization$max_omnibus <- maxima_t
    studentization$comparison_tolerance <- tolerance_t
  }
  structure(c(list(
    method = method, fit = fit, B = B, resolution = 1 / (B + 1),
    maxT_method = if (is.null(studentization)) "raw" else "studentized",
    studentization = studentization,
    maxT_diagnostics = list(raw_max_variant = maximum_variants, raw_max_ties = maximum_ties,
      studentized_max_variant = if (is.null(studentization)) NULL else maximum_t_variants,
      studentized_max_ties = if (is.null(studentization)) NULL else maximum_t_ties),
    p_raw = (1 + exceed) / (B + 1),
    p_omnibus = (1 + exceed_omnibus) / (B + 1),
    p_maxT_omnibus = (1 + exceed_max) / (B + 1),
    inference_scope = "global_null_only",
    inference_status = if (method == "individual_permutation")
      "randomization_global_null" else "experimental_parametric_global_null",
    null_summary = list(
      mean_abs_effect = mean_trait,
      sd_abs_effect = if (B > 1L) sqrt(pmax(m2_trait, 0) / (B - 1)) else mean_trait * NA_real_,
      mean_omnibus = mean_omnibus,
      sd_omnibus = if (B > 1L) sqrt(pmax(m2_omnibus, 0) / (B - 1)) else mean_omnibus * NA_real_,
      max_omnibus = maxima)), metadata), class = "poogwam_calibration")
}

#' Calibrate ridge scores by individual genotype permutation
#'
#' See the installed help page for exchangeability and global-null restrictions.
calibrate_pool_permutation <- function(genotypes, phenotypes, pool_builder,
                                       B = 999, seed = 1, strata = NULL,
                                       fit_args = list(),
                                       maxT = c("raw", "studentized"),
                                       B_pilot = 199, scale_floor = 1e-8,
                                       rng_kind = NULL) {
  B <- .calibration_integer(B, "B")
  maxT <- match.arg(maxT)
  if (maxT == "studentized") {
    B_pilot <- .calibration_integer(B_pilot, "B_pilot", 2)
    .fit_scalar(scale_floor, "scale_floor", strict = TRUE)
    if (scale_floor > 1)
      stop("scale_floor must be in (0, 1].", call. = FALSE)
  }
  fit_args <- .calibration_fit_args(fit_args)
  genotypes <- .calibration_matrix(genotypes, "genotypes")
  phenotypes <- .calibration_matrix(phenotypes, "phenotypes")
  if (any(!genotypes %in% 0:2))
    stop("genotypes must contain observed diploid 0/1/2 dosages.", call. = FALSE)
  if (!setequal(rownames(genotypes), rownames(phenotypes)))
    stop("genotypes and phenotypes must contain identical individual IDs.", call. = FALSE)
  genotypes <- genotypes[rownames(phenotypes), , drop = FALSE]
  if (!is.function(pool_builder))
    stop("pool_builder must be a function of the complete phenotype matrix.", call. = FALSE)
  ids <- rownames(phenotypes)
  if (is.null(strata)) {
    blocks <- list(seq_along(ids))
  } else {
    if (length(strata) != length(ids) || !is.atomic(strata) ||
        is.null(names(strata)) || anyNA(names(strata)) ||
        anyDuplicated(names(strata)) || !setequal(names(strata), ids) ||
        anyNA(strata))
      stop("strata must be a nonmissing vector named by every individual ID.", call. = FALSE)
    strata <- strata[ids]
    blocks <- split(seq_along(ids), as.character(strata))
    if (!length(blocks) || !any(lengths(blocks) > 1L))
      stop("strata permit no nontrivial permutations.", call. = FALSE)
    warning("Blocked permutation assumes genotype rows are exchangeable within strata; ",
            "blocks alone do not establish this assumption.", call. = FALSE)
  }
  restore <- .calibration_rng(seed, rng_kind)
  on.exit(restore(), add = TRUE)
  rng <- attr(restore, "provenance")
  # Y stays fixed, preserving complete multivariate phenotype vectors and their
  # correlations. Conditional on Y, its phenotype-only pooling plan stays fixed.
  design <- pool_builder(phenotypes)
  rng$after_pooling <- .calibration_rng_state()
  dat <- pool_from_individuals(genotypes, phenotypes, design)
  if (!is.null(dat$covariates) && ncol(as.matrix(dat$covariates)) > 0L)
    stop("Permutation calibration does not support covariate adjustment.", call. = FALSE)
  observed_fit <- do.call(fit_pool_ridge, c(list(data = dat), fit_args))
  simulate <- function(b) {
    index <- seq_along(ids)
    for (block in blocks) index[block] <- block[sample.int(length(block))]
    permuted <- genotypes[index, , drop = FALSE]
    rownames(permuted) <- ids
    null_data <- pool_from_individuals(permuted, phenotypes, design)
    do.call(fit_pool_ridge, c(list(data = null_data), fit_args))
  }
  studentization <- NULL
  rng$pilot_start <- rng$pilot_end <- NULL
  if (maxT == "studentized") {
    rng$pilot_start <- .calibration_rng_state()
    studentization <- .calibration_pilot(observed_fit, simulate, B_pilot, scale_floor)
    rng$pilot_end <- .calibration_rng_state()
  }
  rng$main_start <- .calibration_rng_state()
  result <- .calibration_run(observed_fit, simulate, B, "individual_permutation", list(
    seed = seed, strata = strata,
    assumptions = c("Genotype rows exchangeable relative to fixed phenotype vectors within strata.",
                    "The pooling rule uses phenotypes only; solver settings are fixed in advance.",
                    "Exact individual dosages: no sequencing-read observation model is used."),
    limitations = c("Tests the global null of genotype-phenotype independence.",
                    "Per-variant causal nulls, strong FWER and FDR under partial alternatives are not guaranteed.",
                    "No automatic handling of ancestry, kinship, continuous covariates or outcome-dependent tuning.")),
    studentization = studentization)
  rng$main_end <- .calibration_rng_state()
  rng$strategy <- if (maxT == "raw") "single_series" else "successive_pilot_and_main_segments"
  result$rng <- rng
  result$provenance <- list(
    fingerprint_method = "MD5 of uncompressed R serialization version 2",
    input_fingerprint = .calibration_fingerprint(list(genotypes = genotypes,
      phenotypes = phenotypes, strata = strata)),
    settings_fingerprint = .calibration_fingerprint(list(design = design,
      fit = observed_fit$settings, B = B, maxT = maxT,
      B_pilot = if (maxT == "studentized") B_pilot else NULL,
      scale_floor = if (maxT == "studentized") scale_floor else NULL)),
    pool_design = design,
    pool_builder_expression = list(formals = formals(pool_builder), body = body(pool_builder)))
  result
}

#' Experimental bootstrap under a homogeneous pool-count null
calibrate_pool_counts <- function(data, B = 999, seed = 1, fit_args = list(),
                                  rng_kind = NULL) {
  B <- .calibration_integer(B, "B")
  fit_args <- .calibration_fit_args(fit_args)
  if (!inherits(data, "pool_gwas_data") || is.null(data$counts))
    stop("data must be pool_gwas_data with observed, typed counts.", call. = FALSE)
  if (!is.null(data$covariates) && ncol(as.matrix(data$covariates)) > 0L)
    stop("Count bootstrap does not support covariates or population structure.", call. = FALSE)
  counts <- data$counts
  if (length(counts$type) != 1L || !counts$type %in% c("reads", "chromosomes"))
    stop("counts$type must explicitly be 'reads' or 'chromosomes'.", call. = FALSE)
  # Recheck the full public data contract: S3 class membership alone does not
  # establish validity after callers have edited list components. The solver
  # retains a validated, canonically aligned object for compatibility checks.
  validate_pool_data(data)
  if (!is.matrix(counts$alt) || !is.matrix(counts$depth) ||
      !identical(dimnames(counts$alt), dimnames(data$frequencies)) ||
      !identical(dimnames(counts$depth), dimnames(data$frequencies)) ||
      !identical(dim(counts$alt), dim(data$frequencies)) ||
      !identical(dim(counts$depth), dim(data$frequencies)) ||
      any(!is.finite(counts$depth)) || any(counts$depth <= 0) ||
      any(counts$depth != floor(counts$depth)) ||
      any(!is.finite(counts$alt)) || any(counts$alt < 0) ||
      any(counts$alt > counts$depth) || any(counts$alt != floor(counts$alt)))
    stop("Counts require aligned integer alt/depth matrices with positive depths.", call. = FALSE)
  if (!isTRUE(all.equal(unname(data$frequencies),
                       unname(counts$alt / counts$depth), tolerance = 1e-12)))
    stop("frequencies must equal observed alt/depth for count calibration.", call. = FALSE)
  n <- as.numeric(data$pool_sizes[rownames(data$frequencies)])
  if (any(!is.finite(n)) || any(n < 1) || any(n != floor(n)))
    stop("Count bootstrap requires positive integer diploid pool sizes.", call. = FALSE)
  if (counts$type == "chromosomes" && any(counts$depth != 2 * n))
    stop("Chromosome count depths must equal twice each diploid pool size.", call. = FALSE)
  warning("Experimental count bootstrap assumes unrelated diploid individuals, homogeneous ",
          "allele frequencies, independent variants and equal DNA contributions; ",
          "plug-in p-values have no established finite-sample calibration.", call. = FALSE)
  restore <- .calibration_rng(seed, rng_kind)
  on.exit(restore(), add = TRUE)
  rng <- attr(restore, "provenance")
  observed_fit <- do.call(fit_pool_ridge, c(list(data = data), fit_args))
  frequencies <- as.numeric(crossprod(n / sum(n), data$frequencies))
  names(frequencies) <- colnames(data$frequencies)
  dimensions <- dim(data$frequencies)
  dimension_names <- dimnames(data$frequencies)
  simulate <- function(b) {
    latent <- matrix(stats::rbinom(prod(dimensions), size = rep(2 * n, dimensions[2]),
                                  prob = rep(frequencies, each = dimensions[1])),
                     dimensions[1], dimensions[2], dimnames = dimension_names)
    alt <- if (counts$type == "reads") {
      matrix(stats::rbinom(length(latent), size = as.vector(counts$depth),
                          prob = as.vector(latent / (2 * n))),
             dimensions[1], dimensions[2], dimnames = dimension_names)
    } else latent
    null_data <- as_pool_gwas_data(
      frequencies = alt / counts$depth,
      phenotypes = data$phenotypes,
      pool_sizes = data$pool_sizes,
      counts = list(alt = alt, depth = counts$depth, type = counts$type))
    do.call(fit_pool_ridge, c(list(data = null_data), fit_args))
  }
  rng$main_start <- .calibration_rng_state()
  result <- .calibration_run(observed_fit, simulate, B, "pool_count_bootstrap", list(
    seed = seed, null_frequencies = frequencies, count_type = counts$type,
    assumptions = c("One shared allele frequency per variant across every pool.",
                    "Independent unrelated diploid chromosomes and independent variants.",
                    "Read counts, when present, are binomial draws from latent chromosome fractions.",
                    "Observed pool design, phenotype means, sizes and sequencing depths are fixed."),
    limitations = c("Experimental parametric bootstrap with estimated null frequencies; +1 is not proof of calibration.",
                    "No LD, kinship, ancestry, unequal DNA contributions, sequencing errors or overdispersion model.",
                    "No established FDR, strong FWER or partial-null guarantees.")))
  rng$main_end <- .calibration_rng_state()
  rng$strategy <- "single_series"
  result$rng <- rng
  result$provenance <- list(
    fingerprint_method = "MD5 of uncompressed R serialization version 2",
    input_fingerprint = .calibration_fingerprint(data),
    settings_fingerprint = .calibration_fingerprint(list(fit = observed_fit$settings, B = B)))
  result
}

#' Return effects and explicitly labelled association scores
association_table <- function(fit, calibration = NULL,
                              adjustment = c("BH", "BY", "holm", "none")) {
  adjustment <- match.arg(adjustment)
  scores <- .calibration_scores(fit)
  p <- nrow(fit$beta)
  d <- ncol(fit$beta)
  status <- "exploratory"
  scope <- "none"
  raw <- matrix(NA_real_, p, d, dimnames = dimnames(fit$beta))
  omnibus_raw <- max_t <- rep(NA_real_, p)
  resolution <- NA_real_
  max_t_method <- NA_character_
  if (!is.null(calibration)) {
    if (!inherits(calibration, "poogwam_calibration") ||
        !identical(fit$data, calibration$fit$data) ||
        !identical(fit$settings, calibration$fit$settings) ||
        !identical(fit$beta, calibration$fit$beta) ||
        !identical(fit$beta_standardized, calibration$fit$beta_standardized))
      stop("calibration must match this fit's exact data, settings and coefficients.", call. = FALSE)
    if (!identical(dimnames(calibration$p_raw), dimnames(fit$beta)) ||
        !identical(dim(calibration$p_raw), dim(fit$beta)) ||
        !identical(names(calibration$p_omnibus), rownames(fit$beta)) ||
        !identical(names(calibration$p_maxT_omnibus), rownames(fit$beta)))
      stop("Calibration p-values have incompatible variant or trait identifiers.", call. = FALSE)
    all_p <- c(calibration$p_raw, calibration$p_omnibus, calibration$p_maxT_omnibus)
    if (any(!is.finite(all_p)) || any(all_p <= 0 | all_p > 1))
      stop("Calibration p-values must be finite and in (0, 1].", call. = FALSE)
    raw <- calibration$p_raw
    omnibus_raw <- calibration$p_omnibus
    max_t <- calibration$p_maxT_omnibus
    status <- calibration$inference_status
    scope <- calibration$inference_scope
    resolution <- calibration$resolution
    max_t_method <- if (is.null(calibration$maxT_method)) "raw" else calibration$maxT_method
    if (!is.character(max_t_method) || length(max_t_method) != 1L ||
        is.na(max_t_method) || !max_t_method %in% c("raw", "studentized"))
      stop("Calibration maxT_method must be 'raw' or 'studentized'.", call. = FALSE)
  }
  # These are two explicitly distinct families, each adjusted once from raw p.
  adjust <- function(values) stats::p.adjust(values, method = adjustment)
  trait_table <- data.frame(
    variant_id = rep(rownames(fit$beta), times = d),
    trait = rep(colnames(fit$beta), each = p),
    test = "trait", effect = as.vector(fit$beta),
    effect_standardized = as.vector(fit$beta_standardized),
    score = as.vector(scores$trait), p_raw = as.vector(raw),
    p_adjusted = adjust(as.vector(raw)), p_maxT = NA_real_,
    test_family = "all_variant_trait_pairs", stringsAsFactors = FALSE)
  omnibus_table <- data.frame(
    variant_id = rownames(fit$beta), trait = NA_character_, test = "omnibus",
    effect = NA_real_, effect_standardized = NA_real_, score = scores$omnibus,
    p_raw = omnibus_raw, p_adjusted = adjust(omnibus_raw), p_maxT = max_t,
    test_family = "all_variant_omnibus_scores", stringsAsFactors = FALSE,
    row.names = NULL)
  result <- rbind(trait_table, omnibus_table)
  result$adjustment <- adjustment
  result$inference_status <- status
  result$inference_scope <- scope
  result$p_resolution <- resolution
  result$maxT_method <- ifelse(result$test == "omnibus", max_t_method, NA_character_)
  attr(result, "inference_note") <- paste(
    "Adjustment is arithmetic only: validity requires calibrated null p-values and",
    "the chosen method's dependence assumptions. Global-null calibration alone does",
    "not justify causal-variant FDR or strong FWER under partial alternatives.")
  result
}
