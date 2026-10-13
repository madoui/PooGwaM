# Deterministic pooled regression. Frequencies remain on the input allele scale;
# genotype dosage is always 2F. No individual genotype matrix is reconstructed.

.fit_scalar <- function(x, name, lower = 0, strict = FALSE, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      (if (strict) x <= lower else x < lower) ||
      (integer && (x != floor(x) || x > .Machine$integer.max)))
    stop(name, " has an invalid value", call. = FALSE)
  invisible(x)
}

.fit_ids <- function(x, name) {
  if (is.null(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x))
    stop(name, " must be nonempty unique identifiers", call. = FALSE)
  invisible(x)
}

.fit_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x)) ||
      nrow(x) < 1L || ncol(x) < 1L)
    stop(name, " must be a finite numeric matrix", call. = FALSE)
  .fit_ids(rownames(x), paste(name, "row names"))
  .fit_ids(colnames(x), paste(name, "column names"))
  invisible(x)
}

.fit_blocks <- function(p, size) {
  starts <- seq.int(1L, p, by = size)
  lapply(starts, function(first) seq.int(first, min(p, first + size - 1L)))
}

.fit_prepare <- function(data, penalty, penalty_floor, block_size,
                         standardize_traits) {
  if (!inherits(data, "pool_gwas_data"))
    stop("data must inherit from pool_gwas_data", call. = FALSE)
  .fit_scalar(penalty_floor, "penalty_floor", strict = TRUE)
  .fit_scalar(block_size, "block_size", lower = 1, integer = TRUE)
  if (!is.logical(standardize_traits) || length(standardize_traits) != 1L ||
      is.na(standardize_traits))
    stop("standardize_traits must be TRUE or FALSE", call. = FALSE)
  F <- data$frequencies
  Y <- data$phenotypes
  .fit_matrix(F, "frequencies")
  .fit_matrix(Y, "phenotypes")
  if (nrow(F) < 2L || any(F < 0 | F > 1))
    stop("frequencies require at least two pools and values in [0, 1]", call. = FALSE)
  ids <- rownames(F)
  if (!setequal(ids, rownames(Y)))
    stop("phenotypes and frequencies have different pool IDs", call. = FALSE)
  Y <- Y[ids, , drop = FALSE]
  sizes <- data$pool_sizes
  .fit_ids(names(sizes), "pool_sizes names")
  if (!is.numeric(sizes) || !setequal(names(sizes), ids) ||
      any(!is.finite(sizes)) || any(sizes <= 0 | sizes != floor(sizes)))
    stop("pool_sizes must be positive integers with matching pool IDs", call. = FALSE)
  sizes <- sizes[ids]
  # Scaling first avoids an overflowing sum for large finite sizes.
  weights <- sizes / max(sizes)
  weights <- weights / sum(weights)
  sqrt_weights <- sqrt(weights)
  C <- data$covariates
  if (!is.null(C)) {
    .fit_matrix(C, "covariates")
    if (!setequal(rownames(C), ids))
      stop("covariates have different pool IDs", call. = FALSE)
    if ("(Intercept)" %in% colnames(C))
      stop("covariates cannot be named (Intercept)", call. = FALSE)
    C <- C[ids, , drop = FALSE]
  }
  Z <- cbind("(Intercept)" = rep(1, nrow(F)), C)
  rownames(Z) <- ids
  Zw <- Z * sqrt_weights
  zqr <- qr(Zw, tol = 1e-10)
  if (zqr$rank != ncol(Z))
    stop("intercept and covariates must have full column rank", call. = FALSE)
  if (ncol(Z) >= nrow(F))
    stop("at least one residual pool degree of freedom is required", call. = FALSE)
  Q <- qr.Q(zqr)
  y_mean <- colSums(Y * weights)
  y_centered <- sweep(Y, 2L, y_mean, "-")
  is_constant <- vapply(seq_len(ncol(Y)), function(j) all(Y[, j] == Y[1L, j]), FALSE)
  y_centered[, is_constant] <- 0
  # This rescaling avoids squaring very large or very small response units.
  amplitude <- apply(abs(y_centered), 2L, max)
  denominator <- amplitude
  denominator[denominator == 0] <- 1
  trait_sd <- sqrt(colSums(sweep(y_centered, 2L, denominator, "/")^2 * weights)) * amplitude
  if (any(!is.finite(trait_sd)))
    stop("phenotype scale overflows; rescale phenotypes before fitting", call. = FALSE)
  constant_traits <- colnames(Y)[is_constant]
  if (any(trait_sd == 0 & !is_constant))
    stop("phenotype scale underflows; rescale phenotypes before fitting", call. = FALSE)
  trait_sd[trait_sd == 0] <- 1
  working_scales <- if (standardize_traits) trait_sd else rep(1, ncol(Y))
  names(working_scales) <- colnames(Y)
  Ys <- sweep(Y, 2L, working_scales, "/")
  # Remove the weighted mean before projecting; this makes constant traits
  # exactly zero and avoids projecting a potentially very large intercept.
  Yw <- sweep(y_centered, 2L, working_scales, "/") * sqrt_weights
  Yr <- Yw - Q %*% crossprod(Q, Yw)
  blocks <- .fit_blocks(ncol(F), block_size)
  D_raw <- numeric(ncol(F))
  monomorphic <- logical(ncol(F))
  for (j in blocks) {
    f <- F[, j, drop = FALSE]
    D_raw[j] <- if (penalty == "identity") 1 else colSums(2 * f * (1 - f) * weights)
    monomorphic[j] <- apply(f, 2L, function(v) all(v == v[1L]))
  }
  D <- pmax(D_raw, penalty_floor)
  names(D) <- names(D_raw) <- names(monomorphic) <- colnames(F)
  list(F = F, Y = Y, Ys = Ys, Yr = Yr, Z = Z, zqr = zqr, Q = Q,
       weights = weights, sqrt_weights = sqrt_weights, D = D,
       D_raw = D_raw, monomorphic = monomorphic, blocks = blocks,
       trait_sd = trait_sd, working_scales = working_scales,
       constant_traits = constant_traits)
}

.fit_xblock <- function(prepared, j) {
  Xw <- 2 * prepared$F[, j, drop = FALSE] * prepared$sqrt_weights
  Xr <- Xw - prepared$Q %*% crossprod(prepared$Q, Xw)
  Xr[, prepared$monomorphic[j]] <- 0
  # QR roundoff must not invent effects for a variant fully explained by Z.
  numerical_zero <- colSums(Xr^2) <= (20 * .Machine$double.eps)^2 * colSums(Xw^2)
  Xr[, numerical_zero] <- 0
  Xr
}

.fit_linear_predictor <- function(prepared, B) {
  out <- matrix(0, nrow(prepared$F), ncol(B),
                dimnames = list(rownames(prepared$F), colnames(B)))
  for (j in prepared$blocks)
    out <- out + 2 * prepared$F[, j, drop = FALSE] %*% B[j, , drop = FALSE]
  out
}

.fit_residual <- function(prepared, B) {
  out <- -prepared$Yr
  for (j in prepared$blocks)
    out <- out + .fit_xblock(prepared, j) %*% B[j, , drop = FALSE]
  out
}

.fit_gradient <- function(prepared, B, residual, lambda) {
  gradient <- lambda * prepared$D * B
  for (j in prepared$blocks)
    gradient[j, ] <- gradient[j, , drop = FALSE] +
      crossprod(.fit_xblock(prepared, j), residual)
  gradient
}

.fit_assemble <- function(data, prepared, B, settings, diagnostics,
                          objective, method) {
  dimnames(B) <- list(colnames(prepared$F), colnames(prepared$Y))
  beta <- sweep(B, 2L, prepared$working_scales, "*")
  genetic <- .fit_linear_predictor(prepared, beta)
  coefficients <- qr.coef(prepared$zqr,
                         (prepared$Y - genetic) * prepared$sqrt_weights)
  fitted <- prepared$Z %*% coefficients + genetic
  dimnames(fitted) <- dimnames(prepared$Y)
  residuals <- prepared$Y - fitted
  diagnostics$covariate_rank <- prepared$zqr$rank
  diagnostics$residual_pool_df <- nrow(prepared$F) - prepared$zqr$rank
  diagnostics$no_between_pool_variation <- names(prepared$monomorphic)[prepared$monomorphic]
  diagnostics$monomorphic_variants <- names(prepared$monomorphic)[
    prepared$monomorphic & prepared$F[1L, ] %in% c(0, 1)]
  diagnostics$penalty_floored_variants <- names(prepared$D)[prepared$D_raw < settings$penalty_floor]
  diagnostics$constant_traits <- prepared$constant_traits
  diagnostics$objective_scale <- if (settings$standardize_traits)
    "weighted-standardized phenotypes" else "original phenotype units"
  diagnostics$effect_units <- "phenotype units per additional alternate allele"
  structure(list(beta = beta,
                 beta_standardized = sweep(beta, 2L, prepared$trait_sd, "/"),
                 beta_working = B, fitted = fitted, residuals = residuals,
                 intercept = setNames(as.numeric(coefficients[1L, ]), colnames(beta)),
                 covariate_effects = coefficients[-1L, , drop = FALSE],
                 objective = objective, diagnostics = diagnostics,
                 weights = prepared$weights, penalty = prepared$D,
                 trait_scales = prepared$trait_sd, settings = settings,
                 method = method, inference_status = "exploratory", data = data),
            class = c(paste0("poogwam_", method), "poogwam_fit"))
}

#' Fit an exact pooled ridge regression
#' @export
fit_pool_ridge <- function(data, lambda = 1,
                           penalty = c("within_pool", "identity"),
                           penalty_floor = 1e-8, block_size = 1000L,
                           standardize_traits = FALSE) {
  .fit_scalar(lambda, "lambda", strict = TRUE)
  penalty <- match.arg(penalty)
  prep <- .fit_prepare(data, penalty, penalty_floor, block_size, standardize_traits)
  K <- nrow(prep$F)
  kernel <- matrix(0, K, K)
  for (j in prep$blocks) {
    A <- sweep(.fit_xblock(prep, j), 2L, sqrt(prep$D[j]), "/")
    kernel <- kernel + tcrossprod(A)
  }
  system <- kernel + diag(lambda, K)
  if (any(!is.finite(system)))
    stop("ridge system overflows; increase penalty_floor or rescale inputs", call. = FALSE)
  factor <- tryCatch(chol(system), error = function(e)
    stop("ridge system is numerically singular; increase lambda", call. = FALSE))
  dual <- backsolve(factor, forwardsolve(t(factor), prep$Yr))
  B <- matrix(0, ncol(prep$F), ncol(prep$Y),
              dimnames = list(colnames(prep$F), colnames(prep$Y)))
  for (j in prep$blocks)
    B[j, ] <- crossprod(.fit_xblock(prep, j), dual) / prep$D[j]
  residual <- .fit_residual(prep, B)
  gradient <- 2 * .fit_gradient(prep, B, residual, lambda)
  objective <- sum(residual^2) + lambda * sum(prep$D * B^2)
  if (!is.finite(objective) || any(!is.finite(B)))
    stop("nonfinite ridge solution; rescale the input phenotypes", call. = FALSE)
  values <- eigen(kernel, symmetric = TRUE, only.values = TRUE)$values
  rank_tolerance <- max(0, values[1L]) * K * .Machine$double.eps
  condition <- kappa(system, exact = TRUE)
  diagnostics <- list(converged = TRUE, iterations = 1L,
                      max_abs_gradient = max(abs(gradient)),
                      relative_gradient = max(abs(gradient)) / max(1, max(abs(B))),
                      system_condition = condition,
                      effective_design_rank = sum(values > rank_tolerance),
                      numerical_warning = if (condition > 1e12)
                        "ill-conditioned dual system; inspect gradient or increase lambda" else NULL)
  settings <- list(lambda = lambda, penalty = penalty, penalty_floor = penalty_floor,
                   block_size = block_size, standardize_traits = standardize_traits)
  .fit_assemble(data, prep, B, settings, diagnostics, objective, "ridge")
}

#' @export
print.poogwam_fit <- function(x, ...) {
  cat("PooGwaM", x$method, "fit (exploratory)\n")
  cat(nrow(x$fitted), "pools;", nrow(x$beta), "variants;", ncol(x$beta), "traits\n")
  cat("Converged:", x$diagnostics$converged, "; objective:",
      format(x$objective, digits = 6), "\n")
  invisible(x)
}

#' @export
predict.poogwam_fit <- function(object, newdata = NULL, covariates = NULL, ...) {
  if (length(list(...))) stop("unused prediction arguments", call. = FALSE)
  if (is.null(newdata)) {
    if (!is.null(covariates))
      stop("covariates require newdata", call. = FALSE)
    return(object$fitted)
  }
  if (inherits(newdata, "pool_gwas_data")) {
    if (!is.null(covariates))
      stop("pass covariates in newdata or separately, not both", call. = FALSE)
    covariates <- newdata$covariates
    F <- newdata$frequencies
  } else F <- newdata
  .fit_matrix(F, "newdata frequencies")
  if (any(F < 0 | F > 1) || !setequal(colnames(F), rownames(object$beta)))
    stop("newdata must have matching variant IDs and frequencies in [0, 1]", call. = FALSE)
  F <- F[, rownames(object$beta), drop = FALSE]
  required <- rownames(object$covariate_effects)
  Z <- matrix(1, nrow(F), 1L)
  if (length(required)) {
    .fit_matrix(covariates, "prediction covariates")
    if (!setequal(colnames(covariates), required) ||
        !setequal(rownames(covariates), rownames(F)))
      stop("prediction covariates must have matching pool and covariate IDs", call. = FALSE)
    Z <- cbind(Z, covariates[rownames(F), required, drop = FALSE])
  } else if (!is.null(covariates)) {
    stop("this model was fitted without covariates", call. = FALSE)
  }
  out <- Z %*% rbind(object$intercept, object$covariate_effects)
  for (j in .fit_blocks(ncol(F), object$settings$block_size))
    out <- out + 2 * F[, j, drop = FALSE] %*% object$beta[j, , drop = FALSE]
  dimnames(out) <- list(rownames(F), colnames(object$beta))
  out
}
