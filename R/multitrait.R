# Sparse-group elastic net: the nonsmooth penalty acts on variant rows and
# shares selection across traits. The ridge reference remains separable.

.fit_prox <- function(B, group, l1) {
  out <- sign(B) * pmax(abs(B) - l1, 0)
  norms <- sqrt(rowSums(out^2))
  multiplier <- numeric(length(norms))
  positive <- norms > 0
  multiplier[positive] <- pmax(1 - group / norms[positive], 0)
  out * multiplier
}

.fit_smooth <- function(prepared, B, lambda) {
  residual <- .fit_residual(prepared, B)
  list(residual = residual,
       value = 0.5 * (sum(residual^2) + lambda * sum(prepared$D * B^2)))
}

.fit_nonsmooth <- function(B, group, l1) {
  group * sum(sqrt(rowSums(B^2))) + l1 * sum(abs(B))
}

#' Fit a pooled sparse-group multitrait regression
#' @export
fit_pool_multitrait <- function(data, lambda_group = 0.1, lambda_l1 = 0,
                                lambda_ridge = 1,
                                penalty = c("within_pool", "identity"),
                                penalty_floor = 1e-8,
                                standardize_traits = TRUE,
                                max_iter = 10000L, tolerance = 1e-8,
                                block_size = 1000L) {
  .fit_scalar(lambda_group, "lambda_group")
  .fit_scalar(lambda_l1, "lambda_l1")
  .fit_scalar(lambda_ridge, "lambda_ridge", strict = TRUE)
  .fit_scalar(max_iter, "max_iter", lower = 1, integer = TRUE)
  .fit_scalar(tolerance, "tolerance", strict = TRUE)
  penalty <- match.arg(penalty)
  settings <- list(lambda_group = lambda_group, lambda_l1 = lambda_l1,
                   lambda_ridge = lambda_ridge, penalty = penalty,
                   penalty_floor = penalty_floor,
                   standardize_traits = standardize_traits,
                   max_iter = max_iter, tolerance = tolerance, block_size = block_size)
  if (lambda_group == 0 && lambda_l1 == 0) {
    fit <- fit_pool_ridge(data, lambda_ridge, penalty, penalty_floor,
                          block_size, standardize_traits)
    fit$objective <- fit$objective / 2
    fit$settings <- settings
    fit$method <- "multitrait"
    fit$diagnostics$kkt_residual <- fit$diagnostics$max_abs_gradient / 2
    fit$diagnostics$relative_kkt_residual <- fit$diagnostics$relative_gradient / 2
    fit$diagnostics$selected_variants <- rownames(fit$beta)[rowSums(abs(fit$beta)) > 0]
    class(fit) <- c("poogwam_multitrait", "poogwam_fit")
    return(fit)
  }
  prep <- .fit_prepare(data, penalty, penalty_floor, block_size, standardize_traits)
  B <- matrix(0, ncol(prep$F), ncol(prep$Y),
              dimnames = list(colnames(prep$F), colnames(prep$Y)))
  kernel <- matrix(0, nrow(prep$F), nrow(prep$F))
  for (j in prep$blocks) kernel <- kernel + tcrossprod(.fit_xblock(prep, j))
  spectrum <- eigen(kernel, symmetric = TRUE, only.values = TRUE)$values
  L <- max(spectrum[1L], 0) + lambda_ridge * max(prep$D)
  if (!is.finite(L) || L <= 0)
    stop("invalid Lipschitz bound; rescale inputs or regularization", call. = FALSE)
  V <- B
  momentum <- 1
  current <- .fit_smooth(prep, B, lambda_ridge)$value
  history <- numeric(max_iter + 1L)
  history[1L] <- current
  converged <- FALSE
  backtracking <- 0L
  restarts <- 0L
  kkt <- Inf
  for (iteration in seq_len(max_iter)) {
    repeat {
      smooth_v <- .fit_smooth(prep, V, lambda_ridge)
      gradient <- .fit_gradient(prep, V, smooth_v$residual, lambda_ridge)
      accepted <- FALSE
      for (attempt in seq_len(100L)) {
        candidate <- .fit_prox(V - gradient / L, lambda_group / L, lambda_l1 / L)
        candidate[prep$monomorphic, ] <- 0
        smooth_candidate <- .fit_smooth(prep, candidate, lambda_ridge)
        delta <- candidate - V
        majorant <- smooth_v$value + sum(gradient * delta) + 0.5 * L * sum(delta^2)
        slack <- 50 * .Machine$double.eps * max(1, abs(smooth_v$value))
        if (is.finite(smooth_candidate$value) && smooth_candidate$value <= majorant + slack) {
          accepted <- TRUE
          break
        }
        L <- L * 2
        backtracking <- backtracking + 1L
        if (!is.finite(L)) break
      }
      if (!accepted) stop("multitrait line search failed; rescale inputs", call. = FALSE)
      candidate_objective <- smooth_candidate$value +
        .fit_nonsmooth(candidate, lambda_group, lambda_l1)
      if (candidate_objective <= current + slack || momentum == 1) break
      # Monotone restart removes an extrapolation that increases the criterion.
      V <- B
      momentum <- 1
      restarts <- restarts + 1L
    }
    gradient_candidate <- .fit_gradient(prep, candidate, smooth_candidate$residual,
                                        lambda_ridge)
    # A unit-step proximal-gradient mapping is zero exactly at a KKT point.
    mapping <- candidate - .fit_prox(candidate - gradient_candidate,
                                     lambda_group, lambda_l1)
    kkt <- max(abs(mapping))
    relative_kkt <- kkt / max(1, max(abs(candidate)))
    history[iteration + 1L] <- candidate_objective
    previous <- B
    B <- candidate
    current <- candidate_objective
    if (relative_kkt <= tolerance) {
      converged <- TRUE
      break
    }
    next_momentum <- (1 + sqrt(1 + 4 * momentum^2)) / 2
    V <- B + (momentum - 1) / next_momentum * (B - previous)
    momentum <- next_momentum
  }
  if (!converged)
    warning("multitrait solver reached max_iter before its KKT tolerance", call. = FALSE)
  diagnostics <- list(converged = converged, iterations = iteration,
                      kkt_residual = kkt, relative_kkt_residual = relative_kkt,
                      lipschitz_bound = L, backtracking_steps = backtracking,
                      acceleration_restarts = restarts,
                      objective_history = history[seq_len(iteration + 1L)],
                      effective_design_rank = sum(spectrum >
                        max(spectrum[1L], 0) * nrow(kernel) * .Machine$double.eps),
                      selected_variants = rownames(B)[rowSums(abs(B)) > 0])
  .fit_assemble(data, prep, B, settings, diagnostics, current, "multitrait")
}
