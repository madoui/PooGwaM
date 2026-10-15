covariate_fit_fixture <- function() {
  ids <- paste0("pool", seq_len(8))
  F <- matrix(c(.05,.15,.2,.45,.4,.8,.6,.7,
                .12,.15,.75,.32,.1,.45,.65,.25,
                .7,.4,.12,.9,.2,.65,.3,.5,
                .2,.55,.7,.15,.8,.3,.45,.65), 8, 4,
              dimnames = list(ids, paste0("s", 1:4)))
  metadata <- data.frame(age = c(18,25,37,22,45,29,51,34),
    sex = factor(c("F","F","M","F","M","M","F","M"), levels = c("F","M")),
    row.names = ids)
  Y <- cbind(height = c(8,12,18,13,20,15,22,17),
             weight = c(3,7,4,11,9,5,6,12))
  rownames(Y) <- ids
  sizes <- setNames(c(4,6,5,8,3,7,4,9), ids)
  list(F = F, Y = Y, metadata = metadata, sizes = sizes)
}

test_that("formula covariates match an independent unpenalized primal reference", {
  f <- covariate_fit_fixture()
  data <- as_pool_gwas_data(f$F, f$Y, f$sizes,
                            covariates = ~age + sex, metadata = f$metadata)
  C <- cbind(age = f$metadata$age, sexM = as.numeric(f$metadata$sex == "M"))
  rownames(C) <- rownames(f$F)
  expect_equal(data$covariates, C, ignore_attr = TRUE)
  Z <- cbind("(Intercept)" = 1, C)
  X <- cbind(Z, 2 * f$F)
  weights <- f$sizes / sum(f$sizes)
  for (penalty in c("within_pool", "identity")) {
    lambda <- .7
    D <- if (penalty == "identity") rep(1, ncol(f$F)) else
      colSums(2 * f$F * (1 - f$F) * weights)
    regularization <- diag(c(rep(0, ncol(Z)), lambda * D))
    reference <- solve(crossprod(X * sqrt(weights)) + regularization,
                       crossprod(X * sqrt(weights), f$Y * sqrt(weights)))
    fitted <- fit_pool_ridge(data, lambda = lambda, penalty = penalty)
    expect_equal(unname(fitted$beta), unname(reference[-seq_len(ncol(Z)), , drop = FALSE]),
                 tolerance = 1e-10)
    expect_equal(unname(rbind(fitted$intercept, fitted$covariate_effects)),
                 unname(reference[seq_len(ncol(Z)), , drop = FALSE]), tolerance = 1e-10)
    expect_equal(unname(fitted$fitted), unname(X %*% reference), tolerance = 1e-10)
    expect_identical(fitted$covariate_spec, data$covariate_spec)
    expect_identical(fitted$diagnostics$covariate_design$columns, c("age", "sexM"))
    expect_identical(fitted$diagnostics$covariate_design$input_level, "pool")
    expect_equal(fitted$diagnostics$covariate_rank, 3)
    expect_equal(fitted$diagnostics$residual_pool_df, 5)
    expect_length(fitted$diagnostics$covariate_design$constant_columns, 0)
  }
})

test_that("ridge and sparse-group share the same covariate design", {
  f <- covariate_fit_fixture()
  data <- as_pool_gwas_data(f$F, f$Y, f$sizes,
                            covariates = ~age + sex, metadata = f$metadata)
  manual <- as_pool_gwas_data(f$F, f$Y, f$sizes, covariates = data$covariates)
  for (penalty in c("within_pool", "identity")) {
    fit <- fit_pool_multitrait(data, lambda_group = .03, lambda_l1 = .01,
      lambda_ridge = .7, penalty = penalty, tolerance = 1e-10)
    reference <- fit_pool_multitrait(manual, lambda_group = .03, lambda_l1 = .01,
      lambda_ridge = .7, penalty = penalty, tolerance = 1e-10)
    expect_true(fit$diagnostics$converged)
    expect_equal(fit$beta, reference$beta, tolerance = 1e-10)
    expect_equal(fit$covariate_effects, reference$covariate_effects, tolerance = 1e-10)
    expect_equal(predict(fit, data), fit$fitted, tolerance = 1e-10)
    expect_identical(fit$covariate_spec, data$covariate_spec)
    ridge <- fit_pool_ridge(data, lambda = .7, penalty = penalty, standardize_traits = TRUE)
    zero_sparse <- fit_pool_multitrait(data, lambda_group = 0, lambda_l1 = 0,
                                      lambda_ridge = .7, penalty = penalty)
    expect_equal(zero_sparse$beta, ridge$beta, tolerance = 1e-12)
    expect_identical(zero_sparse$covariate_spec, data$covariate_spec)
  }
})

test_that("prediction freezes categorical coding and transformations across new populations", {
  f <- covariate_fit_fixture()
  formula <- ~scale(age) + sex + I(age^2 / 1000)
  data <- as_pool_gwas_data(f$F, f$Y, f$sizes, covariates = formula, metadata = f$metadata)
  fit <- fit_pool_ridge(data, lambda = .7)
  F <- f$F[c(2, 5), , drop = FALSE]
  rownames(F) <- c("new1", "new2")
  metadata <- data.frame(age = c(60, 18),
    sex = factor(c("M", "F"), levels = c("M", "F")), row.names = rownames(F))
  manual <- cbind((metadata$age - mean(f$metadata$age)) / stats::sd(f$metadata$age),
                   as.numeric(metadata$sex == "M"), metadata$age^2 / 1000)
  dimnames(manual) <- list(rownames(F), rownames(fit$covariate_effects))
  expected <- cbind(1, manual) %*% rbind(fit$intercept, fit$covariate_effects) + 2 * F %*% fit$beta
  dimnames(expected) <- list(rownames(F), colnames(fit$beta))
  observed <- predict(fit, F, metadata = metadata)
  expect_equal(observed, expected, tolerance = 1e-10)
  expect_equal(predict(fit, F[1, , drop = FALSE], metadata = metadata[1, , drop = FALSE]),
               expected[1, , drop = FALSE], tolerance = 1e-10)
  expect_equal(predict(fit, F[2:1, 4:1], metadata = metadata), expected[2:1, ], tolerance = 1e-10)
  expect_equal(predict(fit, F, covariates = manual), expected, tolerance = 1e-10)
  new_Y <- f$Y[c(2, 5), , drop = FALSE]
  rownames(new_Y) <- rownames(F)
  new_data <- as_pool_gwas_data(F, new_Y, setNames(c(6, 3), rownames(F)),
    covariates = formula, metadata = metadata)
  # The same formula fitted to a new population would center/scale differently.
  # predict() must encode from raw metadata with the training terms instead.
  expect_false(isTRUE(all.equal(new_data$covariates, manual, check.attributes = FALSE)))
  expect_equal(predict(fit, new_data), expected, tolerance = 1e-10)
  changed_formula <- as_pool_gwas_data(F, new_data$phenotypes, new_data$pool_sizes,
    covariates = ~sex + age, metadata = metadata)
  expect_equal(predict(fit, changed_formula), expected, tolerance = 1e-10)
  unknown <- metadata
  unknown$sex <- c("M", "unknown")
  expect_error(predict(fit, F, metadata = unknown), "level|unknown")
  expect_error(predict(fit, F, metadata = metadata[1, , drop = FALSE]), "IDs|individual|pool")
  expect_error(predict(fit, F), "metadata")
  expect_error(predict(fit, F, covariates = manual, metadata = metadata), "not both")
  expect_error(predict(fit, data, covariates = manual), "not both")
})

test_that("individual formulas are encoded before pooling and replayed before aggregation", {
  individuals <- paste0("i", seq_len(32))
  groups <- setNames(rep(paste0("p", seq_len(8)), each = 4), individuals)
  ages <- c(18,22,24,28, 20,25,30,31, 30,31,32,35, 19,23,29,33,
            40,43,45,50, 23,28,35,38, 43,47,51,59, 29,30,39,44)
  sex <- c(rep("F",4), "F","F","F","M", "F","F","M","M", rep("M",4),
           "M","F","M","M", "F","M","F","F", "M","F","M","F", "M","M","M","F")
  metadata <- data.frame(age = ages, sex = factor(sex, levels = c("F","M")),
                          row.names = individuals)
  G <- matrix((seq_len(32 * 4) + rep(seq_len(32) %/% 3, 4)) %% 3, 32, 4,
              dimnames = list(individuals, paste0("s", 1:4)))
  Y <- cbind(trait = ages / 3 + G[, 1] + sin(seq_len(32)))
  rownames(Y) <- individuals
  data <- pool_from_individuals(G, Y, groups, covariates = ~scale(age) + sex + I(age^2),
                                metadata = metadata)
  C <- cbind((ages - mean(ages)) / stats::sd(ages), as.numeric(sex == "M"), ages^2)
  manual <- rowsum(C, groups, reorder = TRUE) / 4
  dimnames(manual) <- dimnames(data$covariates)
  expect_equal(data$covariates, manual, ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(as.numeric(manual[, 3]), as.numeric(rowsum(ages, groups) / 4)^2)))
  fit <- fit_pool_ridge(data, lambda = .8)
  expect_equal(predict(fit, data), fit$fitted, tolerance = 1e-10)
  expect_identical(fit$diagnostics$covariate_design$input_level, "individual")
  keep <- groups %in% c("p2", "p5", "p7")
  pools <- c("p7", "p2", "p5")
  F <- data$frequencies[pools, , drop = FALSE]
  expect_equal(predict(fit, F, metadata = metadata[keep, , drop = FALSE], clusters = groups[keep]),
               fit$fitted[pools, , drop = FALSE], tolerance = 1e-10)
  expect_error(predict(fit, F, metadata = metadata[keep, , drop = FALSE]), "clusters")
  bad <- data
  bad$pool_sizes[1] <- bad$pool_sizes[1] + 1
  expect_error(predict(fit, bad), "size|count")
  expect_error(predict(fit, F, metadata = metadata[keep, , drop = FALSE],
                        clusters = groups), "IDs|match")
})

test_that("post-pooling rank diagnostics name constant and aliased cofactor columns", {
  f <- covariate_fit_fixture()
  C <- cbind(age = f$metadata$age, age_copy = 2 * f$metadata$age, constant = 1)
  rownames(C) <- rownames(f$F)
  data <- as_pool_gwas_data(f$F, f$Y, f$sizes, covariates = C)
  err <- tryCatch(fit_pool_ridge(data), error = identity)
  expect_s3_class(err, "poogwam_covariate_rank_error")
  expect_match(conditionMessage(err), "full column rank")
  expect_match(conditionMessage(err), "age_copy")
  expect_match(conditionMessage(err), "constant")
  expect_equal(err$covariate_diagnostics$rank, 2)
  expect_setequal(err$covariate_diagnostics$aliased_columns, c("age_copy", "constant"))
  expect_identical(err$covariate_diagnostics$constant_columns, "constant")
  expect_error(fit_pool_multitrait(data), "full column rank")
  ids <- paste0("i", 1:16)
  G <- matrix(rep(c(0,1,2,1), 8), 16, 2,
              dimnames = list(ids, c("a","b")))
  Y <- matrix(1:16, 16, 1, dimnames = list(ids, "trait"))
  metadata <- data.frame(sex = factor(rep(c("F","M"), 8)), row.names = ids)
  pooled <- pool_from_individuals(G, Y, setNames(rep(letters[1:4], each = 4), ids),
                                 covariates = ~sex, metadata = metadata)
  expect_error(fit_pool_ridge(pooled), "constant columns: sexM")
  C <- cbind(a = c(0,1,0), b = c(0,0,1))
  rownames(C) <- rownames(f$F)[1:3]
  saturated <- as_pool_gwas_data(f$F[1:3, ], f$Y[1:3, ], f$sizes[1:3], covariates = C)
  expect_error(fit_pool_ridge(saturated), "residual pool degree of freedom.*3 pools")
})

test_that("legacy numeric covariates and intercept-only formulas remain explicit", {
  f <- covariate_fit_fixture()
  C <- matrix(f$metadata$age, ncol = 1, dimnames = list(rownames(f$F), "age"))
  legacy <- as_pool_gwas_data(f$F, f$Y, f$sizes, covariates = C)
  fit <- fit_pool_ridge(legacy)
  expect_null(fit$covariate_spec)
  expect_equal(predict(fit, legacy), fit$fitted, tolerance = 1e-10)
  expect_equal(predict(fit, f$F, C), fit$fitted, tolerance = 1e-10)
  expect_error(predict(fit, f$F, metadata = f$metadata), "covariate formula")
  formula_data <- as_pool_gwas_data(f$F, f$Y, f$sizes, covariates = ~age, metadata = f$metadata)
  expect_error(predict(fit, formula_data), "without a covariate formula")
  intercept <- as_pool_gwas_data(f$F, f$Y, f$sizes, covariates = ~1, metadata = f$metadata)
  only <- fit_pool_ridge(intercept)
  expect_length(only$covariate_effects, 0)
  expect_equal(predict(only, f$F), only$fitted, tolerance = 1e-10)
  expect_equal(predict(only, intercept), only$fitted, tolerance = 1e-10)
  expect_error(predict(fit, metadata = f$metadata), "require newdata")
  expect_error(predict(fit, f$F, covariates = C, unknown = TRUE), "unused prediction")
})
