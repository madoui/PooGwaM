test_that("dual ridge equals a direct augmented primal reference", {
  for (traits in 1:2) for (covariates in c(FALSE, TRUE)) {
    dat <- fit_fixture(covariates, traits)
    fit <- fit_pool_ridge(dat, lambda = 0.37, block_size = 2L)
    X <- cbind(1, dat$covariates, 2 * dat$frequencies)
    w <- as.numeric(dat$pool_sizes) / sum(dat$pool_sizes)
    penalty <- c(rep(0, ncol(X) - nrow(fit$beta)), .37 * fit$penalty)
    reference <- solve(crossprod(X, w * X) + diag(penalty),
                       crossprod(X, w * dat$phenotypes))
    expect_equal(unname(fit$beta), unname(tail(reference, nrow(fit$beta))), tolerance = 1e-10)
    expect_equal(unname(fit$fitted), unname(X %*% reference), tolerance = 1e-10)
    expect_lt(fit$diagnostics$max_abs_gradient, 1e-10)
    expect_identical(fit$inference_status, "exploratory")
  }
})

test_that("ridge is invariant to allele recoding with an adjusted intercept", {
  dat <- fit_fixture()
  original <- fit_pool_ridge(dat, lambda = .6)
  dat$frequencies[, c(2, 5)] <- 1 - dat$frequencies[, c(2, 5)]
  flipped <- fit_pool_ridge(dat, lambda = .6)
  expected <- original$beta
  expected[c(2, 5), ] <- -expected[c(2, 5), ]
  expect_equal(flipped$beta, expected, tolerance = 1e-10)
  expect_equal(flipped$fitted, original$fitted, tolerance = 1e-10)
  expect_equal(flipped$intercept,
               original$intercept + 2 * colSums(original$beta[c(2, 5), , drop = FALSE]),
               tolerance = 1e-10)
})

test_that("pool, variant and covariate ordering is controlled by IDs", {
  dat <- fit_fixture()
  fit <- fit_pool_ridge(dat)
  shuffled <- dat
  shuffled$phenotypes <- dat$phenotypes[7:1, , drop = FALSE]
  shuffled$pool_sizes <- dat$pool_sizes[c(2:7, 1)]
  shuffled$covariates <- dat$covariates[c(4:7, 1:3), , drop = FALSE]
  shuffled$frequencies <- dat$frequencies[7:1, 6:1, drop = FALSE]
  fit_shuffled <- fit_pool_ridge(shuffled)
  expect_equal(fit_shuffled$beta[rownames(fit$beta), ], fit$beta, tolerance = 1e-10)
  expect_equal(fit_shuffled$fitted[rownames(fit$fitted), ], fit$fitted, tolerance = 1e-10)
  expect_equal(predict(fit, shuffled)[rownames(fit$fitted), ], fit$fitted, tolerance = 1e-10)
  bad <- dat$frequencies
  colnames(bad)[1] <- "unknown"
  expect_error(predict(fit, bad, covariates = dat$covariates), "matching variant")
  bad <- dat$covariates
  rownames(bad)[1] <- "unknown"
  expect_error(predict(fit, dat$frequencies, covariates = bad), "matching pool")
  expect_error(predict(fit, dat$frequencies), "prediction covariates")
})

test_that("blockwise ridge gives the same objective and coefficients", {
  dat <- fit_fixture()
  small <- fit_pool_ridge(dat, penalty = "identity", block_size = 1L)
  large <- fit_pool_ridge(dat, penalty = "identity", block_size = 1000L)
  expect_equal(small$beta, large$beta, tolerance = 1e-12)
  expect_equal(small$objective, large$objective, tolerance = 1e-12)
  expect_equal(predict(small), small$fitted)
})

test_that("constant variants receive zero coefficients and constant traits are explicit", {
  dat <- fit_fixture(FALSE)
  dat$frequencies[, 1:3] <- matrix(rep(c(0, 1, .3), each = 7), 7)
  dat$phenotypes[, 2] <- 4
  fit <- fit_pool_ridge(dat, standardize_traits = TRUE)
  expect_identical(unname(fit$beta[1:3, ]), matrix(0, 3, 2))
  expect_equal(unname(fit$fitted[, 2]), rep(4, 7), tolerance = 1e-12)
  expect_equal(fit$diagnostics$monomorphic_variants, c("v1", "v2"))
  expect_equal(fit$diagnostics$no_between_pool_variation, c("v1", "v2", "v3"))
  expect_equal(fit$diagnostics$constant_traits, "mass")
  expect_true(all(is.finite(fit$beta_standardized)))
})

test_that("trait scaling is reversible and standardized effects are always supplied", {
  dat <- fit_fixture()
  fit <- fit_pool_ridge(dat)
  scaled <- fit_pool_ridge(dat, standardize_traits = TRUE)
  expect_equal(scaled$beta, fit$beta, tolerance = 1e-11)
  expect_equal(scaled$fitted, fit$fitted, tolerance = 1e-11)
  expect_equal(fit$beta_standardized,
               sweep(fit$beta, 2, fit$trait_scales, "/"), tolerance = 1e-12)
  dat$phenotypes[, 1] <- 10 * dat$phenotypes[, 1] + 7
  rescaled <- fit_pool_ridge(dat)
  expect_equal(rescaled$beta[, 1], fit$beta[, 1] * 10, tolerance = 1e-10)
  expect_equal(rescaled$beta_standardized, fit$beta_standardized, tolerance = 1e-10)
})

test_that("ridge rejects unidentified nuisance models and invalid settings", {
  dat <- fit_fixture()
  for (value in list(0, -1, Inf, NA_real_, c(1, 2), "1"))
    expect_error(fit_pool_ridge(dat, lambda = value), "lambda")
  expect_error(fit_pool_ridge(dat, block_size = .5), "block_size")
  expect_error(fit_pool_ridge(dat, penalty_floor = 0), "penalty_floor")
  dat$covariates <- cbind(dat$covariates, duplicate = dat$covariates[, 1])
  expect_error(fit_pool_ridge(dat), "full column rank")
  dat <- fit_fixture(FALSE)
  dat$covariates <- diag(7)[, 1:6]
  dimnames(dat$covariates) <- list(rownames(dat$frequencies), paste0("c", 1:6))
  expect_error(fit_pool_ridge(dat), "residual pool degree")
})

test_that("ridge uses no random numbers", {
  set.seed(101)
  before <- .Random.seed
  fit_pool_ridge(fit_fixture())
  expect_identical(.Random.seed, before)
})

test_that("constant decimal traits and nuisance-explained variants stay exactly zero", {
  dat <- fit_fixture()
  dat$phenotypes[, 2] <- .1
  dat$pool_sizes <- setNames(c(1, 7, 13, 19, 23, 29, 31), rownames(dat$frequencies))
  dat$frequencies[, 1] <- (dat$covariates[, 1] + 2) / 6
  fit <- fit_pool_ridge(dat, standardize_traits = TRUE)
  expect_identical(unname(fit$beta[, 2]), rep(0, nrow(fit$beta)))
  expect_identical(unname(fit$beta[1, ]), c(0, 0))
  expect_equal(fit$diagnostics$constant_traits, "mass")
})

test_that("both fits reject count-frequency contradictions in modified S3 objects", {
  dat <- small_pool_data()
  depth <- matrix(100, nrow(dat$frequencies), ncol(dat$frequencies),
                  dimnames = dimnames(dat$frequencies))
  dat <- as_pool_gwas_data(phenotypes = dat$phenotypes, pool_sizes = dat$pool_sizes,
    counts = list(alt = round(dat$frequencies * depth), depth = depth, type = "reads"))
  dat$frequencies[1, 1] <- dat$frequencies[1, 1] + .1
  expect_error(fit_pool_ridge(dat), "inconsistent with counts")
  expect_error(fit_pool_multitrait(dat), "inconsistent with counts")
  expect_error(poogwam(dat), "inconsistent with counts")
})

test_that("fits retain aligned validated data and earlier QC exclusions", {
  dat <- small_pool_data()
  frequencies <- dat$frequencies
  frequencies[1, 1] <- NA_real_
  dat <- as_pool_gwas_data(frequencies, dat$phenotypes, dat$pool_sizes,
                           missing = "drop_variants")
  dat$phenotypes <- dat$phenotypes[4:1, , drop = FALSE]
  dat$pool_sizes <- dat$pool_sizes[4:1]
  dat$variants <- dat$variants[3:1, , drop = FALSE]
  dat$within_pool_ss <- "ancillary metadata retained"
  expected <- validate_pool_data(dat)$data
  for (fitter in list(fit_pool_ridge, fit_pool_multitrait)) {
    fit <- fitter(dat)
    expect_identical(fit$data, expected)
    expect_identical(fit$data$qc$excluded$variant_id, "snp1")
    expect_equal(fit$data$qc$n_input_variants, 4)
  }
  dat$variants$variant_id[1] <- "unknown"
  expect_error(fit_pool_ridge(dat), "variants IDs")
})
