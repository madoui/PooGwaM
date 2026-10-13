test_that("zero group and sparse penalties recover ridge exactly", {
  dat <- fit_fixture()
  fit <- fit_pool_multitrait(dat, lambda_group = 0, lambda_l1 = 0,
                             lambda_ridge = .4)
  reference <- fit_pool_ridge(dat, lambda = .4, standardize_traits = TRUE)
  expect_equal(fit$beta, reference$beta, tolerance = 1e-12)
  expect_equal(fit$fitted, reference$fitted, tolerance = 1e-12)
  expect_equal(fit$objective, reference$objective / 2, tolerance = 1e-12)
  expect_true(fit$diagnostics$converged)
})

test_that("sparse group optimization meets KKT and decreases its objective", {
  dat <- fit_fixture()
  fit <- fit_pool_multitrait(dat, lambda_group = .03, lambda_l1 = .015,
                             lambda_ridge = .4, tolerance = 1e-10, block_size = 2)
  expect_true(fit$diagnostics$converged)
  expect_lte(fit$diagnostics$relative_kkt_residual, 1e-10)
  expect_true(all(diff(fit$diagnostics$objective_history) <= 1e-12))
  expect_equal(predict(fit, dat), fit$fitted, tolerance = 1e-12)
  # Independent stationarity check: active rows use the analytical group gradient.
  w <- as.numeric(dat$pool_sizes) / sum(dat$pool_sizes)
  error <- (fit$fitted - dat$phenotypes) / rep(fit$trait_scales, each = nrow(fit$fitted))
  gradient <- crossprod(2 * dat$frequencies, error * w) + .4 * fit$penalty * fit$beta_working
  B <- fit$beta_working
  for (j in seq_len(nrow(B))) {
    nonzero <- abs(B[j, ]) > 1e-8
    if (any(nonzero)) {
      group_gradient <- .03 * B[j, ] / sqrt(sum(B[j, ]^2))
      expect_lt(max(abs((gradient[j, ] + group_gradient + .015 * sign(B[j, ]))[nonzero])), 1e-7)
    } else {
      soft <- sign(gradient[j, ]) * pmax(abs(gradient[j, ]) - .015, 0)
      expect_lte(sqrt(sum(soft^2)), .03 + 1e-8)
    }
  }
})

test_that("very strong penalties select the nuisance-only model", {
  dat <- fit_fixture()
  fit <- fit_pool_multitrait(dat, lambda_group = 1e6)
  expect_identical(unname(fit$beta), matrix(0, 6, 2))
  expect_equal(length(fit$diagnostics$selected_variants), 0L)
  Z <- cbind(1, dat$covariates)
  w <- as.numeric(dat$pool_sizes) / sum(dat$pool_sizes)
  expected <- Z %*% solve(crossprod(Z, w * Z), crossprod(Z, w * dat$phenotypes))
  expect_equal(unname(fit$fitted), unname(expected), tolerance = 1e-10)
})

test_that("multitrait is invariant to allele orientation and positive trait units", {
  dat <- fit_fixture()
  original <- fit_pool_multitrait(dat, lambda_group = .03, lambda_l1 = .01,
                                  tolerance = 1e-10)
  flipped_data <- dat
  flipped_data$frequencies[, 1] <- 1 - dat$frequencies[, 1]
  flipped <- fit_pool_multitrait(flipped_data, lambda_group = .03, lambda_l1 = .01,
                                 tolerance = 1e-10)
  beta_expected <- original$beta
  beta_expected[1, ] <- -beta_expected[1, ]
  expect_equal(flipped$beta, beta_expected, tolerance = 1e-8)
  expect_equal(flipped$fitted, original$fitted, tolerance = 1e-8)
  dat$phenotypes[, 2] <- dat$phenotypes[, 2] * 17 + 100
  scaled <- fit_pool_multitrait(dat, lambda_group = .03, lambda_l1 = .01,
                                tolerance = 1e-10)
  expect_equal(scaled$beta_standardized, original$beta_standardized, tolerance = 1e-8)
  expect_equal(scaled$beta[, 2], original$beta[, 2] * 17, tolerance = 1e-8)
})

test_that("single-trait sparse group equals one-dimensional elastic net soft thresholding", {
  F <- matrix(c(.1, .3, .6, .9), 4, 1,
              dimnames = list(paste0("p", 1:4), "v1"))
  Y <- matrix(c(-2, -1, 1, 3), 4, 1, dimnames = list(rownames(F), "y"))
  dat <- as_pool_gwas_data(F, Y, setNames(c(1, 3, 2, 4), rownames(F)))
  fit <- fit_pool_multitrait(dat, lambda_group = .04, lambda_l1 = .02,
                             lambda_ridge = .3, standardize_traits = FALSE,
                             penalty = "identity", tolerance = 1e-12)
  w <- dat$pool_sizes / sum(dat$pool_sizes)
  x <- 2 * F[, 1] - sum(w * 2 * F[, 1])
  y <- Y[, 1] - sum(w * Y[, 1])
  score <- sum(w * x * y)
  expected <- sign(score) * max(abs(score) - .06, 0) / (sum(w * x^2) + .3)
  expect_equal(unname(fit$beta[1, 1]), unname(expected), tolerance = 1e-10)
})

test_that("iteration failure is explicit, settings are checked and RNG is unchanged", {
  dat <- fit_fixture()
  expect_warning(fit <- fit_pool_multitrait(dat, lambda_group = .001,
      lambda_ridge = .001, max_iter = 1, tolerance = 1e-14), "max_iter")
  expect_false(fit$diagnostics$converged)
  expect_error(fit_pool_multitrait(dat, lambda_group = -1), "lambda_group")
  expect_error(fit_pool_multitrait(dat, lambda_ridge = 0), "lambda_ridge")
  expect_error(fit_pool_multitrait(dat, max_iter = 1.5), "max_iter")
  set.seed(203)
  before <- .Random.seed
  fit_pool_multitrait(dat)
  expect_identical(.Random.seed, before)
})
