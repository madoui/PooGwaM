covariate_individual_fixture <- function() {
  ids <- paste0("id", seq_len(8))
  list(G = matrix(rep(c(0, 1, 2, 1), 4), 8, 2,
                  dimnames = list(ids, c("snp1", "snp2"))),
       Y = matrix(seq_len(8), 8, 1, dimnames = list(ids, "trait")),
       clusters = setNames(rep(c("a", "b", "c", "d"), each = 2), ids),
       metadata = data.frame(age = c(1, 3, 2, 6, 4, 8, 5, 9),
         sexe = factor(c("F", "M", "F", "F", "M", "M", "F", "M"), levels = c("F", "M")),
         lot = c("A", "B", "A", "B", "A", "B", "B", "B"), row.names = ids))
}

test_that("individual formulas encode before averaging, including nonlinear terms", {
  x <- covariate_individual_fixture()
  formula <- ~ sexe + age + I(age^2) + age:sexe
  d <- pool_from_individuals(x$G, x$Y, x$clusters, covariates = formula, metadata = x$metadata)
  m <- model.matrix(formula, x$metadata)[, -1, drop = FALSE]
  expected <- rowsum(m, x$clusters) / 2
  expect_equal(d$covariates, expected, ignore_attr = TRUE)
  expect_equal(unname(d$covariates[, "sexeM"]), c(.5, 0, 1, .5))
  expect_equal(d$covariates["a", "I(age^2)"], 5)
  expect_false(d$covariates["a", "I(age^2)"] == d$covariates["a", "age"]^2)
  expect_equal(d$covariates["a", "sexeM:age"], 1.5)
  expect_identical(d$covariate_spec$input_level, "individual")
  expect_identical(d$covariate_clusters, x$clusters)
  expect_identical(validate_pool_data(d)$data$covariate_spec, d$covariate_spec)
  expect_equal(.pg_predict_covariates(d$covariate_spec, x$metadata, rownames(d$frequencies),
                                    x$clusters, d$pool_sizes), d$covariates)
  numeric <- pool_from_individuals(x$G, x$Y, x$clusters, m, TRUE)
  expect_identical(numeric$covariates, d$covariates)
  expect_null(numeric$covariate_spec)
  expect_equal(numeric$individuals$covariates, m, ignore_attr = TRUE)
})

test_that("metadata and cluster ID order never control covariate alignment", {
  x <- covariate_individual_fixture()
  d <- pool_from_individuals(x$G, x$Y, x$clusters, covariates = ~ age + sexe,
                            metadata = x$metadata)
  shuffled <- pool_from_individuals(x$G, x$Y, x$clusters[8:1], covariates = ~ age + sexe,
                                   metadata = x$metadata[c(3:8, 1:2), ])
  expect_equal(d$covariates, shuffled$covariates)
  expected <- d$covariates[4:1, , drop = FALSE]
  expect_equal(.pg_predict_covariates(d$covariate_spec, x$metadata[8:1, ],
    rownames(expected), x$clusters[c(5:8, 1:4)], d$pool_sizes[4:1]), expected)
  bad <- x$metadata; rownames(bad)[1] <- "unknown"
  expect_error(pool_from_individuals(x$G, x$Y, x$clusters, covariates = ~age, metadata = bad), "IDs")
  expect_error(.pg_predict_covariates(d$covariate_spec, x$metadata, rownames(expected)), "IDs|clusters")
  expect_error(.pg_predict_covariates(d$covariate_spec, x$metadata, rownames(expected),
    x$clusters, d$pool_sizes + 1), "sizes")
  expect_error(.pg_predict_covariates(d$covariate_spec, x$metadata, "unknown", x$clusters), "pool IDs")
})

test_that("pool-level formulas preserve the legacy numeric interface", {
  d <- small_pool_data()
  meta <- data.frame(age = 1:4, sexe = factor(c("F", "M", "F", "M")),
                     row.names = rownames(d$frequencies))
  form <- as_pool_gwas_data(d$frequencies, d$phenotypes, d$pool_sizes,
    covariates = ~ sexe + age, metadata = meta[4:1, ])
  expect_identical(form$covariate_spec$input_level, "pool")
  num <- as_pool_gwas_data(d$frequencies, d$phenotypes, d$pool_sizes,
                          covariates = form$covariates)
  expect_identical(num$covariates, form$covariates)
  expect_null(num$covariate_spec)
  expect_identical(validate_pool_data(form)$data$covariates, form$covariates)
  expect_identical(validate_pool_data(d)$data, d)
  expect_error(as_pool_gwas_data(d$frequencies, d$phenotypes, d$pool_sizes,
    covariates = form$covariates, metadata = meta), "requires a covariate formula")
  expect_error(.pg_predict_covariates(form$covariate_spec, meta, rownames(meta),
    clusters = setNames(rownames(meta), rownames(meta))), "Pool-level")
})

test_that("intercept-only formulas retain their specification without an empty matrix", {
  d <- small_pool_data()
  out <- as_pool_gwas_data(d$frequencies, d$phenotypes, d$pool_sizes, covariates = ~1)
  expect_null(out$covariates)
  expect_identical(out$covariate_spec$columns, character())
  expect_null(.pg_predict_covariates(out$covariate_spec, NULL, rownames(d$frequencies)))
  expect_true(validate_pool_data(out)$valid)
  x <- covariate_individual_fixture()
  individual <- pool_from_individuals(x$G, x$Y, x$clusters, covariates = ~1)
  expect_null(individual$covariates)
  expect_true(validate_pool_data(individual)$valid)
  expect_null(.pg_predict_covariates(individual$covariate_spec, NULL,
    rownames(individual$frequencies), x$clusters, individual$pool_sizes))
})

test_that("formula missingness, unsupported calls and environment lookup fail clearly", {
  d <- small_pool_data()
  meta <- data.frame(age = 1:4, row.names = rownames(d$frequencies))
  make <- function(formula, metadata = meta)
    as_pool_gwas_data(d$frequencies, d$phenotypes, d$pool_sizes,
                     covariates = formula, metadata = metadata)
  bad <- meta; bad$age[1] <- NA
  expect_error(make(~age, bad), "missing")
  bad$age[1] <- Inf
  expect_error(make(~age, bad), "nonfinite")
  sex <- c(0, 1, 0, 1)
  expect_error(make(~sex), "must be columns of metadata")
  expect_error(make(age ~ 1), "one-sided")
  expect_error(make(~age - 1), "intercept")
  expect_error(make(~0 + age), "intercept")
  expect_error(make(~offset(age)), "Unsupported")
  expect_error(make(~mean(age)), "Unsupported")
  expect_error(make(~ custom(age)), "Unsupported")
  expect_error(make(~I(scale(age))), "standalone")
  expect_error(make(~log(poly(age, 2))), "standalone")
  expect_error(make(~poly(scale(age), 2)), "standalone")
  expect_error(make(~I(age + missing_external)), "must be columns of metadata")
  expect_error(make(~.), "explicitly")
  expect_error(make(~age, NULL), "requires metadata")
  expect_error(suppressWarnings(make(~log(age - 4))), "missing|nonfinite|NaN")
})

test_that("factor levels and numerical contrasts are frozen for prediction", {
  ids <- paste0("p", 1:4)
  meta <- data.frame(group = factor(c("A", "B", "A", "B"), levels = c("A", "B")),
                     flag = c(TRUE, FALSE, TRUE, FALSE), row.names = ids)
  old <- options(contrasts = c("contr.treatment", "contr.poly"))
  on.exit(options(old))
  trained <- .pg_covariates(~ group + flag, meta, ids, "pool")
  options(contrasts = c("contr.sum", "contr.poly"))
  pred <- .pg_predict_covariates(trained$spec, meta[4:1, ], ids)
  expect_identical(pred, trained$matrix)
  absent <- meta; absent$group <- factor(rep("A", 4)); absent$flag <- rep(FALSE, 4)
  pred <- .pg_predict_covariates(trained$spec, absent, ids)
  expect_equal(unname(pred), matrix(0, 4, 2))
  unknown <- meta; unknown$group <- c("A", "B", "C", "A")
  expect_error(.pg_predict_covariates(trained$spec, unknown, ids), "new level")
  one <- meta; one$group <- factor(rep("A", 4))
  expect_error(.pg_covariates(~group, one, ids, "pool"), "2 or more levels: group")
})

test_that("scale and poly reuse their training parameters for new datasets", {
  ids <- paste0("p", 1:6)
  meta <- data.frame(age = c(2, 3, 5, 9, 11, 15), row.names = ids)
  trained <- .pg_covariates(~scale(age) + poly(age, 2), meta, ids, "pool")
  fresh <- data.frame(age = c(20, 31, 44), row.names = paste0("n", 1:3))
  predicted <- .pg_predict_covariates(trained$spec, fresh, rownames(fresh))
  expect_equal(unname(predicted[, "scale(age)"]),
               (fresh$age - mean(meta$age)) / sd(meta$age))
  expected_poly <- predict(poly(meta$age, 2), newdata = fresh$age)
  expect_equal(unname(predicted[, c("poly(age, 2)1", "poly(age, 2)2")]),
               unname(as.matrix(expected_poly)), ignore_attr = TRUE)
  one <- .pg_predict_covariates(trained$spec, fresh[1, , drop = FALSE], rownames(fresh)[1])
  expect_equal(one, predicted[1, , drop = FALSE])
  expect_false(isTRUE(all.equal(predicted[, "scale(age)"], as.numeric(scale(fresh$age)))))
  expect_equal(.pg_predict_covariates(trained$spec, meta, ids), trained$matrix)
  log_scaled <- .pg_covariates(~scale(log(age)), meta, ids, "pool")
  log_prediction <- .pg_predict_covariates(log_scaled$spec, fresh, rownames(fresh))
  expect_equal(unname(log_prediction[, 1]),
    (log(fresh$age) - mean(log(meta$age))) / sd(log(meta$age)))
  # A caller's helper with a familiar name cannot replace the supported function.
  scale <- function(x) stop("caller scale must not be evaluated")
  guarded <- .pg_covariates(~scale(age), meta, ids, "pool")
  expect_equal(unname(guarded$matrix[, 1]), as.numeric(base::scale(meta$age)))
})

test_that("validation detects stale matrices, metadata and specifications", {
  x <- covariate_individual_fixture()
  d <- pool_from_individuals(x$G, x$Y, x$clusters, covariates = ~age + sexe,
                            metadata = x$metadata)
  bad <- d; bad$covariates[1, 1] <- bad$covariates[1, 1] + 1
  expect_error(validate_pool_data(bad), "inconsistent")
  bad <- d; bad$covariate_metadata$age[1] <- 100
  expect_error(validate_pool_data(bad), "inconsistent")
  bad <- d; bad$covariate_spec$columns <- rev(bad$covariate_spec$columns)
  expect_error(validate_pool_data(bad), "columns")
  bad <- d; bad$covariate_clusters[1] <- "b"
  expect_error(validate_pool_data(bad), "sizes")
  bad <- d; bad$covariate_spec$predvars <- quote(list(age))
  expect_error(validate_pool_data(bad), "Inconsistent")
  bad <- d; bad$covariate_spec$predvars <- NULL
  attr(bad$covariate_spec$terms, "predvars") <- NULL
  expect_error(validate_pool_data(bad), "Invalid covariate_spec")
  bad <- d; bad$covariate_spec <- NULL
  expect_error(validate_pool_data(bad), "requires covariate_spec")
})
