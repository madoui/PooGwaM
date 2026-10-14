calibration_fixture <- function() {
  ids <- paste0("i", seq_len(12))
  g <- cbind(a = rep(c(0, 1, 2), 4), b = rep(c(0, 2), 6), c = c(rep(0, 6), rep(2, 6)))
  y <- cbind(first = seq_len(12), second = c(3, 1, 4, 2, 9, 6, 7, 5, 8, 12, 10, 11))
  rownames(g) <- rownames(y) <- ids
  builder <- function(y) setNames(rep(c("low", "middle", "high"), each = 4), rownames(y))
  list(g = g, y = y, builder = builder)
}

test_that("individual calibration retains traits, IDs, reproducibility and caller RNG", {
  f <- calibration_fixture()
  set.seed(120)
  before <- .Random.seed
  cal <- calibrate_pool_permutation(f$g, f$y, f$builder, B = 23, seed = 5)
  expect_identical(.Random.seed, before)
  again <- calibrate_pool_permutation(f$g, f$y, f$builder, B = 23, seed = 5)
  expect_identical(cal, again)
  expect_identical(dimnames(cal$p_raw), list(colnames(f$g), colnames(f$y)))
  expect_true(all(cal$p_raw >= 1 / 24 & cal$p_raw <= 1))
  expect_true(all(cal$p_maxT_omnibus >= cal$p_omnibus))
  shuffled <- calibrate_pool_permutation(f$g[12:1, ], f$y, f$builder, B = 23, seed = 5)
  expect_identical(cal$p_raw, shuffled$p_raw)
  expect_identical(cal$p_omnibus, shuffled$p_omnibus)
  expect_identical(cal$inference_scope, "global_null_only")
  expect_error(calibrate_pool_permutation(f$g, f$y, function(y) stop("builder failure"),
                                        B = 2, seed = 17), "builder failure")
  expect_identical(.Random.seed, before)
})

test_that("calibration does not leave a new global random state", {
  f <- calibration_fixture()
  had_state <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  saved_state <- if (had_state) get(".Random.seed", envir = .GlobalEnv) else NULL
  restore <- function() {
    if (had_state) assign(".Random.seed", saved_state, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }
  on.exit(restore(), add = TRUE)
  if (had_state) rm(".Random.seed", envir = .GlobalEnv)
  calibrate_pool_permutation(f$g, f$y, f$builder, B = 2, seed = 3)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  restore()
})

test_that("permutation preserves genomic rows and stays inside supplied blocks", {
  f <- calibration_fixture()
  f$g[,] <- rep(c(0, 1, 2), each = 4)
  blocks <- setNames(rep(1:3, each = 4), rownames(f$g))
  expect_warning(cal <- calibrate_pool_permutation(f$g, f$y, f$builder, B = 11,
                                                  strata = blocks[12:1]), "exchangeable")
  expect_true(all(cal$p_raw == 1))
  expect_true(all(cal$p_omnibus == 1))
  expect_error(calibrate_pool_permutation(f$g, f$y, f$builder,
    strata = setNames(seq_len(12), rownames(f$g))), "no nontrivial")
  expect_error(calibrate_pool_permutation(f$g, f$y, f$builder, strata = unname(blocks)), "named")
})

test_that("constant responses and monomorphic variants produce p = 1", {
  f <- calibration_fixture()
  f$y[,] <- 7
  cal <- calibrate_pool_permutation(f$g, f$y, f$builder, B = 9)
  expect_true(all(cal$p_raw == 1))
  expect_true(all(cal$p_omnibus == 1))
  f$g[,] <- 0
  f$y[, 1] <- seq_len(12)
  cal <- calibrate_pool_permutation(f$g, f$y, f$builder, B = 9)
  expect_true(all(cal$p_raw == 1))
  expect_true(all(cal$p_omnibus == 1))
})

test_that("calibration is invariant to small phenotype measurement units", {
  f <- calibration_fixture()
  cal <- calibrate_pool_permutation(f$g, f$y, f$builder, B = 31, seed = 3)
  f$y[, 2] <- f$y[, 2] * 1e-18
  scaled <- calibrate_pool_permutation(f$g, f$y, f$builder, B = 31, seed = 3)
  expect_identical(scaled$p_raw, cal$p_raw)
  expect_identical(scaled$p_omnibus, cal$p_omnibus)
  expect_identical(scaled$p_maxT_omnibus, cal$p_maxT_omnibus)
})

test_that("small exact randomization agrees with enumerated tail probability", {
  g <- matrix(c(0, 0, 2, 2), 4, 1, dimnames = list(letters[1:4], "snp"))
  y <- matrix(c(0, 0, 1, 1), 4, 1, dimnames = list(letters[1:4], "trait"))
  builder <- function(y) setNames(c("A", "A", "B", "B"), rownames(y))
  # Of six equiprobable two-individual allocations, two attain the observed
  # absolute association. This tests the actual randomization mechanism.
  cal <- calibrate_pool_permutation(g, y, builder, B = 399, seed = 12)
  expect_equal(as.numeric(cal$p_raw), 1 / 3, tolerance = 0.08)
  expect_equal(unname(cal$p_omnibus), as.numeric(cal$p_raw))
})

test_that("association tables correct raw values once over declared families", {
  f <- calibration_fixture()
  cal <- calibrate_pool_permutation(f$g, f$y, f$builder, B = 19)
  tab <- association_table(cal$fit, cal, "BH")
  traits <- tab$test == "trait"
  omni <- tab$test == "omnibus"
  expect_equal(tab$p_adjusted[traits], p.adjust(as.vector(cal$p_raw), "BH"))
  expect_equal(tab$p_adjusted[omni], unname(p.adjust(cal$p_omnibus, "BH")))
  expect_equal(tab$p_maxT[omni], unname(cal$p_maxT_omnibus))
  expect_true(all(is.na(tab$p_maxT[traits])))
  tab0 <- association_table(cal$fit)
  expect_true(all(is.na(tab0$p_raw)))
  expect_true(all(is.na(tab0$p_adjusted)))
  expect_true(all(tab0$inference_status == "exploratory"))
  altered <- cal$fit
  altered$settings$lambda <- 4
  expect_error(association_table(altered, cal), "exact data, settings")
  altered <- cal$fit
  altered$data$phenotypes[1, 1] <- altered$data$phenotypes[1, 1] + 1
  expect_error(association_table(altered, cal), "exact data, settings")
  expect_equal(association_table(cal$fit, cal, "none")$p_adjusted, tab$p_raw)
})

test_that("inputs reject misalignment and unsupported calibration settings", {
  f <- calibration_fixture()
  y <- f$y
  rownames(y)[1] <- "other"
  expect_error(calibrate_pool_permutation(f$g, y, f$builder), "identical individual")
  expect_error(calibrate_pool_permutation(f$g, f$y, f$builder, B = 0), "B must")
  expect_error(calibrate_pool_permutation(f$g, f$y, f$builder, seed = NA), "seed must")
  expect_error(calibrate_pool_permutation(f$g, f$y, f$builder,
                                        fit_args = list(data = 1)), "forbidden")
  expect_error(calibrate_pool_permutation(f$g, f$y, f$builder,
                                        fit_args = list(1)), "uniquely named")
})

test_that("typed count bootstrap samples both stages and restores RNG", {
  pool_ids <- c("a", "b", "c")
  alt <- matrix(c(1, 9, 16, 0, 4, 18), 3, 2,
                dimnames = list(pool_ids, c("s1", "s2")))
  depth <- matrix(20, 3, 2, dimnames = dimnames(alt))
  y <- matrix(c(-2, 0, 2), 3, 1, dimnames = list(pool_ids, "trait"))
  sizes <- setNames(rep(10, 3), pool_ids)
  data <- as_pool_gwas_data(alt / depth, y, sizes,
                           counts = list(alt = alt, depth = depth, type = "reads"))
  set.seed(154)
  before <- .Random.seed
  expect_warning(cal <- calibrate_pool_counts(data, B = 23, seed = 13), "Experimental")
  expect_identical(.Random.seed, before)
  expect_warning(again <- calibrate_pool_counts(data, B = 23, seed = 13), "Experimental")
  expect_identical(cal, again)
  expect_true(all(cal$p_raw > 0 & cal$p_raw <= 1))
  expect_equal(cal$null_frequencies, colMeans(alt / depth))
  expect_identical(cal$inference_status, "experimental_parametric_global_null")
  data$counts$type <- "chromosomes"
  expect_warning(chrom <- calibrate_pool_counts(data, B = 23, seed = 13), "Experimental")
  expect_false(identical(cal$null_summary, chrom$null_summary))
  data$covariates <- matrix(seq_len(3), 3, 1)
  expect_error(calibrate_pool_counts(data), "does not support covariates")
  data$covariates <- NULL
  data$counts$type <- NULL
  expect_error(calibrate_pool_counts(data), "explicitly")
})

test_that("count bootstrap aligns named pool sizes and validates modified objects", {
  ids <- c("small", "medium", "large")
  sizes <- setNames(c(2, 5, 10), ids)
  alt <- matrix(c(1, 6, 18, 2, 4, 9), 3, 2,
                dimnames = list(ids, c("v1", "v2")))
  depth <- matrix(rep(2 * sizes, 2), 3, 2, dimnames = dimnames(alt))
  y <- matrix(c(1, 4, 7), 3, 1, dimnames = list(ids, "trait"))
  data <- as_pool_gwas_data(alt / depth, y, sizes,
    counts = list(alt = alt, depth = depth, type = "chromosomes"))
  expect_warning(cal <- calibrate_pool_counts(data, B = 13, seed = 3), "Experimental")
  changed <- data
  changed$pool_sizes <- changed$pool_sizes[3:1]
  expect_warning(reordered <- calibrate_pool_counts(changed, B = 13, seed = 3), "Experimental")
  expect_identical(reordered$p_raw, cal$p_raw)
  expect_identical(reordered$p_omnibus, cal$p_omnibus)
  expect_identical(reordered$null_summary, cal$null_summary)
  expect_identical(reordered$fit$data, validate_pool_data(changed)$data)
  changed$pool_sizes <- unname(changed$pool_sizes)
  expect_error(calibrate_pool_counts(changed, B = 3), "pool_sizes.*IDs")
  changed <- data
  changed$phenotypes <- NULL
  expect_error(calibrate_pool_counts(changed, B = 3), "phenotypes")
})
