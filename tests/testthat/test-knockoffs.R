knockoff_fixture <- function() {
  ids <- sprintf("i%02d", seq_len(24))
  variants <- c("v1", "v2", "pair_column_1", "v1__knockoff", "v5", "v6")
  g <- matrix((rep(seq_len(24), 6) * rep(c(1, 2, 5, 7, 11, 13), each = 24) +
                rep(seq_len(6), each = 24) + rep(seq_len(24) %/% 4, 6)) %% 3,
              24, 6, dimnames = list(ids, variants))
  frequencies <- setNames(c(.2, .35, .45, .1, .6, .8), variants)
  copies <- generate_independent_knockoffs(g, frequencies, seed = 51)
  y <- cbind(first = seq_len(24) + 3 * g[, 1],
             second = sin(seq_len(24) / 3) + g[, 4])
  rownames(y) <- ids
  builder <- function(y) design_pools(y[, "first", drop = FALSE], k = 6, min_size = 2)
  list(g = g, copies = copies, y = y, frequencies = frequencies, builder = builder)
}

test_that("ridge knockoff statistics obey exact single and multiple swap symmetry", {
  f <- knockoff_fixture()
  original <- select_pool_knockoffs(f$g, f$copies, f$y, f$builder,
                                    fit_args = list(lambda = .7, block_size = 3))
  for (swap in list(1L, c(2L, 4L, 6L), seq_len(6))) {
    g <- f$g
    copies <- f$copies
    g[, swap] <- f$copies[, swap]
    copies[, swap] <- f$g[, swap]
    changed <- select_pool_knockoffs(g, copies, f$y, f$builder,
                                     fit_args = list(lambda = .7, block_size = 3))
    expected <- original$W
    expected[swap] <- -expected[swap]
    expect_identical(changed$W, expected)
    expect_identical(changed$fit$beta, original$fit$beta)
    expect_identical(changed$fit$data, original$fit$data)
  }
  expect_s3_class(original, "poogwam_knockoff")
  expect_identical(original$inference_status, "experimental_knockoff")
  expect_identical(original$model_validation, "not_certified")
  expect_identical(original$generator_provenance$generator, "independent_binomial")
  expect_false(any(c("p", "p_value", "p_values", "p_raw") %in% names(original)))
})

test_that("multitrait sparse-group selection also has exact swap symmetry", {
  f <- knockoff_fixture()
  args <- list(lambda_group = .005, lambda_l1 = .002, lambda_ridge = .7,
               tolerance = 1e-10, block_size = 3)
  original <- select_pool_knockoffs(f$g, f$copies, f$y, f$builder,
                                   method = "multitrait", fit_args = args)
  g <- f$g
  copies <- f$copies
  g[, c(1, 3)] <- f$copies[, c(1, 3)]
  copies[, c(1, 3)] <- f$g[, c(1, 3)]
  changed <- select_pool_knockoffs(g, copies, f$y, f$builder,
                                  method = "multitrait", fit_args = args)
  expected <- original$W
  expected[c(1, 3)] <- -expected[c(1, 3)]
  expect_identical(changed$W, expected)
  expect_true(original$fit$diagnostics$converged)
  expected_norm <- sqrt(rowSums(original$beta_standardized^2)) -
    sqrt(rowSums(original$knockoff_beta_standardized^2))
  expect_equal(original$W, expected_norm, tolerance = 1e-14)
  expect_identical(original$inference_unit, "variant_association_with_the_trait_vector")
})

test_that("named alignment and reserved-looking labels do not alter pairs", {
  f <- knockoff_fixture()
  original <- select_pool_knockoffs(f$g, f$copies, f$y, f$builder)
  shuffled <- select_pool_knockoffs(f$g[24:1, 6:1], f$copies[c(2:24, 1), c(3:6, 1:2)],
                                    f$y[c(5:24, 1:4), ], f$builder)
  expect_identical(shuffled$W[names(original$W)], original$W)
  expect_identical(shuffled$fit$beta, original$fit$beta)
  expect_equal(nrow(original$pair_map), ncol(f$g))
  expect_equal(length(unique(c(original$pair_map$original_column,
                               original$pair_map$knockoff_column))), 2 * ncol(f$g))
  expect_setequal(original$pair_map$variant_id, colnames(f$g))
  invalid <- f$copies
  colnames(invalid)[1] <- "absent"
  expect_error(select_pool_knockoffs(f$g, invalid, f$y, f$builder), "matching")
  colnames(invalid)[1] <- colnames(invalid)[2]
  expect_error(select_pool_knockoffs(f$g, invalid, f$y, f$builder), "unique")
  invalid <- f$copies
  invalid[1, 1] <- .5
  expect_error(select_pool_knockoffs(f$g, invalid, f$y, f$builder), "0/1/2")
  invalid[1, 1] <- NA_real_
  expect_error(select_pool_knockoffs(f$g, invalid, f$y, f$builder), "missing")
})

test_that("knockoff+ threshold agrees with an independent exhaustive reference", {
  brute <- function(w, q) {
    candidates <- sort(unique(abs(w[w != 0])))
    accepted <- candidates[vapply(candidates, function(t)
      (1 + sum(w <= -t)) / max(1, sum(w >= t)) <= q, FALSE)]
    if (length(accepted)) min(accepted) else Inf
  }
  vectors <- list(rep(0, 8), c(2, 2, 2, 2, -1, 0), c(3, 2, 1, -3, -2, 0),
                  rep(1, 19), rep(1, 20), c(rep(3, 10), rep(1, 15), -2))
  for (w in vectors) for (q in c(.05, .1, .25, .5, 1))
    expect_identical(.knockoff_threshold(w, q), brute(w, q))
  expect_identical(.knockoff_threshold(rep(1, 19), .05), Inf)
  expect_identical(.knockoff_threshold(rep(1, 20), .05), 1)
  expect_error(.knockoff_threshold(c(1, NA_real_), .1), "finite")
  expect_error(.knockoff_threshold(1, 0), "positive")
})

test_that("identical copies, monomorphs and constant traits yield no false contrasts", {
  f <- knockoff_fixture()
  identical <- select_pool_knockoffs(f$g, f$g, f$y, f$builder,
                                     fit_args = list(block_size = 3))
  expect_identical(unname(identical$W), rep(0, ncol(f$g)))
  expect_identical(identical$threshold, Inf)
  expect_length(identical$selected, 0)
  expect_setequal(identical$diagnostics$identical_pairs, colnames(f$g))
  f$g[, 1] <- f$copies[, 1] <- 0
  with_constant <- select_pool_knockoffs(f$g, f$copies, f$y, f$builder)
  expect_identical(unname(with_constant$W[1]), 0)
  expect_length(with_constant$W, ncol(f$g))
  y <- f$y * 0
  fixed_builder <- function(y) setNames(rep(c("a", "b", "c"), each = 8), rownames(y))
  all_constant <- select_pool_knockoffs(f$g, f$copies, y, fixed_builder)
  expect_identical(unname(all_constant$W), rep(0, ncol(f$g)))
})

test_that("allele metadata is aligned and ambiguous orientations are rejected", {
  f <- knockoff_fixture()
  a <- data.frame(variant_id = colnames(f$g), ref = "A", alt = "G")
  correct <- select_pool_knockoffs(f$g, f$copies, f$y, f$builder,
    allele_metadata = list(genotypes = a, knockoffs = a[6:1, ]))
  expect_identical(correct$allele_metadata$variant_id, colnames(f$g))
  expect_match(correct$allele_validation, "labels_match")
  flipped <- a
  flipped$ref[1] <- "G"
  flipped$alt[1] <- "A"
  expect_error(select_pool_knockoffs(f$g, f$copies, f$y, f$builder,
    allele_metadata = list(genotypes = a, knockoffs = flipped)), "orientation")
  expect_error(select_pool_knockoffs(f$g, f$copies, f$y, f$builder,
    allele_metadata = a), "data frames")
})

test_that("selection restores caller RNG even for a random pool builder or errors", {
  f <- knockoff_fixture()
  random_builder <- function(y) setNames(sample(rep(letters[1:4], each = 6)), rownames(y))
  set.seed(503)
  before <- .Random.seed
  first <- select_pool_knockoffs(f$g, f$copies, f$y, random_builder, seed = 77)
  expect_identical(.Random.seed, before)
  second <- select_pool_knockoffs(f$g, f$copies, f$y, random_builder, seed = 77)
  expect_identical(first, second)
  expect_false(identical(first$rng$start, first$rng$after_pool_builder))
  expect_error(select_pool_knockoffs(f$g, f$copies, f$y,
    function(y) { runif(1); stop("builder failed") }), "builder failed")
  expect_identical(.Random.seed, before)
  expect_error(select_pool_knockoffs(f$g, f$copies, f$y, f$builder,
    fit_args = list(lambda = 0)), "lambda")
  expect_identical(.Random.seed, before)
  expect_error(select_pool_knockoffs(f$g, f$copies, f$y, f$builder,
    fit_args = list(covariates = 1)), "forbidden")
  expect_warning(expect_error(select_pool_knockoffs(f$g, f$copies, f$y, f$builder,
    method = "multitrait", fit_args = list(lambda_group = .001, lambda_ridge = .001,
      max_iter = 1, tolerance = 1e-14)), "converged"), "max_iter")
  expect_identical(.Random.seed, before)
})

test_that("reference generator uses supplied probabilities rather than observed frequencies", {
  f <- knockoff_fixture()
  set.seed(871)
  before <- .Random.seed
  copies <- generate_independent_knockoffs(f$g, f$frequencies, seed = 33)
  expect_identical(.Random.seed, before)
  entirely_different <- f$g * 0
  expected <- generate_independent_knockoffs(entirely_different, f$frequencies, seed = 33)
  expect_identical(as.vector(copies), as.vector(expected))
  ids <- sort(colnames(f$g))
  set.seed(33)
  independent <- matrix(rbinom(length(f$g), 2,
    rep(f$frequencies[ids], each = nrow(f$g))), nrow(f$g), ncol(f$g),
    dimnames = list(sort(rownames(f$g)), ids))
  expect_identical(as.vector(copies), as.vector(independent[rownames(f$g), colnames(f$g)]))
  reordered <- generate_independent_knockoffs(f$g[24:1, 6:1], f$frequencies[6:1], seed = 33)
  expect_identical(as.vector(reordered[rownames(f$g), colnames(f$g)]), as.vector(copies))
  expect_error(generate_independent_knockoffs(f$g, unname(f$frequencies)), "IDs")
  expect_error(generate_independent_knockoffs(f$g, f$frequencies * 4), "probabilities")
  impossible <- f$frequencies
  impossible[1] <- 0
  expect_error(generate_independent_knockoffs(f$g, impossible), "support")
  expect_error(generate_independent_knockoffs(f$g, f$frequencies, seed = -1), "seed")
  expect_identical(attr(copies, "knockoff_provenance")$frequency_source,
                   "externally_supplied_as_known_not_estimated_by_this_function")
  expect_identical(attr(copies, "knockoff_provenance")$output_fingerprint,
                   attr(expected, "knockoff_provenance")$output_fingerprint)
})

test_that("new random state is removed on return and generator model is never certified", {
  f <- knockoff_fixture()
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  previous <- if (had) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had) assign(".Random.seed", previous, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  if (had) rm(".Random.seed", envir = .GlobalEnv)
  generated <- generate_independent_knockoffs(f$g, f$frequencies, seed = 1)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  selected <- select_pool_knockoffs(f$g, generated, f$y, f$builder,
                                    generator_provenance = list(model = "user declared exact"))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_identical(selected$model_validation, "not_certified")
  expect_identical(selected$generator_provenance$model, "user declared exact")
  expect_error(select_pool_knockoffs(f$g, generated, f$y, f$builder,
    generator_provenance = list("unnamed")), "uniquely named")
})

test_that("nonrestorable Box-Muller RNG is rejected without mutation", {
  f <- knockoff_fixture()
  kind <- RNGkind()
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  previous <- if (had) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (had) assign(".Random.seed", previous, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind(normal.kind = "Box-Muller")
  set.seed(17)
  rnorm(1)
  before <- .Random.seed
  expect_error(generate_independent_knockoffs(f$g, f$frequencies), "Box-Muller")
  expect_identical(.Random.seed, before)
  expect_error(select_pool_knockoffs(f$g, f$copies, f$y, f$builder), "Box-Muller")
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind()[2], "Box-Muller")
})
