studentized_fixture <- function() {
  ids <- paste0("person", seq_len(16))
  g <- cbind(rare = c(2, rep(0, 15)),
             common = rep(c(0, 1, 2, 1), 4),
             shifted = rep(c(0, 0, 1, 2), each = 4),
             constant = rep(0, 16))
  y <- cbind(first = seq_len(16), second = rep(c(1, 4, 2, 3), 4))
  rownames(g) <- rownames(y) <- ids
  builder <- function(y) setNames(rep(letters[1:4], each = 4), rownames(y))
  list(g = g, y = y, builder = builder)
}

# Independent dense calculation: no streaming, studentization or RNG helpers
# from the implementation are used. This also pins the legacy raw RNG order.
dense_calibration_reference <- function(f, B, seed, B_pilot = 0L,
                                        scale_floor = 1e-8) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed)
  design <- f$builder(f$y)
  fit <- fit_pool_ridge(pool_from_individuals(f$g, f$y, design))
  observed <- sqrt(rowSums(fit$beta_standardized^2))
  draw <- function() {
    g <- f$g[sample.int(nrow(f$g)), , drop = FALSE]
    rownames(g) <- rownames(f$y)
    fit_pool_ridge(pool_from_individuals(g, f$y, design))
  }
  pilot <- if (B_pilot) replicate(B_pilot,
    sqrt(rowSums(draw()$beta_standardized^2))) else NULL
  fits <- replicate(B, draw(), simplify = FALSE)
  norms <- vapply(fits, function(x) sqrt(rowSums(x$beta_standardized^2)),
                  numeric(ncol(f$g)))
  traits <- simplify2array(lapply(fits, function(x) abs(x$beta)))
  compare <- function(x, value) (1 + sum(x >= value * (1 - 64 * .Machine$double.eps))) / (B + 1)
  p_omnibus <- vapply(seq_along(observed), function(j) compare(norms[j, ], observed[j]), 0.0)
  p_raw <- abs(fit$beta)
  for (j in seq_len(nrow(p_raw))) for (k in seq_len(ncol(p_raw)))
    p_raw[j, k] <- compare(traits[j, k, ], abs(fit$beta[j, k]))
  maximum <- apply(norms, 2, max)
  if (B_pilot) {
    center <- rowMeans(pilot)
    sd <- apply(pilot, 1, stats::sd)
    floor <- max(scale_floor * max(c(center, sd)), .Machine$double.xmin)
    scale <- pmax(sd, floor)
    t_obs <- (observed - center) / scale
    t_null <- sweep(sweep(norms, 1, center, "-"), 1, scale, "/")
    tolerance <- 64 * .Machine$double.eps * pmax(abs(t_obs), observed / scale, center / scale)
    max_t <- apply(t_null, 2, max)
    p_max <- vapply(seq_along(observed), function(j)
      (1 + sum(max_t >= t_obs[j] - tolerance[j])) / (B + 1), 0.0)
  } else {
    p_max <- vapply(observed, function(s) compare(maximum, s), 0.0)
  }
  list(fit = fit, pilot = pilot, norms = norms, p_raw = p_raw,
       p_omnibus = p_omnibus, p_max = p_max, maximum = maximum,
       center = if (B_pilot) center else NULL, scale = if (B_pilot) scale else NULL,
       max_t = if (B_pilot) max_t else NULL, t_obs = if (B_pilot) t_obs else NULL)
}

test_that("studentized streaming matches a dense independent pilot reference", {
  f <- studentized_fixture()
  ref <- dense_calibration_reference(f, B = 31, seed = 14, B_pilot = 23)
  cal <- calibrate_pool_permutation(f$g, f$y, f$builder, B = 31, seed = 14,
    maxT = "studentized", B_pilot = 23)
  expect_equal(cal$p_raw, ref$p_raw)
  expect_equal(unname(cal$p_omnibus), unname(ref$p_omnibus))
  expect_equal(unname(cal$p_maxT_omnibus), unname(ref$p_max))
  expect_equal(cal$studentization$mean_omnibus, ref$center)
  expect_equal(cal$studentization$scale, ref$scale)
  expect_equal(cal$studentization$max_omnibus, ref$max_t)
  expect_equal(cal$studentization$observed_omnibus, ref$t_obs)
  expect_equal(cal$null_summary$mean_omnibus, rowMeans(ref$norms))
  expect_equal(cal$null_summary$sd_omnibus, apply(ref$norms, 1, stats::sd))
  expect_equal(cal$null_summary$max_omnibus, ref$maximum)
  expect_identical(cal$maxT_diagnostics$raw_max_variant,
    rownames(ref$norms)[apply(ref$norms, 2, which.max)])
  expect_true(any(ref$t_obs < 0))
  expect_true(all(cal$p_maxT_omnibus >= cal$p_omnibus))
  expect_true(all(cal$p_maxT_omnibus >= 1 / 32 & cal$p_maxT_omnibus <= 1))
  expect_equal(cal$p_maxT_omnibus["constant"], c(constant = 1))
  expect_identical(cal$maxT_method, "studentized")
  expect_identical(cal$inference_scope, "global_null_only")
  tab <- association_table(cal$fit, cal)
  expect_true(all(tab$maxT_method[tab$test == "omnibus"] == "studentized"))
  expect_true(all(is.na(tab$maxT_method[tab$test == "trait"])))
})

test_that("raw remains the legacy algorithm without pilot RNG consumption", {
  f <- studentized_fixture()
  ref <- dense_calibration_reference(f, B = 19, seed = 7)
  cal <- calibrate_pool_permutation(f$g, f$y, f$builder, B = 19, seed = 7)
  explicit <- calibrate_pool_permutation(f$g, f$y, f$builder, B = 19, seed = 7,
    maxT = "raw", B_pilot = 2, scale_floor = 0.5)
  expect_identical(cal, explicit)
  expect_identical(cal$maxT_method, "raw")
  expect_null(cal$studentization)
  expect_equal(cal$p_raw, ref$p_raw)
  expect_equal(unname(cal$p_omnibus), unname(ref$p_omnibus))
  expect_equal(unname(cal$p_maxT_omnibus), unname(ref$p_max))
  expect_equal(cal$null_summary$mean_omnibus, rowMeans(ref$norms))
  expect_equal(cal$null_summary$sd_omnibus, apply(ref$norms, 1, stats::sd))
  expect_equal(cal$null_summary$max_omnibus, ref$maximum)
  expect_identical(cal$rng$main_start, cal$rng$after_pooling)
})

test_that("RNG kinds and phase states replay studentization and random pooling", {
  f <- studentized_fixture()
  f$builder <- function(y) setNames(sample(rep(letters[1:4], each = 4)), rownames(y))
  old_kind <- RNGkind()
  set.seed(62)
  old_seed <- .Random.seed
  on.exit({do.call(RNGkind, as.list(old_kind)); assign(".Random.seed", old_seed, .GlobalEnv)})
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(100)
  before <- .Random.seed
  cal <- calibrate_pool_permutation(f$g, f$y, f$builder, B = 17, seed = 71,
    maxT = "studentized", B_pilot = 11)
  expect_identical(.Random.seed, before)
  expect_identical(cal$rng$pilot_start, cal$rng$after_pooling)
  expect_identical(cal$rng$main_start, cal$rng$pilot_end)
  expect_false(identical(cal$rng$pilot_start, cal$rng$main_start))
  expect_false(identical(cal$rng$initial_state, cal$rng$after_pooling))
  expect_identical(unname(cal$rng$kind), RNGkind())
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(300)
  before_replay <- .Random.seed
  again <- calibrate_pool_permutation(f$g, f$y, f$builder, B = 17, seed = cal$seed,
    maxT = "studentized", B_pilot = cal$studentization$B, rng_kind = cal$rng$kind)
  expect_identical(again, cal)
  expect_identical(.Random.seed, before_replay)
  expect_identical(RNGkind()[1], "L'Ecuyer-CMRG")
  expect_match(cal$provenance$input_fingerprint, "^[0-9a-f]{32}$")
  expect_match(cal$provenance$settings_fingerprint, "^[0-9a-f]{32}$")
  expect_identical(cal$rng$r_version, R.version.string)
  expect_error(calibrate_pool_permutation(f$g, f$y, function(y) stop("failed pooling"),
    B = 1, rng_kind = cal$rng$kind), "failed pooling")
  expect_identical(.Random.seed, before_replay)
  expect_error(calibrate_pool_permutation(f$g, f$y, f$builder, B = 1,
    rng_kind = c("invalid", "Inversion", "Rejection")), "not a valid")
  expect_identical(.Random.seed, before_replay)
})

test_that("studentized scores preserve units and individual alignment", {
  f <- studentized_fixture()
  cal <- calibrate_pool_permutation(f$g, f$y, f$builder, B = 17, seed = 4,
    maxT = "studentized", B_pilot = 19)
  f$y[, 2] <- f$y[, 2] * 1e-150
  scaled <- calibrate_pool_permutation(f$g[16:1, ], f$y, f$builder, B = 17, seed = 4,
    maxT = "studentized", B_pilot = 19)
  expect_identical(scaled$p_omnibus, cal$p_omnibus)
  expect_identical(scaled$p_maxT_omnibus, cal$p_maxT_omnibus)
  expect_equal(scaled$studentization$scale, cal$studentization$scale)
})

test_that("degenerate pilot scales do not assert that unobserved events are impossible", {
  fake_fit <- function(scores) {
    beta <- matrix(scores, ncol = 1, dimnames = list(c("a", "b"), "trait"))
    structure(list(beta = beta, beta_standardized = beta), class = "poogwam_fit")
  }
  fit <- fake_fit(c(2, 1))
  expect_warning(pilot <- PooGwaM:::.calibration_pilot(fit,
    function(b) fake_fit(c(0, 0)), B_pilot = 3, scale_floor = 1e-8), "does not establish")
  expect_identical(pilot$status, "all_zero_pilot_unit_scale")
  expect_equal(unname(pilot$scale), c(1, 1))
  cal <- PooGwaM:::.calibration_run(fit,
    function(b) fake_fit(if (b == 1) c(3, 0) else c(0, 0)),
    B = 3, method = "individual_permutation", metadata = list(), studentization = pilot)
  expect_equal(unname(cal$p_maxT_omnibus), c(0.5, 0.5))
  tiny <- PooGwaM:::.calibration_pilot(fake_fit(c(1e-300, 0)),
    function(b) fake_fit(c(b * 1e-300, 0)), B_pilot = 3, scale_floor = 1e-300)
  expect_true(all(is.finite(tiny$scale) & tiny$scale > 0))
  expect_equal(tiny$effective_floor, .Machine$double.xmin)
  expect_equal(unname(tiny$sd_omnibus[1]) / 1e-300, 1)
  expect_true(all(is.finite(PooGwaM:::.calibration_studentize(c(a = 1e-300, b = 0), tiny))))
  near_overflow <- PooGwaM:::.calibration_pilot(fake_fit(c(11, 10)),
    function(b) fake_fit(c(10, 10)), B_pilot = 3, scale_floor = 1e-309)
  bounded <- PooGwaM:::.calibration_run(fake_fit(c(11, 10)),
    function(b) fake_fit(c(10, 10)), B = 3, method = "individual_permutation",
    metadata = list(), studentization = near_overflow)
  expect_true(all(is.finite(bounded$studentization$comparison_tolerance)))
  expect_equal(unname(bounded$p_maxT_omnibus), c(1 / 4, 1))
  f <- studentized_fixture()
  f$g[,] <- 0
  expect_warning(all_zero <- calibrate_pool_permutation(f$g, f$y, f$builder,
    B = 3, maxT = "studentized", B_pilot = 3), "fixed unit scale")
  expect_true(all(all_zero$p_maxT_omnibus == 1))
})

test_that("a fixed pilot transform controls the complete tiny permutation orbit", {
  ids <- letters[1:4]
  g <- cbind(a = c(0, 0, 2, 2), b = c(2, 0, 0, 0), c = c(0, 1, 2, 1))
  y <- matrix(c(0, 0, 1, 1), ncol = 1, dimnames = list(ids, "trait"))
  rownames(g) <- ids
  design <- setNames(c("low", "low", "high", "high"), ids)
  permute <- function(x) {
    if (length(x) == 1) return(list(x))
    unlist(lapply(seq_along(x), function(i)
      lapply(permute(x[-i]), function(rest) c(x[i], rest))), recursive = FALSE)
  }
  fits <- lapply(permute(seq_len(4)), function(index) {
    shuffled <- g[index, , drop = FALSE]
    rownames(shuffled) <- ids
    fit_pool_ridge(pool_from_individuals(shuffled, y, design))
  })
  # Conditional validity must hold for a fixed admissible pilot realization,
  # not merely when averaging over a favorable set of pilots.
  pilot <- PooGwaM:::.calibration_pilot(fits[[1]], function(b) fits[[b]],
    B_pilot = 5, scale_floor = 1e-8)
  min_p <- vapply(seq_along(fits), function(i) {
    others <- setdiff(seq_along(fits), i)
    cal <- PooGwaM:::.calibration_run(fits[[i]], function(b) fits[[others[b]]],
      B = length(others), method = "individual_permutation", metadata = list(),
      studentization = pilot)
    min(cal$p_maxT_omnibus)
  }, 0.0)
  for (alpha in seq_len(length(fits)) / length(fits))
    expect_lte(mean(min_p <= alpha), alpha + 1e-14)
})

test_that("nonserializable RNG caches are rejected without changing the caller", {
  f <- studentized_fixture()
  old_kind <- RNGkind()
  set.seed(78)
  old_seed <- .Random.seed
  on.exit({do.call(RNGkind, as.list(old_kind)); assign(".Random.seed", old_seed, .GlobalEnv)})
  RNGkind(normal.kind = "Box-Muller")
  set.seed(19)
  ignored <- stats::rnorm(1)
  before <- .Random.seed
  expect_error(calibrate_pool_permutation(f$g, f$y, f$builder, B = 1), "normal cache")
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind()[2], "Box-Muller")
  second <- stats::rnorm(1)
  set.seed(19)
  expect_identical(c(ignored, second), stats::rnorm(2))
  RNGkind(normal.kind = "Inversion")
  set.seed(23)
  before <- .Random.seed
  expect_error(calibrate_pool_permutation(f$g, f$y, f$builder, B = 1,
    rng_kind = c("Mersenne-Twister", "Box-Muller", "Rejection")), "normal cache")
  expect_identical(.Random.seed, before)
  expect_error(calibrate_pool_permutation(f$g, f$y, f$builder, B = 1,
    rng_kind = c("user-supplied", "Inversion", "Rejection")), "cannot be captured")
  expect_identical(.Random.seed, before)
})

test_that("studentization validates controls and restores RNG after pilot failure", {
  f <- studentized_fixture()
  expect_error(calibrate_pool_permutation(f$g, f$y, f$builder,
    maxT = "studentized", B_pilot = 1), "B_pilot")
  expect_error(calibrate_pool_permutation(f$g, f$y, f$builder,
    maxT = "studentized", scale_floor = 0), "scale_floor")
  expect_error(calibrate_pool_permutation(f$g, f$y, f$builder,
    maxT = "studentized", scale_floor = 2), "scale_floor")
  expect_error(calibrate_pool_permutation(f$g, f$y, f$builder,
    rng_kind = "Mersenne-Twister"), "rng_kind")
  set.seed(99)
  before <- .Random.seed
  testthat::local_mocked_bindings(.calibration_pilot = function(...) stop("pilot failed"),
                                .package = "PooGwaM")
  expect_error(calibrate_pool_permutation(f$g, f$y, f$builder, B = 3,
    maxT = "studentized", B_pilot = 3), "pilot failed")
  expect_identical(.Random.seed, before)
})

test_that("count bootstrap records replayable RNG provenance without studentization", {
  f <- studentized_fixture()
  data <- pool_from_individuals(f$g, f$y, f$builder(f$y))
  old_kind <- RNGkind()
  set.seed(63)
  old_seed <- .Random.seed
  on.exit({do.call(RNGkind, as.list(old_kind)); assign(".Random.seed", old_seed, .GlobalEnv)})
  expect_warning(cal <- calibrate_pool_counts(data, B = 7, seed = 12), "Experimental")
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(2026)
  before <- .Random.seed
  expect_warning(replay <- calibrate_pool_counts(data, B = 7, seed = cal$seed,
    rng_kind = cal$rng$kind), "Experimental")
  expect_identical(cal, replay)
  expect_identical(.Random.seed, before)
  expect_identical(cal$maxT_method, "raw")
  expect_null(cal$studentization)
  expect_identical(cal$rng$initial_state, cal$rng$main_start)
  expect_false(identical(cal$rng$main_start, cal$rng$main_end))
  expect_identical(cal$rng$source_mode, "namespace")
  expect_identical(cal$rng$package_version, as.character(getNamespaceVersion("PooGwaM")))
  expect_match(cal$provenance$input_fingerprint, "^[0-9a-f]{32}$")
})
