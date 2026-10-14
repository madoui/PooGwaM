manhattan_fixture <- function() {
  data.frame(variant_id = c("b", "a", "c", "d"), chromosome = c("chr10", "chr2", "X", "chr2"),
    position = c(200, 100, 300, 400), p_maxT = c(.04, .2, 1, .005),
    score = c(2, 1, 0, 4), W = c(-2, 0, 3, -1))
}

test_that("Manhattan coordinates are aligned by ID with explicit chromosome labels", {
  x <- manhattan_fixture()
  map <- x[4:1, c("variant_id", "chromosome", "position")]
  out <- manhattan_data(x[c("variant_id", "p_maxT")], map)
  expect_identical(out$chromosome_labels, c("chr2", "chr10", "X"))
  expect_identical(out$data$SNP, c("a", "d", "b", "c"))
  expect_equal(out$data$BP, c(100, 400, 200, 300))
  expect_equal(out$data$P, c(.2, .005, .04, 1))
  expect_true(out$logp)
  expect_identical(manhattan_data(x, chromosome_order = c("X", "chr10", "chr2"))$chromosome_labels,
                   c("X", "chr10", "chr2"))
  expect_error(manhattan_data(x, chromosome_order = "X"), "exactly once")
  expect_false(manhattan_data(x, value = "score")$logp)
  expect_equal(sort(manhattan_data(x, value = "W")$data$P), sort(x$W))
})

test_that("Manhattan plotting never silently substitutes or fabricates p-values", {
  x <- manhattan_fixture()
  x$p_maxT[1] <- NA_real_
  expect_error(manhattan_data(x), "Missing coordinates or values")
  expect_warning(out <- manhattan_data(x, na_action = "omit"), "1 variants omitted")
  expect_identical(out$omitted_variant_ids, "b")
  x$p_maxT[1] <- 0
  expect_error(manhattan_data(x), "zero values")
  x$p_maxT[1] <- Inf
  expect_error(manhattan_data(x), "finite")
  x <- manhattan_fixture()
  x$position[1] <- 1.5
  expect_error(manhattan_data(x), "base-pair integers")
  expect_error(manhattan_data(manhattan_fixture(), manhattan_fixture()[-1, ]), "Missing coordinates")
  expect_error(manhattan_data(rbind(manhattan_fixture(), manhattan_fixture()[1, ])), "unique")
})

test_that("one trait or omnibus is chosen without pooling testing families", {
  x <- manhattan_fixture()
  a <- transform(x, test = "trait", trait = "height", p_raw = p_maxT)
  b <- transform(x, test = "trait", trait = "mass", p_raw = p_maxT / 2)
  c <- transform(x, test = "omnibus", trait = NA_character_, p_raw = p_maxT / 3)
  combined <- rbind(a, b, c)
  expect_error(manhattan_data(combined, value = "p_raw", test = "trait"), "Choose one trait")
  out <- manhattan_data(combined, value = "p_raw", test = "trait", trait = "mass")
  expect_equal(sort(out$data$P), sort(b$p_raw))
  expect_equal(sort(manhattan_data(combined, value = "p_raw")$data$P), sort(c$p_raw))
  expect_error(manhattan_data(combined, trait = "height"), "only used")
})

test_that("qqman renders p-values, signed scores and a single named chromosome", {
  skip_if_not_installed("qqman")
  target <- tempfile(fileext = ".pdf")
  grDevices::pdf(target)
  on.exit({grDevices::dev.off(); unlink(target)}, add = TRUE)
  x <- manhattan_fixture()
  out <- plot_manhattan(x, threshold = .05, highlight = "d")
  expect_identical(out$backend, "qqman")
  expect_identical(out$threshold, .05)
  expect_silent(plot_manhattan(x, value = "W"))
  expect_identical(plot_manhattan(x, value = "W", threshold = 0)$threshold, 0)
  x$chromosome <- "X"
  expect_silent(plot_manhattan(x))
  x$p_maxT <- 1
  expect_silent(plot_manhattan(x))
  expect_error(plot_manhattan(x, logp = FALSE), "explicit plotting")
  expect_error(plot_manhattan(x, threshold = 0), "threshold")
  expect_error(plot_manhattan(x, highlight = "missing"), "plotted variant IDs")
})

test_that("CMplot renders through the optional backend without writing files", {
  skip_if_not_installed("CMplot")
  target <- tempfile(fileext = ".pdf")
  grDevices::pdf(target)
  on.exit({grDevices::dev.off(); unlink(target)}, add = TRUE)
  out <- plot_manhattan(manhattan_fixture(), backend = "CMplot", threshold = .05)
  expect_identical(out$backend, "CMplot")
  out <- plot_manhattan(manhattan_fixture(), backend = "CMplot", value = "W")
  expect_false(out$logp)
  expect_identical(plot_manhattan(manhattan_fixture(), backend = "CMplot",
                                value = "W", threshold = 0)$threshold, 0)
  single <- manhattan_fixture()[1, ]
  expect_error(plot_manhattan(single, backend = "CMplot"), "backend='qqman'")
  duplicate_position <- manhattan_fixture()[1:2, ]
  duplicate_position$chromosome <- "X"
  duplicate_position$position <- 100
  expect_error(plot_manhattan(duplicate_position, backend = "CMplot"), "single chromosome")
  expect_silent(plot_manhattan(single))
  expect_silent(plot_manhattan(duplicate_position))
})
