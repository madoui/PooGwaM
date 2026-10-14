test_that("null and all-positive truth have defined confusion tables", {
  p <- setNames(c(.001, .01, .5, .8), letters[1:4])
  zero <- setNames(rep(FALSE, 4), letters[1:4])
  x <- evaluate_associations(p, zero, "p_raw", adjustment = "BH")
  expect_equal(x$summary$FP, 2)
  expect_equal(x$summary$TN, 2)
  expect_equal(x$summary$FDP, 1)
  expect_true(is.na(x$summary$recall))
  expect_true(is.na(x$summary$AUROC))
  expect_true(is.na(x$summary$average_precision))
  y <- evaluate_associations(p, !zero, "p_raw")
  expect_equal(y$summary$recall, .5)
  expect_equal(y$summary$FDP, 0)
  expect_equal(y$summary$average_precision, 1)
  expect_true(is.na(y$summary$AUROC))
  expect_true(is.na(y$summary$false_positive_rate))
  z <- evaluate_associations(p, zero, "score", top_n = 0)
  expect_equal(z$summary$FDP, 0)
  expect_equal(z$summary$n_selected, 0)
})

test_that("scores, p-values and adjustment are explicit and cannot be mixed", {
  p <- setNames(c(.01, .025, .5, .9), letters[1:4])
  truth <- setNames(c(TRUE, TRUE, FALSE, FALSE), letters[1:4])
  expect_error(evaluate_associations(p, truth), "input_type")
  raw <- evaluate_associations(p, truth, "p_raw", alpha = .05)
  adjusted <- evaluate_associations(setNames(raw$results$p_adjusted, names(p)), truth,
                                    "p_adjusted", alpha = .05)
  expect_equal(raw$summary, adjusted$summary)
  expect_equal(raw$summary$n_selected, 2)
  expect_error(evaluate_associations(p, truth, "p_adjusted", adjustment = "BH"), "never readjusted")
  expect_error(evaluate_associations(p * 10, truth, "p_raw"), "\\[0,1\\]")
  scores <- evaluate_associations(-p, truth, "score")
  expect_true(all(is.na(scores$results$p_raw)))
  expect_true(all(is.na(scores$results$selected)))
  expect_error(evaluate_associations(p, truth, "score", adjustment = "BH"), "only allowed")
})

test_that("AUC and average precision respect tied scores and explicit IDs", {
  score <- matrix(c(3, 2, 2, 1, 1, 2, 3, 4), 4, 2,
                  dimnames = list(letters[1:4], c("height", "weight")))
  truth <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), 4, 2, dimnames = dimnames(score))
  x <- evaluate_associations(score, truth, "score", top_n = 2)
  y <- evaluate_associations(score[c(4, 3, 1, 2), 2:1], truth[4:1, 2:1], "score", top_n = 2)
  expect_equal(x$summary$AUROC, c(.875, 1))
  expect_equal(x$summary$average_precision, c(5/6, 1))
  expect_equal(x$summary$n_selected, c(3, 2))
  expect_equal(x$summary[2:1, -1], y$summary[, -1], ignore_attr = TRUE)
  flipped <- evaluate_associations(-score, truth, "score", higher_better = FALSE, top_n = 2)
  expect_equal(x$summary, flipped$summary)
  bad <- truth
  rownames(bad)[1] <- "missing"
  expect_error(evaluate_associations(score, bad, "score"), "matching")
  expect_error(evaluate_associations(unname(score), truth, "score"), "IDs")
})

test_that("global and per-trait p adjustment use the announced family", {
  p <- matrix(c(.01, .02, .3, .4), 2, 2, dimnames = list(c("a", "b"), c("x", "y")))
  truth <- p < .05
  x <- evaluate_associations(p, truth, "p_raw", family = "global")
  y <- evaluate_associations(p, truth, "p_raw", family = "trait")
  expect_equal(x$results$p_adjusted, as.numeric(p.adjust(p, "BH")))
  expect_equal(y$results$p_adjusted, as.numeric(apply(p, 2, p.adjust, method = "BH")))
  expect_false(identical(x$results$p_adjusted, y$results$p_adjusted))
})

test_that("one variant and completely tied scores keep meaningful conventions", {
  one <- evaluate_associations(c(s1 = .2), c(s1 = TRUE), "score", top_n = 1)
  expect_equal(one$summary$TP, 1)
  expect_equal(one$summary$average_precision, 1)
  expect_true(is.na(one$summary$AUROC))
  score <- setNames(rep(2, 4), letters[1:4])
  truth <- setNames(c(TRUE, FALSE, TRUE, FALSE), letters[1:4])
  tied <- evaluate_associations(score, truth, "score", top_n = 1)
  expect_equal(tied$summary$AUROC, .5)
  expect_equal(tied$summary$average_precision, .5)
  expect_equal(tied$summary$n_selected, 4)
  expect_equal(tied$summary$recall, 1)
})

test_that("AUROC pair counts do not overflow on large balanced classes", {
  truth <- setNames(rep(c(FALSE, TRUE), each = 50000L), paste0("v", seq_len(100000L)))
  for (direction in c(1, -1, 0)) {
    scores <- setNames(direction * as.numeric(truth), names(truth))
    expect_warning(result <- evaluate_associations(scores, truth, "score"), NA)
    expect_equal(result$summary$AUROC, (direction + 1) / 2)
    expect_true(is.finite(result$summary$average_precision))
  }
})
