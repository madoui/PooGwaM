test_that("quantiles handle per-trait counts and more than 52 bins", {
  Y <- cbind(a=1:120,b=rep(1:3,each=40)); rownames(Y)<-paste0("id",1:120)
  p <- design_pools(Y,k=c(3,2),min_size=1)
  expect_true(all(!is.na(p$clusters)))
  q <- design_pools(Y[,1,drop=FALSE],k=60,min_size=2)
  expect_equal(length(q$pool_sizes),60)
  expect_true(all(q$pool_sizes==2))
  expect_false(any(grepl("NA",levels(q$clusters))))
  expect_equal(quantile_clustering(Y,c(3,2)),p$clusters)
})

test_that("ties are kept intact and missing phenotypes rejected", {
  Y <- matrix(rep(1:5,each=6),ncol=1,dimnames=list(paste0("id",1:30),"trait"))
  p <- design_pools(Y,k=4,min_size=2)
  expect_true(all(vapply(split(as.character(p$clusters),Y[,1]),function(z)length(unique(z))==1,TRUE)))
  b <- p$bounds[order(p$bounds$lower),]
  expect_true(all(head(b$upper,-1)<tail(b$lower,-1)))
  bad<-Y;bad[1,1]<-NA
  expect_error(design_pools(bad),"missing")
  bad[,1]<-1
  expect_error(design_pools(bad),"fewer than two")
})

test_that("balanced partition is deterministic, bounded and ID invariant", {
  Y <- cbind(a=1:100,b=sin(1:100));rownames(Y)<-paste0("id",1:100)
  a <- design_pools(Y,k=7,method="balanced",min_size=8)
  b <- design_pools(Y[100:1,,drop=FALSE],k=7,method="balanced",min_size=8)
  expect_equal(length(a$pool_sizes),7)
  expect_true(all(a$pool_sizes>=8))
  expect_identical(as.character(a$clusters),as.character(b$clusters[names(a$clusters)]))
  expect_error(design_pools(Y,k=20,method="balanced",min_size=8),"infeasible")
  expect_error(design_pools(Y,k=c(2,3,4)),"per trait")
  twelve <- Y[1:12,1,drop=FALSE]
  feasible <- design_pools(twelve,k=3,method="balanced",min_size=4)
  expect_equal(as.numeric(feasible$pool_sizes),rep(4,3))
})

test_that("balanced pooling recovers a feasible tied design after a greedy dead end", {
  Y <- matrix(rep(1:4, c(4, 2, 2, 4)), ncol = 1,
              dimnames = list(paste0("id", 1:12), "trait"))
  design <- design_pools(Y, k = 3, method = "balanced", min_size = 3)
  expect_equal(as.numeric(design$pool_sizes), rep(4, 3))
  expect_identical(design$search$method, "backtracking")
  expect_gt(design$search$states_visited, 1)
  expect_true(all(vapply(split(as.character(design$clusters), Y[, 1]),
                        function(group) length(unique(group)) == 1L, TRUE)))
  reordered <- design_pools(Y[12:1, , drop = FALSE], k = 3,
                            method = "balanced", min_size = 3)
  expect_identical(as.character(design$clusters),
                   as.character(reordered$clusters[names(design$clusters)]))
  expect_error(design_pools(Y, k = 3, method = "balanced", min_size = 3,
                            search_budget = 1), "feasibility is unresolved")
  expect_error(design_pools(Y, k = 3, method = "balanced", search_budget = 0),
               "search_budget")
  impossible <- Y
  impossible[, 1] <- rep(1:2, each = 6)
  expect_error(design_pools(impossible, k = 3, method = "balanced", min_size = 3),
               "No feasible recursive")
})

test_that("small tied one-trait designs agree with independently enumerated partitions", {
  # In one dimension, every admissible partition has cuts between contiguous
  # tie groups. Enumerate those partitions directly, independently of the
  # production recursive splitting procedure.
  configurations <- expand.grid(rep(list(1:3), 4))
  expected <- actual <- logical()
  for (row in seq_len(nrow(configurations))) for (k in 2:4) {
    sizes <- as.numeric(configurations[row, ])
    boundaries <- utils::combn(1:3, k - 1L)
    feasible <- any(apply(boundaries, 2L, function(cut) {
      groups <- cut(seq_along(sizes), breaks = c(0, cut, 4), labels = FALSE)
      all(as.numeric(rowsum(sizes, groups)) >= 2)
    }))
    Y <- matrix(rep(1:4, sizes), ncol = 1,
                dimnames = list(paste0("id", seq_len(sum(sizes))), "trait"))
    result <- tryCatch(design_pools(Y, k, "balanced", min_size = 2),
                       error = function(e) e)
    expected <- c(expected, feasible)
    actual <- c(actual, inherits(result, "pool_design"))
    if (inherits(result, "pool_design"))
      expect_true(length(result$pool_sizes) == k && all(result$pool_sizes >= 2))
  }
  expect_identical(actual, expected)
})
