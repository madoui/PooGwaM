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
