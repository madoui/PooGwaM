test_that("pool, variant, phenotype and count alignment is explicit", {
  d <- small_pool_data()
  C <- cbind(batch = c(0,0,1,1)); rownames(C) <- rownames(d$frequencies)
  a <- as_pool_gwas_data(d$frequencies, d$phenotypes[4:1, , drop=FALSE],
                        d$pool_sizes[4:1], covariates=C[4:1,,drop=FALSE])
  expect_identical(a$phenotypes, d$phenotypes)
  expect_identical(a$pool_sizes, d$pool_sizes)
  expect_equal(a$covariates, C)
  expect_error(as_pool_gwas_data(unname(d$frequencies), d$phenotypes, d$pool_sizes), "IDs")
  bad <- d$phenotypes; rownames(bad)[1] <- "unknown"
  expect_error(as_pool_gwas_data(d$frequencies,bad,d$pool_sizes), "IDs")
  bad <- d$frequencies; colnames(bad)[1] <- colnames(bad)[2]
  expect_error(as_pool_gwas_data(bad,d$phenotypes,d$pool_sizes), "unique")
  expect_true(validate_pool_data(d)$valid)
})

test_that("count type and QC prevent fabricated precision", {
  d <- small_pool_data(); F <- d$frequencies
  depth <- matrix(100, nrow(F), ncol(F), dimnames=dimnames(F))
  counts <- list(alt=round(F*depth),depth=depth,type="reads")
  x <- as_pool_gwas_data(NULL,d$phenotypes,d$pool_sizes,counts=counts)
  expect_equal(x$frequencies,F)
  counts$type <- NULL
  expect_error(as_pool_gwas_data(F,d$phenotypes,d$pool_sizes,counts=counts), "explicit type")
  counts$type <- "chromosomes"
  expect_error(as_pool_gwas_data(F,d$phenotypes,d$pool_sizes,counts=counts), "twice")
  counts$type <- "reads"; counts$depth[1,1] <- 0; counts$alt[1,1] <- 0
  expect_error(as_pool_gwas_data(NULL,d$phenotypes,d$pool_sizes,counts=counts), "zero depth")
  x <- as_pool_gwas_data(NULL,d$phenotypes,d$pool_sizes,counts=counts,missing="drop_variants")
  expect_equal(ncol(x$frequencies),3)
  expect_identical(x$qc$excluded$variant_id,"snp1")
  expect_equal(x$counts$alt, counts$alt[,-1,drop=FALSE])
  expect_error(as_pool_gwas_data(F,d$phenotypes,d$pool_sizes,min_depth=1), "requires counts")
  expect_identical(d$qc$no_between_pool_variation,"snp4")
  expect_error(as_pool_gwas_data(F+.5,d$phenotypes,d$pool_sizes), "\\[0,1\\]")
})

test_that("individual pooling conserves alleles and aligns all rows", {
  G <- matrix(c(0,1,2,1,0,2,1,2),4,2,
              dimnames=list(paste0("id",1:4),c("a","b")))
  Y <- cbind(y=c(1,3,4,8),z=c(5,4,2,1)); rownames(Y)<-rownames(G)
  cl <- setNames(c("lo","lo","hi","hi"),rownames(G))
  d <- pool_from_individuals(G,Y[4:1,,drop=FALSE],cl[c(2,1,4,3)])
  expect_equal(colSums(d$counts$alt),colSums(G))
  expect_equal(d$phenotypes["lo",],c(y=2,z=4.5))
  expect_equal(unname(d$pool_sizes),c(2,2))
  expect_identical(d$counts$type,"chromosomes")
  expect_equal(compute_group_MAFs(G,cl),d$frequencies)
  expect_error(pool_from_individuals(G,Y,unname(cl)), "IDs")
  expect_error(pool_from_individuals(G/2,Y,cl), "0/1/2")
  s <- pool_from_individuals(G,Y,cl,retain_individuals=TRUE)
  expect_equal(s$individuals$genotypes,G)
  validated <- validate_pool_data(s)
  expect_equal(validated$data$individuals,s$individuals)
  expect_equal(validated$data$within_pool_ss,s$within_pool_ss)
  expect_identical(validated$data$qc,validated$qc)
})

test_that("missingness never silently changes sample denominators", {
  d <- small_pool_data(); F <- d$frequencies; F[1,1]<-NA
  expect_error(as_pool_gwas_data(F,d$phenotypes,d$pool_sizes), "Missing")
  q <- as_pool_gwas_data(F,d$phenotypes,d$pool_sizes,missing="drop_variants")
  expect_equal(q$pool_sizes,d$pool_sizes)
  expect_identical(validate_pool_data(q)$data$qc$excluded,q$qc$excluded)
  Y <- d$phenotypes; Y[1,1]<-NA
  expect_error(as_pool_gwas_data(d$frequencies,Y,d$pool_sizes), "missing")
  F[1,1] <- Inf
  expect_error(as_pool_gwas_data(F,d$phenotypes,d$pool_sizes,missing="drop_variants"), "nonfinite")
})
