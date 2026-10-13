fit_fixture <- function(covariates = TRUE, traits = 2L) {
  F <- matrix((sin(seq_len(42L) * 1.7) + 1) / 2, 7L, 6L,
              dimnames = list(paste0("p", 1:7), paste0("v", 1:6)))
  Y <- cbind(height = c(3, -1, 4, 0, 2, 6, -2),
             mass = c(10, 7, 2, 8, -1, 4, 5))[, seq_len(traits), drop = FALSE]
  rownames(Y) <- rownames(F)
  C <- if (covariates) matrix(c(-2, 1, 0, 3, -1, 2, 4), 7L, 1L,
                             dimnames = list(rownames(F), "batch")) else NULL
  as_pool_gwas_data(F, Y, setNames(c(2, 5, 3, 8, 11, 4, 7), rownames(F)),
                    covariates = C)
}
