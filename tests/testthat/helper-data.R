small_pool_data <- function() {
  F <- matrix(c(.1,.3,.6,.9, .3,.5,.7,.9, .2,.25,.4,.65,
                .5,.5,.5,.5), 4, 4,
              dimnames = list(paste0("pool", 1:4), paste0("snp", 1:4)))
  Y <- cbind(height = c(-2,-1,1,2), mass = c(3,1,0,-1))
  rownames(Y) <- rownames(F)
  as_pool_gwas_data(F, Y, setNames(c(7,11,13,19), rownames(F)))
}
