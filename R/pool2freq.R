#' Getting Pooled allele frequencies
#'
#' The following function returns the PAF of a paricular SNP
#'  in all the pools. The output is a vector whose ith element is PAF of
#'  the snp in ith pool. Argument 'clusters_genotype' should be the
#'  output of the 'cls_genotype' function. The next argument 'snp' should be
#'  a character string denoting the SNP, for e.g., 'SNP_111', 'SNP_176' etc.
#'
#' @param clusters_genotype a list of clusters with genotype.
#' @param clusters a list of clusters
#' @param phenotype a PhenotypeSimulator output
#'
#' @return a list of clusters with allele frequencies
#'
#' @examples
#' sim = sim2snp(100, 100, 10)
#' pools = phen2pool(sim, 3, 3)
#' GenotypesCluster = trugen2clust (pools, sim)
#' alleleFreq = pool2freq (pools, GenotypesCluster, sim)
#' @export
pool2freq <- function(clusters,clusters_genotype, phenotype){
  pool_allele_freq1 <- function(clusters_genotype, snp){
    k = length(clusters_genotype)
    a=rep(0,k)
    for(i in 1:k){
      a[i]=sum(clusters_genotype[[i]][,snp])/(2*(nrow(clusters_genotype[[i]])))
    }
    return(a)
  }
  SNPs = phenotype[["setup"]][["id_snps"]] # denote the IDs of the SNPs.

  # the following object 'allele_freq1' is pooled allele frequency matrix
  # whose rows denote the pools, and columns denote the SNPs.
  allele_freq1 = matrix(NA, nrow=length(clusters), ncol=length(SNPs))
  for(i in 1:length(SNPs)){
    allele_freq1[,i] = pool_allele_freq1(clusters_genotype = clusters_genotype, snp=SNPs[i])
  }
  colnames(allele_freq1)=SNPs
  return(allele_freq1)

}
