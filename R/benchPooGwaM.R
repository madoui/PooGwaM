#' benchPooGwaM
#'
#' benchmarking performance of PooGwaM on simulated data
#' plots ROC curve, with AUC and FDP from threshold
#'
#' @param pvalues data frame of phenotypes
#' @param causal.true numeric vector of the clusters
#' @param threshold float default=0.4
#'
#' @return list benchmarking metrics
#'
#' @examples
#' sim = PhenoSim (500, 100, 10, 0.6, 2)
#' causal.true = sim$causal.true
#' K = 4
#' clusters = quantile_clustering (data.frame(sim$phenotypes), K)
#' Freq= compute_group_MAFs(sim$genotypes,as.factor(clusters))
#' test.stat<-GWalpha(sim$phenotypes, d=1, clusters, Freq)
#' benchPooGwaM (test.stat, causal.true)
#'
#'
#' @importFrom ROCR performance prediction
#' @importFrom  graphics par
#' @importFrom stats p.adjust
#' @export

benchPooGwaM<-function(pvalues,causal.true,threshold=0.4){

  fdp <- function(cont_table) {

    # Extract values from contingency table
    TN <- cont_table[1,1]
    FN <- cont_table[1,2]
    FP <- cont_table[2,1]
    TP <- cont_table[2,2]

    # Calculate FDP
    if ((FP + TP) == 0) {
      fdp <- 0
    } else {
      fdp <- FP / (FP + TP)
    }

    return(fdp)
  }

  pred <- prediction( 1-pvalues, causal.true)
  par(mfrow=c(1,2))

  # ROC curve
  perf <- performance(pred, "tpr", "fpr")
  #plot(perf,colorize=TRUE, lwd= 3,main= "ROC curve")
  # Precision Recall curve
  perf <- performance(pred, "prec", "rec")
  #plot(perf, colorize=TRUE, lwd= 3,main= "Precision/Recall")


  # Confusion table for given theshold
  causal.estimated=factor(p.adjust(pvalues,"BH")< threshold, levels=c(FALSE,TRUE))
  test.table<-table(causal.estimated,causal.true)
  #print(test.table<-table(causal.estimated,causal.true))

  #print(paste("FDP",fdp(test.table)),digits=3)
  # AUC
  #print(paste("AUC: ",round(performance(pred,"auc")@y.values[[1]],digits=3)))

  return (list("TN" = test.table[1,1],
          "FN" = test.table[1,2],
          "FP" = test.table[2,1],
          "TP" = test.table[2,2],
          "AUC" = round( performance(pred,"auc")@y.values[[1]], digits = 3 ),
          "FDP" = round( fdp(test.table), digits = 3 ) ) )
}
