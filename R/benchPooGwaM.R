#' benchPooGwaM
#'
#' benchmarking performance of PooGwaM on simulated data
#'
#' @param pvalues data frame of phenotypes
#' @param causal.true numeric vector of the clusters
#' @param threshold float default=0.1
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
#' @importFrom ROCR performance prediction
#' @importFrom stats p.adjust
#' @export

benchPooGwaM<-function( pvalues, causal.true, threshold=0.1 ){
  causal.estimated=factor(p.adjust(pvalues,"BH")< threshold, levels=c(FALSE,TRUE))
  test.table<-table(causal.estimated,causal.true)
  TN <- test.table[1,1]
  FN <- test.table[1,2]
  FP <- test.table[2,1]
  TP <- test.table[2,2]

  # Calculate FDP
  if ((FP + TP) == 0) {
    FDP <- 0
  } else {
    FDP <- FP / (FP + TP)
  }
  FDP = round(FDP, digits = 3)
  pred <- prediction( 1-pvalues, causal.true)
  AUC = round( performance(pred,"auc")@y.values[[1]], digits = 3 )

  # calculate Precision
  if(TP+FP == 0){
    Precision = 0
  }
  else{
    Precision = round(TP/(TP+FP), digits = 3)
  }

  # calculate Recall
  if(TP+FN == 0){
    Recall = 0
  }
  else{
    Recall = round(TP/(TP+FN), digits = 3)
  }

  # Calculate F1-score
  if(Precision+Recall == 0){
    F1 = 0
  }
  else{
    F1 = round((2*Precision*Recall)/(Precision+Recall), digits = 3)
  }

  return (list( "TN" = TN,
                "FN" = FN,
                "FP" = FP,
                "TP" = TP,
                "AUC" = AUC,
                "FDP" = FDP,
                "Recall" = Recall,
                "Precision" = Precision,
                "F1_score" = F1) )
}
