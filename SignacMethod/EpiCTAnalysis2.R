#
#   Code to generate confusion Matrices, auROC, auPR and MCC details
#
#
#   Read in Ground Truth & Results
#
#install.packages("MLmetrics")
library(MLmetrics)
#install.packages("pROC")
library(pROC)
GroundTruth<-read.csv("RNAObjmeta2.csv")
ATACpredictions<-read.csv("ATACObjPred2.csv")
#
# Code added to allow filtering out items with outlier QC data
#
#GroundTruth<-GroundTruth[GroundTruth$percent.mt<5,]
GroundTruth<-GroundTruth[GroundTruth$nFeature_RNA>200,]
filteredAT<-ATACpredictions[ATACpredictions$X %in% GroundTruth$X,]
ATACpredictions<-filteredAT

ATACpredictions<-ATACpredictions[ATACpredictions$nCount_peaks>3000,]
#ATACpredictions<-ATACpredictions[ATACpredictions$nCount_peaks<50000,]
#
# code added in case of different row counts
# code to try filtering by differnt QA
#
#
filteredGT<-GroundTruth[GroundTruth$X %in% ATACpredictions$X,]
GroundTruth<-filteredGT
#
#Function To be used by each Cell Type
#

plotROC<-function(CellType,SeurPred,SeurTruth)
{
  ACTvbl<-gsub(" ", ".", CellType)
  y_truth<-SeurTruth$predicted.id==CellType
  y_pred<-SeurPred$predicted.id==CellType
  print(CellType)
  print(ConfusionMatrix(y_pred,y_truth))
  evalstr=paste("y_prob<-SeurPred$prediction.score.",ACTvbl,sep="")
  eval(parse(text=evalstr))
  roc_obj <- roc(y_truth, y_prob)

  
  
  
}
mcc <- function (actual, predicted)
{
  # handles zero denominator and verflow error on large-ish products in denominator.
  #
  # actual = vector of true outcomes, 1 = Positive, 0 = Negative
  # predicted = vector of predicted outcomes, 1 = Positive, 0 = Negative
  # function returns MCC
  
  TP <- sum(actual == 1 & predicted == 1)
  TN <- sum(actual == 0 & predicted == 0)
  FP <- sum(actual == 0 & predicted == 1)
  FN <- sum(actual == 1 & predicted == 0)
  
  sum1 <- TP+FP; sum2 <-TP+FN ; sum3 <-TN+FP ; sum4 <- TN+FN;
  denom <- as.double(sum1)*sum2*sum3*sum4 # as.double to avoid overflow error on large products
  
  if (any(sum1==0, sum2==0, sum3==0, sum4==0)) {
    denom <- 1
  }
  
  mcc <- ((TP*TN)-(FP*FN)) / sqrt(denom)
  return(mcc)
}

AnalysisCellList<-c(
                     "astrocyte",
#                     "Bergmann glial cell",
                     "central nervous system macrophage",
                     "endothelial cell",
                    "fibroblast",
                     "leukocyte",
                     "neuron",
                     "oligodendrocyte",
                     "oligodendrocyte precursor cell",
                     "pericyte"
                     
)
LegendCellList<-c(
  "Astrocyte",
#  "Bergmann glial cell",
  "CNS macrophage",
  "Endothelial cell",
  "Fibroblast",
  "Leukocyte",
  "Neuron",
  "Oligodendrocyte",
  "Oligo' precursor cell",
  "Pericyte"
  
)
ResultsPerType<-data.frame(LegendCellList)

colors <- c("darkred","darkgreen",1,2,3,4,5,6,7,8,9)
AnalysisCellType<-AnalysisCellList[1]
P1<-plotROC(AnalysisCellType,ATACpredictions,GroundTruth)
plot(P1, 
     main="ROC Curve by Cell Type",lty=1,
     col = colors[1])
ResultsPerType$auROC[1]<-P1$auc


for (i in 2:length(AnalysisCellList)) {
  AnalysisCellType<-AnalysisCellList[i]
  P1<-plotROC(AnalysisCellType,ATACpredictions,GroundTruth)
  lines(P1,col=colors[i])
  ResultsPerType$auROC[i]<-P1$auc

}


legend("bottomright", 
       legend = LegendCellList, 
       col = colors, 
       lty = 1, lwd=3,
       cex = 0.8)

#
# Now do the same for a Precision Recall curve
#
# Install the PRROC package
#install.packages("PRROC")
library(PRROC)

par(pty="s")
par(mar=c(5,4,4,13)+0.1)
i=1
AnalysisCellType<-AnalysisCellList[i]

# Assuming truth and pred are your ground truth labels and prediction scores
truth <- GroundTruth$predicted.id==AnalysisCellType
ACTvbl<- gsub(" ", ".", AnalysisCellType)
evalstr=paste("pred<-ATACpredictions$prediction.score.",ACTvbl,sep="")
eval(parse(text=evalstr))
# Compute PR curve for 1st cell type
pr <- pr.curve(scores.class0 = pred, weights.class0 = truth,curve=TRUE)

# Plot PR curve
plot(pr, 
     main="Precision Recall Curve by Cell Type",lty=1,lwd=2,auc.main=FALSE,
     col = colors[i])
print(paste(AnalysisCellType," AUC:",pr$auc.integral))
binPred<-ATACpredictions$predicted.id==AnalysisCellType
MCC<-mcc(truth,binPred)
print(paste("MCC:",MCC))
ResultsPerType$auPR[i]<-pr$auc.integral
ResultsPerType$MCC[i]<-MCC
ResultsPerType$MeanPredScore[i]<-mean(ATACpredictions$prediction.score.max[binPred])
ResultsPerType$MinPredScore[i]<-min(ATACpredictions$prediction.score.max[binPred])
ResultsPerType$MaxPredScore[i]<-max(ATACpredictions$prediction.score.max[binPred])
CM<-ConfusionMatrix(binPred,truth)
if(ncol(CM)==2) {
  if(nrow(CM)==1) {
    print(paste("NOT SURE FP:",CM[1,2],"TP:",CM[1,1]))
    ResultsPerType$TN[i]<-0
    ResultsPerType$FN[i]<-0
    ResultsPerType$TP[i]<-CM[1,2]
    ResultsPerType$FP[i]<-CM[1,1] 
    
  }
  else {
    print(paste("TN:",CM[1,1],"FN:",CM[2,1],"FP:",CM[1,2],"TP:",CM[2,2]))
    ResultsPerType$TN[i]<-CM[1,1]
    ResultsPerType$FN[i]<-CM[2,1]
    ResultsPerType$FP[i]<-CM[1,2]
    ResultsPerType$TP[i]<-CM[2,2] 
  }
}
if(ncol(CM)==1)  {
  # code for no true positives etc  1 col only!    
  print(paste("TN:",CM[1,1],"FN:",CM[2,1],"FP:",0,"TP:",0))
  ResultsPerType$TN[i]<-CM[1,1]
  ResultsPerType$FN[i]<-CM[2,1]
  ResultsPerType$FP[i]<-0
  ResultsPerType$TP[i]<-0       
  
  
}

#

#
# Now do the rest
#
for (i in 2:length(AnalysisCellList)) {
  AnalysisCellType<-AnalysisCellList[i]
  truth <- GroundTruth$predicted.id==AnalysisCellType
  ACTvbl<- gsub(" ", ".", AnalysisCellType)
  evalstr=paste("pred<-ATACpredictions$prediction.score.",ACTvbl,sep="")
  eval(parse(text=evalstr))
  pr <- pr.curve(scores.class0 = pred, weights.class0 = truth,curve=TRUE)
  lines(pr$curve,col=colors[i],lty=1,lwd=2)
  print(paste(AnalysisCellType," AUC:",pr$auc.integral))
  binPred<-ATACpredictions$predicted.id==AnalysisCellType
  MCC<-mcc(truth,binPred)
  print(paste("MCC:",MCC))
  ResultsPerType$auPR[i]<-pr$auc.integral
  ResultsPerType$MCC[i]<-MCC
  ResultsPerType$MeanPredScore[i]<-mean(ATACpredictions$prediction.score.max[binPred])
  ResultsPerType$MinPredScore[i]<-min(ATACpredictions$prediction.score.max[binPred])
  ResultsPerType$MaxPredScore[i]<-max(ATACpredictions$prediction.score.max[binPred])
  CM<-ConfusionMatrix(binPred,truth)
  if(ncol(CM)==2) {
    if(nrow(CM)==1) {
      print(paste("NOT SURE FP:",CM[1,2],"TP:",CM[1,1]))
      ResultsPerType$TN[i]<-0
      ResultsPerType$FN[i]<-0
      ResultsPerType$TP[i]<-CM[1,2]
      ResultsPerType$FP[i]<-CM[1,1] 
      
    }
    else {
      print(paste("TN:",CM[1,1],"FN:",CM[2,1],"FP:",CM[1,2],"TP:",CM[2,2]))
      ResultsPerType$TN[i]<-CM[1,1]
      ResultsPerType$FN[i]<-CM[2,1]
      ResultsPerType$FP[i]<-CM[1,2]
      ResultsPerType$TP[i]<-CM[2,2] 
    }
  }
  else {
# code for no true positives etc  1 col only!    
      print(paste("TN:",CM[1,1],"FN:",CM[2,1],"FP:",0,"TP:",0))
      ResultsPerType$TN[i]<-CM[1,1]
      ResultsPerType$FN[i]<-CM[2,1]
      ResultsPerType$FP[i]<-0
      ResultsPerType$TP[i]<-0       
      
      
    }
  
    
  }


#Tidy legend

legend("bottomright", 
#       inset = c(-0.4,0),
       legend = LegendCellList, 
       title = "Cell Type Key",
       col = colors, 
       lty = 1, lwd=3,ncol=1,
       cex = 0.8)

#save results in CSV file
write.csv(ResultsPerType,file="CellTypingResults2.csv")
ResultsPerType
