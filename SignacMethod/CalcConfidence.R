GroundTruth3<-read.csv("RNAObjmeta3.csv")
ATACpredictions3<-read.csv("ATACObjPred3.csv")
GroundTruth2<-read.csv("RNAObjmeta2.csv")
ATACpredictions2<-read.csv("ATACObjPred2.csv")
GroundTruth1<-read.csv("RNAObjmeta.csv")
ATACpredictions1<-read.csv("ATACObjPred.csv")

AnalysisCellList<-c(
  "astrocyte",
  "Bergmann glial cell",
  "central nervous system macrophage",
  "endothelial cell",
  "fibroblast",
  "leukocyte",
  "neuron",
  "oligodendrocyte",
  "oligodendrocyte precursor cell",
  "pericyte"
  )
ConfPerType<-data.frame(AnalysisCellList)  
# Assuming truth and pred are your ground truth labels and prediction scores
i=1

for (i in 1:length(AnalysisCellList)) {
  CellType=AnalysisCellList[i]
  relevant<-ATACpredictions1$predicted.id==CellType
  y_prob<-0
  if(sum(relevant==TRUE)>0) {
     newdata<-ATACpredictions1[relevant,]
    ACTvbl<-gsub(" ", ".", CellType)
    evalstr=paste("y_prob<-newdata$prediction.score.",ACTvbl,sep="")
    eval(parse(text=evalstr))
    ystr<-paste(CellType,"mean:",mean(y_prob), "max:",max(y_prob),"min:",min(y_prob),"median:",median(y_prob))
    print(ystr)
    
  
  ConfPerType$mean[i]<-mean(y_prob)
  ConfPerType$max[i]<-max(y_prob)
  ConfPerType$min[i]<-min(y_prob)
  ConfPerType$median[i]<-median(y_prob)
  }
  else {
    ConfPerType$mean[i]<-0
    ConfPerType$max[i]<-0
    ConfPerType$min[i]<-0
    ConfPerType$median[i]<-0  
    
  }
}
write.csv(ConfPerType,"ATACPredConf1.csv")

i=1

for (i in 1:length(AnalysisCellList)) {
  CellType=AnalysisCellList[i]
  relevant<-GroundTruth1$predicted.id==CellType
  y_prob<-0
  if(sum(relevant==TRUE)>0) {
    newdata<-GroundTruth1[relevant,]
    ACTvbl<-gsub(" ", ".", CellType)
    evalstr=paste("y_prob<-newdata$prediction.score.",ACTvbl,sep="")
    eval(parse(text=evalstr))
    ystr<-paste(CellType,"mean:",mean(y_prob), "max:",max(y_prob),"min:",min(y_prob),"median:",median(y_prob))
    print(ystr)
    
    
    ConfPerType$mean[i]<-mean(y_prob)
    ConfPerType$max[i]<-max(y_prob)
    ConfPerType$min[i]<-min(y_prob)
    ConfPerType$median[i]<-median(y_prob)
  }
  else {
    ConfPerType$mean[i]<-0
    ConfPerType$max[i]<-0
    ConfPerType$min[i]<-0
    ConfPerType$median[i]<-0  
    
  }
}
write.csv(ConfPerType,"GTPredConf1.csv")
i=1

for (i in 1:length(AnalysisCellList)) {
  CellType=AnalysisCellList[i]
  relevant<-ATACpredictions2$predicted.id==CellType
  y_prob<-0
  if(sum(relevant==TRUE)>0) {
    newdata<-ATACpredictions2[relevant,]
    ACTvbl<-gsub(" ", ".", CellType)
    evalstr=paste("y_prob<-newdata$prediction.score.",ACTvbl,sep="")
    eval(parse(text=evalstr))
    ystr<-paste(CellType,"mean:",mean(y_prob), "max:",max(y_prob),"min:",min(y_prob),"median:",median(y_prob))
    print(ystr)
    
    
    ConfPerType$mean[i]<-mean(y_prob)
    ConfPerType$max[i]<-max(y_prob)
    ConfPerType$min[i]<-min(y_prob)
    ConfPerType$median[i]<-median(y_prob)
  }
  else {
    ConfPerType$mean[i]<-0
    ConfPerType$max[i]<-0
    ConfPerType$min[i]<-0
    ConfPerType$median[i]<-0  
    
  }
}
write.csv(ConfPerType,"ATACPredConf2.csv")

i=1

for (i in 1:length(AnalysisCellList)) {
  CellType=AnalysisCellList[i]
  relevant<-GroundTruth2$predicted.id==CellType
  y_prob<-0
  if(sum(relevant==TRUE)>0) {
    newdata<-GroundTruth2[relevant,]
    ACTvbl<-gsub(" ", ".", CellType)
    evalstr=paste("y_prob<-newdata$prediction.score.",ACTvbl,sep="")
    eval(parse(text=evalstr))
    ystr<-paste(CellType,"mean:",mean(y_prob), "max:",max(y_prob),"min:",min(y_prob),"median:",median(y_prob))
    print(ystr)
    
    
    ConfPerType$mean[i]<-mean(y_prob)
    ConfPerType$max[i]<-max(y_prob)
    ConfPerType$min[i]<-min(y_prob)
    ConfPerType$median[i]<-median(y_prob)
  }
  else {
    ConfPerType$mean[i]<-0
    ConfPerType$max[i]<-0
    ConfPerType$min[i]<-0
    ConfPerType$median[i]<-0  
    
  }
}
write.csv(ConfPerType,"GTPredConf2.csv")

i=1

for (i in 1:length(AnalysisCellList)) {
  CellType=AnalysisCellList[i]
  relevant<-ATACpredictions3$predicted.id==CellType
  y_prob<-0
  if(sum(relevant==TRUE)>0) {
    newdata<-ATACpredictions3[relevant,]
    ACTvbl<-gsub(" ", ".", CellType)
    evalstr=paste("y_prob<-newdata$prediction.score.",ACTvbl,sep="")
    eval(parse(text=evalstr))
    ystr<-paste(CellType,"mean:",mean(y_prob), "max:",max(y_prob),"min:",min(y_prob),"median:",median(y_prob))
    print(ystr)
    
    
    ConfPerType$mean[i]<-mean(y_prob)
    ConfPerType$max[i]<-max(y_prob)
    ConfPerType$min[i]<-min(y_prob)
    ConfPerType$median[i]<-median(y_prob)
  }
  else {
    ConfPerType$mean[i]<-0
    ConfPerType$max[i]<-0
    ConfPerType$min[i]<-0
    ConfPerType$median[i]<-0  
    
  }
}
write.csv(ConfPerType,"ATACPredConf3.csv")

i=1

for (i in 1:length(AnalysisCellList)) {
  CellType=AnalysisCellList[i]
  relevant<-GroundTruth3$predicted.id==CellType
  y_prob<-0
  if(sum(relevant==TRUE)>0) {
    newdata<-GroundTruth3[relevant,]
    ACTvbl<-gsub(" ", ".", CellType)
    evalstr=paste("y_prob<-newdata$prediction.score.",ACTvbl,sep="")
    eval(parse(text=evalstr))
    ystr<-paste(CellType,"mean:",mean(y_prob), "max:",max(y_prob),"min:",min(y_prob),"median:",median(y_prob))
    print(ystr)
    
    
    ConfPerType$mean[i]<-mean(y_prob)
    ConfPerType$max[i]<-max(y_prob)
    ConfPerType$min[i]<-min(y_prob)
    ConfPerType$median[i]<-median(y_prob)
  }
  else {
    ConfPerType$mean[i]<-0
    ConfPerType$max[i]<-0
    ConfPerType$min[i]<-0
    ConfPerType$median[i]<-0  
    
  }
}
write.csv(ConfPerType,"GTPredConf3.csv")
