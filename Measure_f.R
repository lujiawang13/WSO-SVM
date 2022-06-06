MSEi_f<-function(i,predictedLabels,yl,yr){
  predictedLabel=predictedLabels[i]
  yli=yl[i]
  yri=yr[i]
  if(predictedLabel<yli){
    MSEi=yli-predictedLabel
  }else if(predictedLabel>yri){
    MSEi=predictedLabel-yri
  }else{
    MSEi=0
  }
  return(MSEi)
}

MSE_f<-function(predictedLabels,Data){
  yl=Data[,'yl']
  yr=Data[,'yr']
  MSE<-sapply(1:length(predictedLabels),MSEi_f,predictedLabels,yl,yr)
  MSE<-mean(MSE)
  return(MSE)
}

loss01i_f<-function(i,predictedLabels,yl,yr){
  predictedLabel=predictedLabels[i]
  yli=yl[i]
  yri=yr[i]
  if(predictedLabel<yli | predictedLabel>yri){
    loss01i=1
  }else{
    loss01i=0
  }
  return(loss01i)
}

loss01_f<-function(predictedLabels,Data){
  yl=Data[,'yl']
  yr=Data[,'yr']
  loss01<-sapply(1:length(predictedLabels),loss01i_f,predictedLabels,yl,yr)
  loss01=mean(loss01)
  return(loss01)
}

diagAcc<-function(predictedLabels,Data,K){
  y=Data[,'y']
  predictedLabels
  if(length(levels(factor(y)))<K){y=factor(y, levels =1:K)}
  if(length(levels(factor(predictedLabels)))<K){predictedLabels=factor(predictedLabels, levels = 1:K)}
  diagAcc=mean(diag(table(y=y,predictedLabels))/table(y))
  return(diagAcc)
}