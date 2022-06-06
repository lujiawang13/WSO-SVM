Predict_f<-function(TestData,model){
  pred_info <- OnePredict_f(TestData,model)
  Labels<-LabelsPredict_f(pred_info)
  return(Labels)
}
OnePredict_f<-function(TestData,model){
  features=model$features
  xTrain_scale=model$xTrain_scale
  xTest<- TestData[,features]
  xTest<- scale(xTest, attr(xTrain_scale, "scaled:center"), attr(xTrain_scale, "scaled:scale"))
  
  b1=model$b1
  b2=model$b2
  wx<-wx_f(xTest,model)
  predY1=wx+b1
  predY2=wx+b2
  return(list(predY1=predY1,predY2=predY2,wx=wx,b_chain=c(b1,b2)))
}
# W*x in prediction
wxInd_f<-function(xTestInd,model){
  # xlist=model$xlist
  xTrain=model$X
  # ylist<-model$ylist
  yTrain<-model$y
  alphas<-model$alphas
  sigmaYf<-model$sigmaYf
  sigmaYl<-model$sigmaYl
  kernelname<-model$kernelname
  kv<-apply(xTrain,1,kernel_f,xTestInd,sigmaYf,sigmaYl,kernelname)
  wx<-sum(alphas*yTrain*kv)
  return(wx)
}
wx_f<-function(xTest,model){
  wx=apply(xTest,1,wxInd_f,model)
  return(wx)
}
##################### Label identification
OneLabel_f<-function(wxi,b_chain){
  if(any((wxi+b_chain)>=0)){
    Label<-which((wxi+b_chain)>=0)[1]
  }else{
    Label<-length(b_chain)+1
  }
  return(Label)
}
LabelsPredict_f<-function(pred_info){
  b_chain=pred_info$b_chain
  wx=pred_info$wx
  Labels<-sapply(wx, OneLabel_f,b_chain)
  return(Labels)
}