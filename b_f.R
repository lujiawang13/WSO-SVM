b_f<-function(model,bmethod){
  if(bmethod=='maxminb'){
    b_info<-bmaxmin_f(model)
  }
  if(bmethod=='meanb'){
    b_info<-bmean_f(model)
  }
  return(b_info)
}

# #################################
brelated0_f<-function(indx,xTrain0,yTrain0,
                      xTrain,yTrain,sigmaYf,sigmaYl,kernelname,alphas){
  x0=xTrain0[indx,]
  y0=yTrain0[indx]
  kv0<-apply(xTrain,1,kernel_f,x0,sigmaYf,sigmaYl,kernelname)
  brelated<-sum(alphas*yTrain*kv0)
  return(brelated)
}
# 
# bmean_f<-function(model){
#   xlist<-model$xlist
#   ylist<-model$ylist
#   alphaslist<-model$alphas_list
#   yTrain<-model$y
#   xTrain=model$X
#   alphas<-model$alphas
#   sigmaYf<-model$sigmaYf
#   sigmaYl<-model$sigmaYl
#   kernelname<-model$kernelname
#   #
#   xTrain0<-do.call("rbind",xlist[1:2])
#   yTrain0<-unlist(ylist[1:2])
#   alphas0<-unlist(alphaslist[1:2])
#   indices=which(alphas0>0 & alphas0<model$C1)
#   # b1_chain<-c()
#   # for(indx in indices){
#   #   y0=yTrain0[indx]
#   #   x0=xTrain0[indx,]
#   #   kv0<-apply(xTrain,1,kernel_f,x0,sigmaYf,sigmaYl,kernelname)
#   #   b1=y0-sum(alphas*yTrain*kv0)
#   #   b1_chain<-c(b1_chain,b1)
#   # }
#   b1_chain<-sapply(indices,brelated0_f,xTrain0,yTrain0,
#                    xTrain,yTrain,sigmaYf,sigmaYl,kernelname,alphas)
#   b1=mean(yTrain0[indices]-b1_chain)
#   ###############################################################
#   xTrain0<-do.call("rbind",xlist[3:4])
#   yTrain0<-unlist(ylist[3:4])
#   alphas0<-unlist(alphaslist[3:4])
#   indices=which(alphas0>0 & alphas0<model$C2)
#   # b2_chain<-c()
#   # for(indx in indices){
#   #   y0=yTrain0[indx]
#   #   x0=xTrain0[indx,]
#   #   kv0<-apply(xTrain,1,kernel_f,x0,sigmaYf,sigmaYl,kernelname)
#   #   b2=y0-sum(alphas*yTrain*kv0)
#   #   b2_chain<-c(b2_chain,b2)
#   # }
#   b2_chain<-sapply(indices,brelated0_f,xTrain0,yTrain0,
#                    xTrain,yTrain,sigmaYf,sigmaYl,kernelname,alphas)
#   b2=mean(yTrain0[indices]-b2_chain)
#   return(list(b1=b1,b2=b2))
# }

################################
brelated_f<-function(indx,
                     xTrain,yTrain,sigmaYf,sigmaYl,kernelname,alphas){
  x0=xTrain[indx,]
  y0=yTrain[indx]
  kv0<-apply(xTrain,1,kernel_f,x0,sigmaYf,sigmaYl,kernelname)
  brelated<-sum(alphas*yTrain*kv0)
  return(brelated)
}
bmaxmin_f<-function(model){
  xlist<-model$xlist
  ylist<-model$ylist
  alphaslist<-model$alphas_list
  xTrain=model$X
  yTrain<-model$y
  alphas<-model$alphas
  sigmaYf<-model$sigmaYf
  sigmaYl<-model$sigmaYl
  kernelname<-model$kernelname
  #
  xTrain0<-do.call("rbind",xlist[1:2])
  yTrain0<-unlist(ylist[1:2])
  alphas0<-unlist(alphaslist[1:2])
  b1_chain<-sapply(1:length(yTrain0),brelated0_f,xTrain0,yTrain0,
                   xTrain,yTrain,sigmaYf,sigmaYl,kernelname,alphas)
  b1=-(max(b1_chain[yTrain0==-1])+min(b1_chain[yTrain0==1]))/2
  #
  xTrain0<-do.call("rbind",xlist[3:4])
  yTrain0<-unlist(ylist[3:4])
  alphas0<-unlist(alphaslist[3:4])
  indices=which(alphas0>0 & alphas0<model$C2)
  b2_chain<-sapply(1:length(yTrain0),brelated0_f,xTrain0,yTrain0,
                   xTrain,yTrain,sigmaYf,sigmaYl,kernelname,alphas)
  b2=-(max(b2_chain[yTrain0==-1])+min(b2_chain[yTrain0==1]))/2
  # b2=-(max(b2_chain[yTrain0==1])+min(b2_chain[yTrain0==-1]))/2
  return(list(b1=b1,b2=b2))
}

bmean_f<-function(model){
  xlist<-model$xlist
  ylist<-model$ylist
  alphaslist<-model$alphas_list
  yTrain<-model$y
  xTrain=model$X
  alphas<-model$alphas
  sigmaYf<-model$sigmaYf
  sigmaYl<-model$sigmaYl
  kernelname<-model$kernelname
  #
  alphaslist0<-list()
  alphaslist0[[1]]<-rep(0,length(alphaslist[[1]]))
  alphaslist0[[2]]<-rep(0,length(alphaslist[[2]]))
  alphaslist0[[3]]<-rep(0,length(alphaslist[[3]]))
  alphaslist0[[4]]<-rep(0,length(alphaslist[[4]]))
  #
  alphaslist1<- alphaslist0
  alphaslist1[[1]]<- alphaslist[[1]]
  alphaslist2<- alphaslist0
  alphaslist2[[2]]<- alphaslist[[2]]
  #
  alphas1<-unlist(alphaslist1)
  alphas2<-unlist(alphaslist2)
  indices=c(which(alphas1>0 & alphas1<model$C1[1]),which(alphas2>0 & alphas2<model$C1[2]))
  b1_chain<-sapply(indices,brelated_f,
                   xTrain,yTrain,sigmaYf,sigmaYl,kernelname,alphas)
  b1=mean(yTrain[indices]-b1_chain)
  ###############################################################
  alphaslist0<-alphaslist
  alphaslist0[[1]]<-rep(0,length(alphaslist[[1]]))
  alphaslist0[[2]]<-rep(0,length(alphaslist[[2]]))
  alphas0<-unlist(alphaslist0)
  indices=which(alphas0>0 & alphas0<model$C2)
  b2_chain<-sapply(indices,brelated_f,
                   xTrain,yTrain,sigmaYf,sigmaYl,kernelname,alphas)
  b2=mean(yTrain[indices]-b2_chain)
  return(list(b1=b1,b2=b2))
}
