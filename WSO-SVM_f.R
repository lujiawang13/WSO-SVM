WSOSVM_f<-function(C1,C2,TrainData,features,kernelname='linear',Algorithm='quadratic',solver="optiSolve",bmethod="meanb",sigmaYf=0,sigmaYl=1){
  ####################################################### data preprocess
  if((any(colnames(TrainData)=="y") & any(colnames(TrainData)=="yl") & any(colnames(TrainData)=="yr"))==FALSE){
    stop('The training data must have columns y, yl, and yr.')}
  xTrain<-TrainData[,features]
  xTrain_scale<-scale(xTrain);
  xTrain=xTrain_scale
  Data=data.frame(xTrain,y=TrainData[,'y'],yl=TrainData[,'yl'],yr=TrainData[,'yr'])
  Datapreprocess<-Datapreprocess_f(Data)
  xlist=Datapreprocess$xlist
  # ylist=Datapreprocess$ylist
  if(grepl('auto',kernelname)){
    weights<-featureweights_f(C1,C2,xlist)
    sigmaYl<-weights
  }
  
  model<-WSOmodel_f(C1,C2,xlist,kernelname=kernelname,Algorithm=Algorithm,solver=solver,bmethod=bmethod,sigmaYf,sigmaYl)
  model[['xTrain_scale']]=xTrain_scale
  model[['features']]=features
  return(model)
}
WSOmodel_f<-function(C1,C2,xlist,kernelname='mkl',Algorithm='smo',solver="optiSolve",bmethod="meanb",sigmaYf=0,sigmaYl=1){
  ################ data store
  X=do.call("rbind",xlist)
  p=ncol(X)
  n1=nrow(xlist[[1]])
  n2=nrow(xlist[[2]])
  n3=nrow(xlist[[3]])
  n4=nrow(xlist[[4]])
  n12=n1+n2
  n34=n3+n4
  n=n12+n34
  ylist<-list()
  ylist[[1]]=rep(1,n1)
  ylist[[2]]=rep(-1,n2)
  ylist[[3]]=rep(1,n3)
  ylist[[4]]=rep(-1,n4)
  y<-unlist(ylist)
  Info=list(X=X,y=y,ylist=ylist,C1=C1,C2=C2,
            n1=n1,n2=n2,n3=n3,n4=n4,
            n12=n12,n34=n34,n=n,
            sigmaYf=sigmaYf,sigmaYl=sigmaYl)

  if(kernelname=='mkl'){
    Klist<-lapply(1:p,iK_f,X,sigmaYl)
    Info[['Klist']]=Klist
    model<-MKLorderSVM_f(C1,C2,xlist,kernelname,Algorithm,solver,bmethod,Info=Info)
  }else{
    Klist=NULL
    Info[['Klist']]=Klist
    model<-WSOAlgorithm_f(C1,C2,xlist,kernelname,Algorithm,solver,bmethod,Info=Info)
  }
  return(model)
}

WSOAlgorithm_f<-function(C1,C2,xlist,kernelname,Algorithm,solver,bmethod,Info){
  if(Algorithm=='smo3'){
    model=WSOSMO3_f(C1,C2,xlist,kernelname,solver=solver,bmethod,Info=Info)
  }
  if(Algorithm=='smo2'){
    model=WSOSMO2_f(C1,C2,xlist,kernelname,solver=solver,bmethod,Info=Info)
  }
  if(Algorithm=='quadratic'){
    model=WSOquadratic_f(C1,C2,xlist,kernelname,solver,bmethod,Info=Info)
  }
  return(model)
}


WSOquadratic_f<-function(C1,C2,xlist,kernelname,solver,bmethod,Info){
  sigmaYf=Info$sigmaYf
  sigmaYl=Info$sigmaYl
  X=Info$X
  y=Info$y
  ylist=Info$ylist
  Klist=Info$Klist
  n1=Info$n1
  n2=Info$n2
  n3=Info$n3
  n4=Info$n4
  n12=Info$n12
  n34=Info$n34
  n=Info$n
  
  ###################################################
  K<-K_f(X,X,sigmaYf,sigmaYl,kernelname,Klist)+0.0
  y_matrix<-diag(y)
  Dmat=y_matrix%*%K%*%y_matrix+0.0001*diag(nrow(K))
  dvec=rep(1,nrow(Dmat))
  
  ########### solve the dual form
  if(solver=="lpSolve"){
    Amat0=rbind(y,diag(nrow(X)),-diag(nrow(X)),
                c(rep(1,nrow(xlist[[1]])),rep(-1,nrow(xlist[[2]])),rep(0,n34)))
    Amat=t(Amat0)
    bvec=c(0,rep(0,nrow(X)),
           c(-rep(C1,n12),-rep(C2,n34)),
           0)
    solve_info=solve.QP(Dmat, dvec, Amat, bvec, meq=1, factorized=FALSE)
    alphas=solve_info$solution
  }
  if(solver=="quadprog" | solver=="optiSolve"){
    H=Dmat
    c=-dvec
    A=rbind(y,c(rep(1,nrow(xlist[[1]])),rep(-1,nrow(xlist[[2]])),rep(0,n34)))
    b=rep(0,2)
    r=c(0,999999)
    l=rep(0,nrow(X))
    u=c(rep(C1,n12),
        rep(C2,n34))
    solve_info=ipop(c, H, A, b, l, u, r, sigf = 7, maxiter = 40, margin = 0.05,
                    bound = 10, verb = 0)
    alphas=attributes(solve_info)$primal
  }
  
  j=1
  alphas_list<-list()
  for(i in 1:length(xlist)){
    alphas_list[[i]]<-alphas[j:(j+nrow(xlist[[i]])-1)]
    j=j+nrow(xlist[[i]])
    # print(j)
  }
  ############################
  model=list(alphas=as.numeric(alphas),alphas_list=alphas_list,C1=C1,C2=C2,sigmaYf=sigmaYf,sigmaYl=sigmaYl,ylist=ylist,xlist=xlist,
             X=X,y=y,Klist=Klist,kernelname=kernelname)
  b_info<-b_f(model,bmethod)
  b1<-b_info$b1
  b2<-b_info$b2
  model[["b1"]]=b1
  model[["b2"]]=b2
  return(model)
}



########################################### alphas smo3
WSOSMO3_f <- function(C1,C2,xlist,kernelname,
                        solver="optiSolve",bmethod,Info) {
  tol=1e-4
  niter=50
  ###############################
  sigmaYf=Info$sigmaYf
  sigmaYl=Info$sigmaYl
  X=Info$X
  y=Info$y
  ylist=Info$ylist
  Klist=Info$Klist
  n12=Info$n12
  n34=Info$n34
  n=Info$n
  # K<-K_f(X,X,sigmaYf,sigmaYl,kernelname,Klist)+0.0
  
  # variables
  alphas <- rep(0, n)
  
  iter=0
  while (iter < niter) {
    iter=iter+1
    print(paste('alphas',iter))
    alphas0=alphas
    for (i in 1:n){
      # randomly select another two instances
      # tIndx <- c(sample((1:m)[-i],2,replace = FALSE),i)
      if(i<=n12){A1=(1:n12)[-i];A2=(n12+1):m}else{A1=(1:n12);A2=((n12+1):m)[-(i-n12)]}
      j1=sample(A1,1,replace = FALSE)
      j2=sample(A2,1,replace = FALSE)
      tIndx <- c(i,j1,j2)
      t1=sort(tIndx)[1]
      t2=sort(tIndx)[2]
      t3=sort(tIndx)[3]
      
      
      if((all(y[c(t1,t2,t3)]==1) & sum((alphas*y)[-c(t1,t2,t3)])>0) | (all(y[c(t1,t2,t3)]==-1) & sum((alphas*y)[-c(t1,t2,t3)])<0)) next
      if(sum(c(t1,t2,t3)<=n12)>2 | sum(c(t1,t2,t3)>n12)>2) next
      
      ## compute alphas 
      y_matrix<-diag(y[c(t1,t2,t3)])
      if(kernelname=='mkl'){Klist<-lapply(1:p,iK_f,X[c(t1,t2,t3),],sigmaYl)+0.0}else{Klist=NULL}
      Ksub<-K_f(X[c(t1,t2,t3),],X[c(t1,t2,t3),],sigmaYf,sigmaYl,kernelname,Klist)+0.0
      Dmat=y_matrix%*%Ksub%*%y_matrix+0.0001*diag(3)
      kvt1=apply(X,1,kernel_f,X[t1,],sigmaYf,sigmaYl,kernelname)
      kvt2=apply(X,1,kernel_f,X[t2,],sigmaYf,sigmaYl,kernelname)
      kvt3=apply(X,1,kernel_f,X[t3,],sigmaYf,sigmaYl,kernelname)
      dvec=rep(1,nrow(Dmat))-c(y[t1]*sum((alphas*y*kvt1)[-c(t1,t2,t3)]),
                               y[t2]*sum((alphas*y*kvt2)[-c(t1,t2,t3)]),y[t3]*sum((alphas*y*kvt3)[-c(t1,t2,t3)]))
      
      #################################################
      tsleft=sum(c(t1,t2,t3)<=n12)
      if(tsleft>=2){
        orderconstraint=rep(0,3)
        orderconstraint[1:tsleft]=y[c(t1,t2,t3)[1:tsleft]]
        rightvalue=-sum((alphas*y)[-c(t1,t2,t3)])
        if(solver=="lpSolve"){
          Amat0=rbind(y[c(t1,t2,t3)],diag(3),-diag(3),
                      as.numeric(orderconstraint))
          Amat=t(Amat0)
          bvec=c(rightvalue,rep(0,3),
                 c(-rep(C1,tsleft),-rep(C2,3-tsleft)),
                 -sum((alphas[1:n12]*y[1:n12])[-(c(t1,t2,t3)[1:tsleft])]))
          check=has_error(solve.QP(Dmat, dvec, Amat, bvec, meq=1, factorized=FALSE))
          if(check==TRUE) next
          solve_info=solve.QP(Dmat, dvec, Amat, bvec, meq=1, factorized=FALSE)
          alphas3=solve_info$solution
        }
        ###############################################################################
        if(solver=="optiSolve"){
          a=-dvec
          Amat=rbind(y[c(t1,t2,t3)],
                     as.numeric(orderconstraint))
          bvec=c(rightvalue,
                 -sum((alphas[1:n12]*y[1:n12])[-(c(t1,t2,t3)[1:tsleft])]))
          mycop <- cop(f = quadfun(Q=Dmat, a=a, d=0, id=1:3),
                       lb = lbcon(rep(0,3), id=1:3),
                       ub = ubcon(c(rep(C1,tsleft),rep(C2,3-tsleft)), id=1:3),
                       lc =  lincon(A=Amat, d=rep(0, nrow(Amat)), dir=c("==",">="), val=bvec,
                                    id=1:ncol(Amat), use=rep(TRUE,nrow(Amat)), name=1:nrow(Amat)))
          solve_info <- solvecop(mycop, solver="cccp", quiet=TRUE)
          alphas3=solve_info$x
        }
        #######################################################################
        if(solver=="quadprog"){
          H=Dmat
          # c=c(y[t1]*sum((alphas*y*kvt1)[-c(t1,t2,t3)]),y[t2]*sum((alphas*y*kvt2)[-c(t1,t2,t3)]),y[t3]*sum((alphas*y*kvt3)[-c(t1,t2,t3)]))-rep(1,nrow(Dmat))
          c=-dvec
          A=rbind(y[c(t1,t2,t3)],orderconstraint)
          b=c(-sum((alphas*y)[-c(t1,t2,t3)]),-sum((alphas[1:n12]*y[1:n12])[-(c(t1,t2,t3)[1:tsleft])])) #-0.000001
          r=c(0,999999)
          l=rep(0,3)
          u=c(rep(C1,tsleft),
              rep(C2,3-tsleft))
          solve_info=ipop(c, H, A, b, l, u, r, sigf = 7, maxiter = 40, margin = 0.05,
                          bound = 10, verb = 0)
          alphas3check=attributes(solve_info)$primal
        }
      }else{ #tsleft<2
        orderconstraint=rep(0,3)
        orderconstraint[(tsleft+1):3]=-y[c(t1,t2,t3)[(tsleft+1):3]]
        rightvalue=-sum((alphas*y)[-c(t1,t2,t3)])
        if(solver=="lpSolve"){
          Amat0=rbind(y[c(t1,t2,t3)],diag(3),-diag(3),
                      as.numeric(orderconstraint))
          Amat=t(Amat0)
          bvec=c(rightvalue,rep(0,3),
                 c(-rep(C1,tsleft),-rep(C2,3-tsleft)),
                 sum((alphas*y)[-c(1:n12,c(t1,t2,t3)[(tsleft+1):3])]))
          check=has_error(solve.QP(Dmat, dvec, Amat, bvec, meq=1, factorized=FALSE))
          if(check==TRUE) next
          solve_info=solve.QP(Dmat, dvec, Amat, bvec, meq=1, factorized=FALSE)
          alphas3=solve_info$solution
        }
        ###############################################################################
        if(solver=="optiSolve"){
          a=-dvec
          Amat=rbind(y[c(t1,t2,t3)],
                     as.numeric(orderconstraint))
          bvec=c(rightvalue,
                 sum((alphas*y)[-c(1:n12,c(t1,t2,t3)[(tsleft+1):3])]))
          mycop <- cop(f = quadfun(Q=Dmat, a=a, d=0, id=1:3),
                       lb = lbcon(rep(0,3), id=1:3),
                       ub = ubcon(c(rep(C1,tsleft),rep(C2,3-tsleft)), id=1:3),
                       lc =  lincon(A=Amat, d=rep(0, nrow(Amat)), dir=c("==",">="), val=bvec,
                                    id=1:ncol(Amat), use=rep(TRUE,nrow(Amat)), name=1:nrow(Amat)))
          solve_info <- solvecop(mycop, solver="cccp", quiet=TRUE)
          alphas3=solve_info$x
        }
        #######################################################################
        if(solver=="quadprog"){
          H=Dmat
          # c=c(y[i]*sum((alphas*y*kvi)[-c(t1,t2,t3)]),y[j]*sum((alphas*y*kvj)[-c(t1,t2,t3)]),y[t]*sum((alphas*y*kvt)[-c(t1,t2,t3)]))-rep(1,nrow(Dmat))
          c=-dvec
          A=rbind(y[c(t1,t2,t3)],orderconstraint)
          b=c(-sum((alphas*y)[-c(t1,t2,t3)]),sum((alphas[-c(1:n12)]*y[-c(1:n12)])[-(c(t1,t2,t3)[(tsleft+1):3])])-0.000001)
          r=c(0,999999)
          l=rep(0,3)
          u=c(rep(C1,tsleft),
              rep(C2,3-tsleft))
          solve_info=ipop(c, H, A, b, l, u, r, sigf = 7, maxiter = 40, margin = 0.05,
                          bound = 10, verb = 0)
          alphas3check=attributes(solve_info)$primal
        }
      }
      # if(sum(abs(alphas3check-alphas3))>0.001){print("error for check alphas")}
      
      alphas[c(t1,t2,t3)] <- alphas3
      
    }
    if(sum(abs(alphas-alphas0))<tol) break
  } # End while
  
  
  j=1
  alphas_list<-list()
  for(i in 1:length(xlist)){
    alphas_list[[i]]<-alphas[j:(j+nrow(xlist[[i]])-1)]
    j=j+nrow(xlist[[i]])
    # print(j)
  }
  ############################
  model=list(alphas=as.numeric(alphas),alphas_list=alphas_list,C1=C1,C2=C2,sigmaYf=sigmaYf,sigmaYl=sigmaYl,ylist=ylist,xlist=xlist,
             X=X,y=y,Klist=Klist,kernelname=kernelname)
  b_info<-b_f(model,bmethod)
  b1<-b_info$b1
  b2<-b_info$b2
  model[["b1"]]=b1
  model[["b2"]]=b2
  return(model)
}


########################################### alphas smo2
WSOSMO2_f <- function(C1,C2,xlist,kernelname,
                        solver="optiSolve",bmethod,Info) {
  tol=1e-4
  niter=50
  ##########################################
  sigmaYf=Info$sigmaYf
  sigmaYl=Info$sigmaYl
  X=Info$X
  y=Info$y
  ylist=Info$ylist
  Klist=Info$Klist
  n12=Info$n12
  n34=Info$n34
  n=Info$n
  K<-K_f(X,X,sigmaYf,sigmaYl,kernelname,Klist)+0.0
  
  # variables
  alphas <- rep(0, n)
  
  iter=0
  while (iter < niter) {
    iter=iter+1
    print(paste('alphas',iter))
    alphas0=alphas
    for (i in 1:n){
      # randomly select another two instances
      # tIndx <- c(sample((1:m)[-i],2,replace = FALSE),i)
      if(i<=n12){A1=(1:n12)[-i];A2=(n12+1):m}else{A1=(1:n12);A2=((n12+1):m)[-(i-n12)]}
      j1=sample(A1,1,replace = FALSE)
      j2=sample(A2,1,replace = FALSE)
      tIndx <- c(i,j1,j2)
      # tsleft=sum(tIndx<=n12)
      if(i<=n12){t3=j2;t1=i;t2=j1}
      if(i>n12){t3=i;t1=j1;t2=j2}
      
      
      if((all(y[c(t1,t2,t3)]==1) & sum((alphas*y)[-c(t1,t2,t3)])>0) | (all(y[c(t1,t2,t3)]==-1) & sum((alphas*y)[-c(t1,t2,t3)])<0)) next
      if(sum(c(t1,t2,t3)<=n12)>2 | sum(c(t1,t2,t3)>n12)>2) next
      
      ## compute alphas 
      A=K[t1,t1]-2*K[t1,t3]+K[t3,t3]
      B=K[t2,t2]-2*K[t2,t3]+K[t3,t3]
      D=y[t1]*y[t2]*(K[t1,t2]-K[t1,t3]-K[t2,t3]+K[t3,t3])
      alphasy=alphas*y
      alphasysum4=sum(alphasy[-c(t1,t2,t3)])
      
      Dmat=matrix(c(A,D,D,B),nrow=2,ncol=2)+0.000001*diag(2)
      E=y[t1]*(alphasysum4*(K[t3,t3]-K[t1,t3])+sum((alphasy*(K[t1,]-K[t3,]))[-c(t1,t2,t3)])+y[t3])-1
      G=y[t2]*(alphasysum4*(K[t3,t3]-K[t2,t3])+sum((alphasy*(K[t2,]-K[t3,]))[-c(t1,t2,t3)])+y[t3])-1
      dvec=-c(E,G)
      ################################################
      alpha3eq_1=-c(y[t1]*y[t3],y[t2]*y[t3])
      alpha3eq_2=-alpha3eq_1
      alpha3value=alphasysum4*y[t3]
      ################
      if(i<=n12){
        orderconstraint=y[c(t1,t2)]
        orderrightvalue=-sum(alphasy[(1:n12)[-c(t1,t2)]])
        C12=C1
        C3=C2 # the third one is beta not alpha
        if(y[t3]==1 & (-alphasysum4)<orderrightvalue) next
        if(y[t3]==-1 & (C3-alphasysum4)<orderrightvalue) next
      }else{
        orderconstraint=-y[c(t1,t2)]
        orderrightvalue=sum(alphasy[-c(1:n12,t1,t2)])
        C12=C2
        C3=C1
        if(y[t3]==1 & (C3+alphasysum4)<orderrightvalue) next
        if(y[t3]==-1 & (alphasysum4)<orderrightvalue) next
      }
      #########################################
      
      if(solver=="lpSolve"){
        Amat0=rbind(alpha3eq_1,alpha3eq_2,diag(2),-diag(2),
                    orderconstraint)
        Amat=t(Amat0)
        bvec=c(alpha3value,-C3-alpha3value,rep(0,2),
               -rep(C12,2),
               orderrightvalue)
        check=has_error(solve.QP(Dmat, dvec, Amat, bvec, meq=0, factorized=FALSE))
        if(check==TRUE) next
        solve_info=solve.QP(Dmat, dvec, Amat, bvec, meq=0, factorized=FALSE)
        alphas12=solve_info$solution
      }
      ###############################################################################
      if(solver=="optiSolve"){
        a=-dvec
        Amat=rbind(alpha3eq_1,alpha3eq_2,
                   orderconstraint)
        bvec=c(alpha3value,-C3-alpha3value,
               orderrightvalue)
        mycop <- cop(f = quadfun(Q=Dmat, a=a, d=0, id=1:2),
                     lb = lbcon(rep(0,2), id=1:2),
                     ub = ubcon(rep(C12,2), id=1:2),
                     lc =  lincon(A=Amat, d=rep(0, nrow(Amat)), dir=rep(">=",3), val=bvec,
                                  id=1:ncol(Amat), use=rep(TRUE,nrow(Amat)), name=1:nrow(Amat)))
        solve_info <- solvecop(mycop, solver="cccp", quiet=TRUE)
        if(solve_info$status=="unknown") next
        alphas12=solve_info$x
      }
      #######################################################################
      if(solver=="quadprog"){
        H=Dmat
        # c=c(y[t1]*sum((alphas*y*kvt1)[-c(t1,t2,t3)]),y[t2]*sum((alphas*y*kvt2)[-c(t1,t2,t3)]),y[t3]*sum((alphas*y*kvt3)[-c(t1,t2,t3)]))-rep(1,nrow(Dmat))
        c=-dvec
        A=rbind(alpha3eq_1,orderconstraint)
        b=c(alpha3value,orderrightvalue)
        r=c(C3+alpha3value,99999)
        l=rep(0,2)
        u=rep(C12,2)
        solve_info=ipop(c, H, A, b, l, u, r, sigf = 7, maxiter = 40, margin = 0.05,
                        bound = 10, verb = 0)
        alphas12=attributes(solve_info)$primal
      }
      
      # if(sum(abs(alphas3check-alphas3))>0.001){print("error for check alphas")}
      alphas3=(-sum(alphas12*y[c(t1,t2)])-alphasysum4)*y[t3]
      
      alphas[c(t1,t2)] <- alphas12
      alphas[t3] <- alphas3
      
    } # end i in 1:n
    if(sum(abs(alphas-alphas0))<tol) break
  } # End while
  
  
  j=1
  alphas_list<-list()
  for(i in 1:length(xlist)){
    alphas_list[[i]]<-alphas[j:(j+nrow(xlist[[i]])-1)]
    j=j+nrow(xlist[[i]])
    # print(j)
  }
  ############################
  model=list(alphas=as.numeric(alphas),alphas_list=alphas_list,C1=C1,C2=C2,sigmaYf=sigmaYf,sigmaYl=sigmaYl,ylist=ylist,xlist=xlist,
             X=X,y=y,Klist=Klist,kernelname=kernelname)
  b_info<-b_f(model,bmethod)
  b1<-b_info$b1
  b2<-b_info$b2
  model[["b1"]]=b1
  model[["b2"]]=b2
  return(model)
}



iK_f<-function(i,X,sigmaYl){
  m=nrow(X)
  V=X[,i]
  Xm<-replicate(m, V)
  iK<-exp(-(Xm-t(Xm))^2/sigmaYl/2)+0.0
  return(iK)
}
alphasyK_f<-function(K,alphasy){
  alphasyK=alphasy%*%K%*%alphasy
  return(alphasyK)
}
# wMKL_f<-function(model,alphas,sigmaYl){
#   X=model$X
#   y=model$y
#   Klist=model$Klist
#   p=ncol(X)
#   alphasy=alphas*y
#   objf<-lapply(Klist, alphasyK_f,alphasy)
#   objective.in<-unlist(objf)
#   const.mat=rbind(rep(1,p),diag(p))
#   const.dir=c("==",rep(">=",p))
#   const.rhs=c(1,rep(0,p))
#   solve_info <- lp (direction = "min", objective.in=objective.in, 
#                     const.mat=const.mat, const.dir=const.dir, const.rhs=const.rhs)
#   sigmaYf=solve_info$solution
#   return(list(sigmaYf=sigmaYf))
# }

# Addconstraint_f<-function(X,y,Klist,model){
#   alphas=model$alphas
#   p=ncol(X)
#   alphasy=alphas*y
#   objf<-lapply(Klist, alphasyK_f,alphasy)
#   s<-unlist(objf)
#   Addconst.mat<-c(s,-1)
#   Addconst.dir<-"<="
#   Addconst.rhs<-0
#   return(list(Addconst.mat=Addconst.mat,Addconst.dir=Addconst.dir,
#               Addconst.rhs=Addconst.rhs))
# }
# 
# 
# MKLorderSVM_f<-function(C1,C2,xlist,sigmaYf,sigmaYl,kernelname,Algorithm,solver,Info,
#                         tol=1e-4,niter=50){
#   X=Info$X
#   y=Info$y
#   Klist=Info$Klist
#   p=ncol(X)
#   sigmaYf=rep(1/p,p)
#   sigmaYf0=sigmaYf
#   #####################
#   objective.in<-c(rep(0,p),1)
#   const.mat=rbind(c(rep(1,p),0),cbind(diag(p),rep(0,p)))
#   const.dir=c("==",rep(">=",p))
#   const.rhs=c(1,rep(0,p))
#   #######################
#   iter=0
#   while (iter < niter) {
#     iter=iter+1
#     print(paste("MKL",iter))
#     model<-WSOAlgorithm_f(C1,C2,xlist,sigmaYf,sigmaYl,kernelname,Algorithm,solver,Info)
#     Addconst_info<-Addconstraint_f(X,y,Klist,model)
#     Addconst.mat=Addconst_info$Addconst.mat
#     Addconst.dir=Addconst_info$Addconst.dir
#     Addconst.rhs=Addconst_info$Addconst.rhs
#     const.mat<-rbind(const.mat,Addconst.mat)
#     const.dir<-c(const.dir,Addconst.dir)
#     const.rhs<-c(const.rhs,Addconst.rhs)
#     solve_info <- lp (direction = "min", objective.in=objective.in,
#                       const.mat=const.mat, const.dir=const.dir, const.rhs=const.rhs)
#     sigmaYf=solve_info$solution
#     sigmaYf=sigmaYf[-length(sigmaYf)]
#     # if(abs(1-mean(sigmaYf/sigmaYf0))<tol) break
#     if(sum(abs(sigmaYf-sigmaYf0))<tol) break
#     sigmaYf0=sigmaYf
#   }
#   model[["wMKL"]]=sigmaYf
#   return(model)
# }

###########################################################
Addconstraint_f<-function(X,y,Klist,model){
  alphas=model$alphas
  p=ncol(X)
  alphasy=alphas*y
  objf<-lapply(Klist, alphasyK_f,alphasy)
  s<-unlist(objf)
  Addconst.mat<-c(s,-1)
  Addconst.dir<-">="
  Addconst.rhs<-0
  return(list(Addconst.mat=Addconst.mat,Addconst.dir=Addconst.dir,
              Addconst.rhs=Addconst.rhs))
}


MKLorderSVM_f<-function(C1,C2,xlist,kernelname,Algorithm,solver,bmethod,Info){
  tol=1e-4
  niter=1000
  ############################
  sigmaYf=Info$sigmaYf
  sigmaYl=Info$sigmaYl
  X=Info$X
  y=Info$y
  Klist=Info$Klist
  p=ncol(X)
  sigmaYf=rep(1/p,p)
  #####################
  objective.in<-c(rep(0,p),1)
  const.mat=rbind(c(rep(1,p),0),cbind(diag(p),rep(0,p)))
  const.dir=c("==",rep(">=",p))
  const.rhs=c(1,rep(0,p))
  #######################
  iter=0
  while (iter < niter) {
    iter=iter+1
    print(paste("MKL",iter))
    sigmaYf0=sigmaYf
    model<-WSOAlgorithm_f(C1,C2,xlist,kernelname,Algorithm,solver,bmethod,Info)
    Addconst_info<-Addconstraint_f(X,y,Klist,model)
    Addconst.mat=Addconst_info$Addconst.mat
    Addconst.dir=Addconst_info$Addconst.dir
    Addconst.rhs=Addconst_info$Addconst.rhs
    const.mat<-rbind(const.mat,Addconst.mat)
    const.dir<-c(const.dir,Addconst.dir)
    const.rhs<-c(const.rhs,Addconst.rhs)
    solve_info <- lp (direction = "max", objective.in=objective.in,
                      const.mat=const.mat, const.dir=const.dir, const.rhs=const.rhs)
    sigmaYf=solve_info$solution
    sigmaYf=sigmaYf[-length(sigmaYf)]
    # if(abs(1-mean(sigmaYf/sigmaYf0))<tol) break
    if(sum(abs(sigmaYf-sigmaYf0))<tol) break
  }
  model[["sigmaYf"]]=sigmaYf
  # print(sigmaYf)
  return(model)
}


Datapreprocess_f<-function(Data){
  y=Data[,'y']
  yl=Data[,'yl']
  yr=Data[,'yr']
  x=Data[,!(names(Data) %in% c('y','yl','yr'))]
  K=max(y,na.rm = TRUE)
  
  xlist<-list()
  ylist<-list()
  for(k in 1:(K-1)){
    indx_left<-which(yr==k)
    indx_right<-which(yl==(k+1))
    
    xlist[[k*2-1]]<-x[indx_left,]
    xlist[[k*2]]<-x[indx_right,]
    ylist[[k*2-1]]<-rep(1,length(indx_left))
    ylist[[k*2]]<-rep(-1,length(indx_right))
  }
  
  return(list(xlist=xlist,ylist=ylist))
}

featureweights_f<-function(C1,C2,xlist){
  kernelname='linear'
  Algorithm='quadratic'
  solver="lpSolve"
  # solver="optiSolve"
  # bmethod='maxminb'
  bmethod='meanb'
  model<-WSOmodel_f(C1,C2,xlist,kernelname=kernelname,Algorithm=Algorithm,solver=solver,bmethod=bmethod,
                              sigmaYf=0,sigmaYl=1)
  weights=wlinear_f(model)
  
  return(weights)
}
# calculate w
wilinearInd_f<-function(i,alphas,yTrain,xTrain){
  wilinear<-alphas[i]*yTrain[i]*xTrain[i,]
  return(wilinear)
}
wlinear_f<-function(model){
  xTrain=model$X
  yTrain<-model$y
  alphas<-model$alphas
  wlinear_list=lapply(1:nrow(xTrain),wilinearInd_f,alphas,yTrain,xTrain)
  wlinear=Reduce("+", wlinear_list)  
  
  return(wlinear)
}
