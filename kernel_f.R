kernel_f<-function(Xnew,X,sigmaYf=NULL,sigmaYl=NULL,kernelname){
  Xnew=as.numeric(Xnew)
  X=as.numeric(X)
  if(kernelname=='mkl' | kernelname=='fixedmkl'){
    kernel=sum(sigmaYf*exp(-(Xnew-X)^2/sigmaYl/2))
  }else if(kernelname=='ANOVAmkl'){
    k=3
    d=1
    kernel=sum(sigmaYf*exp(-(Xnew^k-X^k)^2/sigmaYl/2)^d)
  }else if(kernelname=='meanmkl'){
    kernel=mean(exp(-(Xnew-X)^2/sigmaYl/2))
  }else if(kernelname=='linear'){
    kernel=X%*%Xnew
    if(length(kernel)!=1){print('error in X*Xnew')}
  }else if(kernelname=='meanlinear'){
    # kernel=X%*%Xnew/length(X)
    kernel=mean(X*Xnew)
  }else if(kernelname=='meanabsdiff'){
    kernel=1/mean(abs(Xnew-X))
  }else if(kernelname=='rbf'){
    kernel=exp(-mean((Xnew-X)^2)/sigmaYl/2)
  }else if(kernelname=='laplace'){
    kernel=exp(-mean(abs(Xnew-X))/sigmaYl)
  }else if(kernelname=='wrbf'){
    if(length(sigmaYl)<=1){print('error in wrbf weights')}
    kernel=exp(-sum(sigmaYl*(Xnew-X)^2))
  }else if(kernelname=='wabslaplace'){
    if(length(sigmaYl)<=1){print('error in wabslaplace weights')}
    kernel=exp(-sum(sigmaYl*abs(Xnew-X)))
  }else if(kernelname=='wlaplace'){
    if(length(sigmaYl)<=1){print('error in wlaplace weights')}
    kernel=exp(-abs(sum(sigmaYl*(Xnew-X))))
  }else if(grepl('fixedwsquaredrbf',kernelname)){ # not work
    if(length(sigmaYl)<=1){print('error in fixedwsquaredrbf weights')}
    kernel=exp(-sum(sigmaYl*(Xnew-X)^2))
  }else if(grepl('fixedwrbf',kernelname)){
    if(length(sigmaYl)<=1){print('error in fixedwrbf weights')}
    kernel=exp(-(sum(sigmaYl*(Xnew-X)))^2)
  }else if(grepl('fixedwabslaplace',kernelname)){ # not work
    if(length(sigmaYl)<=1){print('error in fixedwabslaplace weights')}
    kernel=exp(-sum(sigmaYl*abs(Xnew-X)))
  }else if(grepl('fixedwlaplace',kernelname)){
    if(length(sigmaYl)<=1){print('error in fixedwlaplace weights')}
    kernel=exp(-abs(sum(sigmaYl*(Xnew-X))))
  }else if(kernelname=='Rational-Quadratic'){
    alpha=1
    # alpha=2
    kernel=sigmaYf*(1+mean((Xnew-X)^2)/(alpha*sigmaYl*2))^(-alpha)
  }else if(grepl('Matern',kernelname)){
    sigmaYf=1
    d=mean(abs(Xnew-X))
    # d=sum(abs(Xnew-X))
    # d=mean((Xnew-X)^2)
    if(kernelname=='Matern12'){
      kernel=sigmaYf^2*exp(-d/sigmaYl)
    }else if(kernelname=='Matern32'){
      kernel=sigmaYf^2*(1+sqrt(3)*d/sigmaYl)*exp(-sqrt(3)*d/sigmaYl)
    }else if(kernelname=='Matern52'){
      kernel=sigmaYf^2*(1+sqrt(5)*d/sigmaYl+5*d^2/(3*sigmaYl^2))*exp(-sqrt(5)*d/sigmaYl)
    }
  }else if(kernelname=='xy2'){
    d=2
    kernel=(sum(Xnew*X)+1)^d
  }else if(kernelname=='xy3'){
    d=3
    kernel=(sum(Xnew*X)+1)^d
  }else if(kernelname=='xy4'){
    d=4
    kernel=(sum(Xnew*X)+1)^d
  }else if(kernelname=='sigmoid'){
    kernel=tanh(sum(Xnew*X)*sigmaYl+sigmaYf)
  }
  return(kernel)
}

K_f <- function(Xnew, X,sigmaYf,sigmaYl,kernelname,Klist=NULL){
  if(kernelname=='mkl'){
    p=ncol(X)
    K=0
    for(i in 1:p){
      K<-K+Klist[[i]]*sigmaYf[i]
    }
  }else if(kernelname=='linear'){
    K=as.matrix(X)%*%as.matrix(t(Xnew))
  }else{
    m <- nrow(Xnew)
    n <- nrow(X)
    K <- matrix(NA,nrow=m,ncol=n)
    kernelCol_f<-function(xj,Xnew, X,sigmaYf,sigmaYl,kernelname){
      z=apply(Xnew,1,kernel_f,X[xj,],sigmaYf,sigmaYl,kernelname)
      return(z)
    }
    K=sapply(1:n,kernelCol_f,Xnew, X,sigmaYf,sigmaYl,kernelname)
  }
  return(K)
}
