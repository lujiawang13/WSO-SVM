######   Solve Weakly Supervised Ordinal Support Vector Machine (WSO-SVM) (Wang et al. 2022)
#   Training Input:
#   C1: tuning parameter for 1 vs (2,3).
#   C2: tuning parameter for (1,2) vs 3.
#   TrainData: Training data. It must have columns y, yl, and yr, which y is the precise label, and it is NA if imprecise; and yl/yr are the lower/upper bounds of the interval labels.
#   features: predictors.
#   Algorithm: To solve the primal WSO-SVM, including Algorithm='quadratic' (solving dual form), and 'smo3', 'smo2'.
#   solver: To solve the quadratic programming in dual WSO-SVM, including solver="lpSolve", "optiSolve" and "quadprog".
#   bmethod: To calculate the ordinal parameters b1 and b2 (b1<=b2). bmethod='meanb' or 'maxminb'.
#   kernelname: 'linear', 'laplace', 'rbf', polynomial kernels, Matern kernels, 'fixedwlaplace_auto', 'fixedwrbf_auto',...

#   Testing Input:
#   TestData: Testing data.

#   Training Output:
#   model: trained models.

#   Testing Output:
#   Labels: Predicted labels.

#   Measurement Out:
#   MAE: Mean Absolute Error.
#   loss01: Average of (y_pred!=y_true), which is equivalent to 1-Accuracy.
#   diagAcc: Average accuracy for each class.


# Written by Lujia Wang @ Georgia Tech

# library(quadprog)
# library(lpSolve)
# library(optiSolve)
rm(list=ls(all=TRUE));
# path of the code functions
CodePath="C:\\"
source(paste(CodePath,"kernel_f.R",sep=""))
source(paste(CodePath,"WSO-SVM_f.R",sep=""))
source(paste(CodePath,"b_f.R",sep=""))
source(paste(CodePath,"Predict_f.R",sep=""))
source(paste(CodePath,"Measure_f.R",sep=""))

# path of training data
Datapath=CodePath
TrainData=read.csv(paste(Datapath,'TrainData.csv',sep=''))
TestData=read.csv(paste(Datapath,'TestData.csv',sep=''))
features<-c("sex", "length", 'diameter', 'height', 'whole_weight', 'shucked_wieght', 'viscera_wieght', 'shell_weight' )

####################################################
Algorithm='quadratic'
solver="lpSolve"
# solver="optiSolve"
# bmethod='maxminb'
bmethod='meanb'     
C1=1.5
C2=1.5
kernelname='linear'
# kernelname='laplace'
# kernelname='rbf'
# kernelname='fixedwlaplace_auto'
# kernelname='fixedwrbf_auto'

#######################
model<-WSOSVM_f(C1,C2,TrainData,features,kernelname,Algorithm,solver,bmethod)

# prediction
Labels <- Predict_f(TestData,model)

##########################################
MAE=MSE_f(predictedLabels=Labels,Data=TestData)
MAE
loss01=loss01_f(predictedLabels=Labels,Data=TestData)
loss01


K=max(TrainData[,'yr'])
diagAcc<-diagAcc(predictedLabels=Labels,Data=TestData,K)
diagAcc
print(table(y=TestData[,'y'],Labels))
