# R code for kernel-based function-on-scalar k-nn method for electricity data

# Load necessary packages
library(readxl)
library(pracma)
library(pdist)

# Get response
electricity=read_xlsx('C:/Downloads/electricity.xlsx') # path of the electricity.xlsx file
Y_raw=t(electricity[,-1]) # make response as a (sample size,number of observed  times) matrix
n=nrow(Y_raw)
T=ncol(Y_raw)

# Pre-smoothing step
time_vector=seq(0,1,length=T) # equally-spaced observed times re-scailed on [0,1]
eval_points=277 # number of evaluation time points
eval_vector=seq(0,1,length=eval_points) # vector of evaluation time points
h_add=0.1 # (minimum bandwidth,minimum bandwidth+h_add) will be the range for candidate bandwidths 
h_length=101 # number of candidate bandwidths
pre_smoothing_h=c() # its i-th component will be an optimal bandwidth for pre-smoothing the i-th curve
Y=matrix(,n,eval_points) # pre-smoothed response matrix
for(i in 1:n)
{
  pre_smoothing_h[i]=optimal_h_loocv(time_vector,Y_raw[i,],h_add,h_length)$h_optimal
  for(t in 1:eval_points)
  {
    Y[i,t]=nw(eval_vector[t],time_vector,Y_raw[i,],pre_smoothing_h[i])
  }
}

# Get predictor
X1=read_xlsx('C:/Users/u0124228/Downloads/temperature.xlsx')[,2] # path of the temperature.xlsx file
X2=read_xlsx('C:/Users/u0124228/Downloads/cloudiness.xlsx')[,2] # path of the cloudiness.xlsx file
X=cbind(X1,X2)
d=ncol(X)

# Re-scale X for better result
for(j in 1:d)
{
  X[,j]=(X[,j]-min(X[,j]))/(max(X[,j])-min(X[,j]))
}

# Kernel function for the kernel-based functional k-nn
K=function(x)
{
  3/4*(1-x^2)*dunif(x,-1,1)*2
}

# Vectorized trapzoidal integration
int<-function(x,y)
{
  index = 2:length(x)
  ((x[index] - x[index-1]) %*% (y[index,] + y[index-1,])) / 2
}

# Function for prediction for the kernel-based functional k-nn
# x: target x (matrix whose column is the number of covariates)
# X: observed X (matrix of size (sample size,the number of covariates))
# Y: observed Y (matrix of size (sample size,the number of evaluation time points))
# k: the number of nearest neighbors
predict_knn=function(x,X,Y,k)
{
  N=nrow(x)
  n=nrow(X)
  T=ncol(Y)
  predict=matrix(,nrow=N,ncol=T)
  K_values=matrix(,nrow=N,ncol=n)
  h=c()
  for(p in 1:N)
  {
    h[p]=sort(as.matrix(pdist(x[p,],X)))[k]+0.001
    K_values[p,]=K(as.matrix(pdist(x[p,],X))/h[p])
    below=sum(K_values[p,])
    upper=colSums(K_values[p,]*Y)
    predict[p,]=upper/below
  }
  return(predict)
}

# Function for cross-validatory optimal k for the kernel-based functional k-nn
# X: observed X (matrix of size (sample size,the number of covariates))
# Y: observed Y (matrix of size (sample size,the number of evaluation time points))
# time_vector: time vector for Y
# nfolds: the number of folds for cross-validation
# k_vec: vector of candidate k
optimal_k=function(X,Y,time_vector,nfolds,k_vec)
{
  n=nrow(X)
  d=ncol(X)
  s=sample(n)
  X=matrix(X[s,],ncol=d)
  Y=Y[s,]
  folds=cut(1:n,breaks=nfolds,labels=FALSE)
  error.fold=c()
  error.h=c()
  for(j in 1:length(k_vec))
  {
    for(k in 1:nfolds)
    {
      X.training=X[-which(folds==k),]
      X.test=X[which(folds==k),]
      Y.training=Y[-which(folds==k),]
      Y.test=Y[which(folds==k),]
      Y.test.hat=predict_knn(matrix(X.test,ncol=d),X.training,Y.training,k_vec[j])
      if(length(which(folds==k))>1) error.fold[k]=sum(int(time_vector,(t(Y.test)-t(Y.test.hat))^2))
      if(length(which(folds==k))==1) error.fold[k]=trapz(time_vector,(Y.test-Y.test.hat)^2)
    }
    error.h[j]=sum(error.fold)
  }
  return(k_vec[which(error.h==min(error.h))])
}

Y_hat=matrix(,n,eval_points)
error=c()
k_selected=c()

# Get ASPE
for(i in 1:n)
{
  print(i)
  k_selected[i]=optimal_k(X[-i,],Y[-i,],eval_vector,nfolds=10,k_vec=c(1:30))
  Y_hat[i,]=predict_knn(matrix(X[i,],1,d),X[-i,],Y[-i,],k_selected[i])
  error[i]=trapz(eval_vector,(Y[i,]-Y_hat[i,])^2)
}

mean(error)
k_selected