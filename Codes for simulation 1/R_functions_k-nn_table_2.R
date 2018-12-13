# Kernel-based functional k-nn estimator for table 2

# Load necessary package
library(pracma)

# Function for calculating the norm of density space
density_norm=function(d1,d2,time_vector)
{
  T=length(time_vector)
  first_integral=c()
  for(t in 1:T)
  {
    first_integral[t]=trapz(time_vector,(log(d1/d1[t])-log(d2/d2[t]))^2)
  }
  trapz(time_vector,first_integral)/2
}

# Kernel function for the kernel-based functional k-nn
K=function(x)
{
  3/4*(1-x^2)*dunif(x,-1,1)*2
}

# Function for prediction for the kernel-based k-nn for density response
# x: target x (matrix whose column is the number of covariates)
# X: observed X (matrix of size (sample size,the number of covariates))
# Y: observed Y (matrix of size (sample size,the number of evaluation time points))
# time_vector: vector of evaluation time points
# k: the number of nearest neighbors
predict_knn_density=function(x,X,Y,time_vector,k)
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
    upper=colProds(Y^K_values[p,])
    temporary=upper^{1/below}
    predict[p,]=temporary/trapz(time_vector,temporary)
  }
  return(predict)
}

# Function for cross-validatory optimal k for the kernel-based k-nn for density response
# X: observed X (matrix of size (sample size,the number of covariates))
# Y: observed Y (matrix of size (sample size,the number of evaluation time points))
# time_vector: vector of evaluation time points
# nfolds: the number of folds for cross-validation
# k_vec: vector of candidate k
optimal_k_density=function(X,Y,time_vector,nfolds,k_vec)
{
  n=nrow(X)
  d=ncol(X)
  s=sample(n)
  X=matrix(X[s,],ncol=d)
  Y=Y[s,]
  folds=cut(1:n,breaks=nfolds,labels=FALSE)
  error.k=c()
  for(j in 1:length(k_vec))
  {
    error.fold=rep(0,nfolds)
    for(k in 1:nfolds)
    {
      X.training=X[-which(folds==k),]
      X.test=X[which(folds==k),]
      Y.training=Y[-which(folds==k),]
      Y.test=Y[which(folds==k),]
      Y.test.hat=predict_knn_density(matrix(X.test,ncol=d),X.training,Y.training,time_vector,k_vec[j])
      for(i in 1:length(which(folds==k)))
      {
        error.fold[k]=error.fold[k]+density_norm(Y.test[i,],Y.test.hat[i,],time_vector)
      }
      if(length(which(folds==k))==1) error.fold[k]=density_norm(Y.test,Y.test.hat,time_vector)
    }
    error.k[j]=sum(error.fold)
  }
  return(k_vec[which(error.k==min(error.k))])
}

# Function for results
# X: array for predictor of dimension (number of monte-carlo simulation(M) by training sample size(n) by dimension of predictor(d))
# Y: array for response of dimension (number of monte-carlo(M) by training sample size(n) by number of evaluation time points(T))
# test_X: array for predictor of dimension (number of monte-carlo simulation(M) by test sample size(N) by dimension of predictor(d))
# test_Y: array for response of dimension (number of monte-carlo(M) by test sample size(N) by number of evaluation time points(T))
# time_vector: vector of evaluation time points
# k_vec: vector of candidiate k
simulation_knn=function(X,Y,test_X,test_Y,time_vector,k_vec)
{
  M=dim(X)[1]
  n=dim(X)[2]
  d=dim(X)[3]
  N=dim(test_X)[2]
  T=length(time_vector)
  
  k_selected=c()
  Y_hat=array(,dim=c(M,N,T))
  error=matrix(,nrow=M,ncol=N)
  
  for(j in 1:M)
  {
    print(j)
    k_selected[j]=optimal_k_density(X[j,,],Y[j,,],time_vector,nfolds=10,k_vec)
    Y_hat[j,,]=predict_knn_density(test_X[j,,],X[j,,],Y[j,,],time_vector,k_selected[j])
    for(i in 1:N)
    {
      error[j,i]=density_norm(test_Y[j,i,],Y_hat[j,i,],time_vector)
    }
  }
  
  return(list(k_selected=k_selected,Y_hat=Y_hat,error=error))
}

result=simulation_knn(data$X,data$Y,data$test_X,data$test_Y,seq(-0.5,0.5,length=101),k_vec=1:30)
