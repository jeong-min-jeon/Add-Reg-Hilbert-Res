# Functional n-w for table 2

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

# Kernel function for the nw
K=function(x)
{
  3/4*(1-x^2)*dunif(x,-1,1)*2
}

# Function for prediction for the n-w for density response
# x: target x (matrix whose column is the number of covariates)
# X: observed X (matrix of size (sample size,the number of covariates))
# Y: observed Y (matrix of size (sample size,the number of evaluation time points))
# time_vector: vector of evaluation time points
# h: bandwidth
predict_nw_density=function(x,X,Y,time_vector,h)
{
  N=nrow(x)
  n=nrow(X)
  T=ncol(Y)
  predict=matrix(0,nrow=N,ncol=T)
  K_values=matrix(0,nrow=N,ncol=n)
  for(p in 1:N)
  {
    K_values[p,]=K(as.matrix(pdist(x[p,],X))/h)
    below=sum(K_values[p,])
    upper=colProds(Y^K_values[p,])
    temporary=upper^{1/below}
    predict[p,]=temporary/trapz(time_vector,temporary)
  }
  return(predict)
}

# Function for cross-validatory optimal h for the nw for density response
# X: observed X (matrix of size (sample size,the number of covariates))
# Y: observed Y (matrix of size (sample size,the number of evaluation time points))
# time_vector: vector of evaluation time points
# nfolds: the number of folds for cross-validation
# h_add and h_length: seq(min_h,min_h+h_add,length=h_length) will be the set of candidate bandwidths
# for some small bandwidth min_h which makes kernel smoothing possible
optimal_h_density=function(X,Y,time_vector,nfolds,h_add,h_length)
{
  n=nrow(X)
  d=ncol(X)
  s=sample(n)
  X=X[s,]
  Y=Y[s,]
  folds=cut(1:n,breaks=nfolds,labels=FALSE)
  distance=c()
  for(k in 1:nfolds)
  {
    X.training=X[-which(folds==k),]
    X.test=X[which(folds==k),]
    new.distance=c()
    for(p in 1:nrow(X.test))
    {
      new.distance[p]=min(as.matrix(pdist(matrix(X.test[p,],ncol=d),X.training)))
    }
    distance[k]=max(new.distance)
  }
  min_h=max(distance)+0.001
  h_vector=seq(min_h,min_h+h_add,length=h_length)
  error.h=c()
  for(j in 1:h_length)
  {
    error.fold=rep(0,nfolds)
    for(k in 1:nfolds)
    {
      X.training=X[-which(folds==k),]
      X.test=X[which(folds==k),]
      Y.training=Y[-which(folds==k),]
      Y.test=Y[which(folds==k),]
      Y.test.hat=predict_nw_density(matrix(X.test,ncol=d),X.training,Y.training,time_vector,h_vector[j])
      for(i in 1:length(which(folds==k)))
      {
        error.fold[k]=error.fold[k]+density_norm(Y.test[i,],Y.test.hat[i,],time_vector)
      }
    }
    error.h[j]=sum(error.fold)
  }
  return(list(h_optimal=min(h_vector[which.min(error.h)]),h_initial=min_h))
}

# Function for results
# X: array for predictor of dimension (number of monte-carlo simulation(M) by training sample size(n) by dimension of predictor(d))
# Y: array for response of dimension (number of monte-carlo(M) by training sample size(n) by number of evaluation time points(T))
# test_X: array for predictor of dimension (number of monte-carlo simulation(M) by test sample size(N) by dimension of predictor(d))
# test_Y: array for response of dimension (number of monte-carlo(M) by test sample size(N) by number of evaluation time points(T))
# time_vector: vector of evaluation time points
# h_add and h_length: seq(min_h,min_h+h_add,length=h_length) will be the set of candidate bandwidths
# for some small bandwidth min_h which makes kernel smoothing possible
simulation_nw=function(X,Y,test_X,test_Y,time_vector,h_add,h_length)
{
  M=dim(X)[1]
  n=dim(X)[2]
  d=dim(X)[3]
  N=dim(test_X)[2]
  T=length(time_vector)
  
  h_selected=c()
  Y_hat=array(,dim=c(M,N,T))
  error=matrix(,nrow=M,ncol=N)
  
  for(j in 1:M)
  {
    print(j)
    h_selected[j]=optimal_h_density(X[j,,],Y[j,,],time_vector,nfolds=10,h_add,h_length)$h_optimal
    Y_hat[j,,]=predict_nw_density(test_X[j,,],X[j,,],Y[j,,],time_vector,h_selected[j])
    for(i in 1:N)
    {
      error[j,i]=density_norm(test_Y[j,i,],Y_hat[j,i,],time_vector)
    }
  }
  
  return(list(h_selected=h_selected,Y_hat=Y_hat,error=error))
}

result=simulation_nw(data$X,data$Y,data$test_X,data$test_Y,seq(-0.5,0.5,length=101),h_add=0.2,h_length=201)
