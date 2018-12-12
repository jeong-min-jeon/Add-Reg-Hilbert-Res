# R functions for proposed method for simplex response

# Load necessary packages
library(pdist)
library(matrixStats)
library(fields)

# Function for obtaining candidate bandwidths
# X: observed X (matrix of size (sample size,number of covariates))
# add: number of additional grid points for numerical integration
# cv: number of folds for cross-validation
# h_add and h_length: ncol(X)-dimensional vectors
# seq(min_h[j],min_h[j]+h_add[j],length=h_length[j]) will be the candidate set of j-th bandwidth
# for some small bandwidth min_h[j]
# resulting h_candidate_matrix is the matrix whose j-th column is the candidate set of j-th bandwidth
get_h_candidate=function(X,add,cv,h_add,h_length)
{
  n=nrow(X)
  d=ncol(X)
  grid_add=seq(0,1,length=add)
  s=sample(1:n)
  X=X[s,]
  folds=new_cut(n,cv)
  distance=matrix(,nrow=cv,ncol=d)
  for(k in 1:cv)
  {
    X_training=X[-which(folds==k),]
    X_test=X[which(folds==k),]
    for(j in 1:d)
    {
      if(is.vector(X_test)) grid_vector=c(X_test[j],grid_add)
      else grid_vector=c(X_test[,j],grid_add)
      distance[k,j]=max(rowMins(rdist(grid_vector,X_training[,j])))
    }
  }
  min_h=c()
  h_candidate_matrix=matrix(,nrow=max(h_length),ncol=d)
  for(j in 1:d)
  {
    min_h[j]=max(distance[,j])+0.001
    h_candidate_matrix[1:h_length[j],j]=seq(min_h[j],min_h[j]+h_add[j],length=h_length[j])
  }
  return(list(h_candidate_matrix=h_candidate_matrix,s=s))
}

# Function giving index of splited data for cross-validation
# n: sample size
# cv: number of folds for cross-validation
# only the last fold has a different number of observations
new_cut=function(n,cv)
{
  folds=c()
  size=floor(n/cv)
  size_last=n-(cv-1)*size
  for(k in 1:cv)
  {
    if(k<cv) folds[((k-1)*size+1):(k*size)]=k
    if(k==cv) folds[((cv-1)*size+1):n]=k
  }
  return(folds)
}

# Function for obtaining additional grid points for numerical integration
# x: X_test matrix taking values in [0,1]^d
# we use sort(X_test[,j],seq(0,1,length=add)) as the grid points for numerical integration for j-th integral
integration.grid=function(x,add)
{
  add_vector=seq(0,1,length=add)
  added_x=matrix(,add+nrow(x),ncol(x))
  for(j in 1:ncol(x))
  {
    added_x[,j]=sort(c(add_vector,x[,j]))
  }
  added_x
}

# Function for CBS
# X_training and Y_training are usual
# add: number of additional grid points for numerical integration
# max_iteration: maximum iteration number for Bochner smooth bacfkfitting algorithm
# epsilon: convergence criterion for Bochner smooth backfitting algorithm
# max_cbs_iteration: maximum iteration number for CBS algorithm
# cv: number of folds for cross-validation
# h_add and h_length: ncol(X_training)-dimensional vectors
# seq(min_h[j],min_h[j]+h_add[j],length=h_length[j]) will be the candidate set of j-th bandwidth
# for some small bandwidth min_h[j]
CBS_simplex=function(X_training,Y_training,add=11,max_iteration=50,epsilon=10^-4
                ,max_cbs_iteration=20,cv=10,h_add=NULL,h_length=NULL)
{
  n=nrow(X_training)
  d=ncol(X_training)
  T=ncol(Y_training)
  
  if(is.null(h_add)) h_add=rep(0.5,d)
  if(is.null(h_length)) h_length=as.integer(rep(101,d))
  else h_length=as.integer(h_length)
  add=as.integer(add)
  max_iteration=as.integer(max_iteration)
  max_cbs_iteration=as.integer(max_cbs_iteration)
  cv=as.integer(cv)
  cv_test_size=as.integer(floor(n/cv))
  
  result=get_h_candidate(X_training,add,cv,h_add,h_length)
  s=result$s
  X_training=X_training[s,]
  Y_training=Y_training[s,]
  h_candidate_matrix=result$h_candidate_matrix
  h_candidate_matrix[is.na(h_candidate_matrix)]=0
  
  result=.Fortran('CBS_continuous_simplex',n=n,d=d,T=T,X_training=X_training,Y_training=Y_training,
                  h_vector=h_candidate_matrix,add=add,max_iteration=max_iteration,epsilon=epsilon,
                  ngrid_h=max(h_length),cv=cv,cv_test_size=cv_test_size,true_cd=as.integer(0),h_length=h_length,
                  max_cd_iteration=max_cbs_iteration,optimal_h=h_candidate_matrix[1,],
                  total_iteration=array(as.integer(0),c(max_cbs_iteration,d,max(h_length),cv)))
  return(list(optimal_h=result$optimal_h,CBS_iteration=result$true_cd,SBF_iteration=result$total_iteration))
}

# Function for fitting Bochner smooth backfitting
# X_test, X_training, Y_training are usual
# add: number of additional grid points for numerical integration
# max_iteration: maximum iteration number for Bochner smooth bacfkfitting algorithm
# epsilon: convergence criterion for Bochner smooth backfitting algorithm
# h: ncol(X_training)-dimensional vector of bandwidths
SBF_simplex=function(X_test,X_training,Y_training,add=11,max_iteration=50,epsilon=10^-4,h)
{
  n=nrow(X_training)
  d=ncol(X_training)
  T=ncol(Y_training)
  if(is.matrix(X_test)) N=nrow(X_test)
  else N=as.integer(1)
  X_test=matrix(X_test,nrow=N,ncol=d)
  
  add=as.integer(add)
  max_iteration=as.integer(max_iteration)
  h=as.double(h)
  
  g=integration.grid(X_test,add)
  ngrid=nrow(g)
  actual_iteration=as.integer(0)
  yhat=matrix(0,N,T)
  mhat=array(0,dim=c(d,ngrid,T))
  
  result=.Fortran('SBF_continuous_simplex',n0=N,n=n,d=d,ngrid=ngrid,T=T,X0=X_test,X=X_training,Y=Y_training,g=g,h=h,
                  max_iteration=max_iteration,epsilon=epsilon,actual_iteration=actual_iteration,
                  yhat=yhat,mhat=mhat)
  
  return(list(Y_hat=result$yhat,m_hat=result$mhat,SBF_iteration=result$actual_iteration))
  if(any(result$actual_iteration==max_iteration)) print('Some smooth backfitting algorithm did not converge.')
}
