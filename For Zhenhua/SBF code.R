dyn.load('C:/Users/u0124228/Desktop/SBF_continuous_pointwise.dll')

# Function for obtaining additional grid points for numerical integration
# X_test: N times d matrix whose elements take values in [0,1]
# We use sort(X_test[,j],seq(0,1,length=add)) as the grid points for numerical integration for jth integral
integration.grid=function(X_test,add)
{
  add_vector=seq(0,1,length=add)
  added_X_test=matrix(,add+nrow(X_test),ncol(X_test))
  for(j in 1:ncol(X_test))
  {
    added_X_test[,j]=sort(c(add_vector,X_test[,j]))
  }
  added_X_test
}

# Function for fitting smooth backfitting
# X_test: N times d matrix whose elements take values in [0,1]
# X_training: n times d matrix whose elements take values in [0,1]
# Y_training: n-dimensional vector
# add: number of additional grid points for numerical integration
# max_iteration: maximum iteration number for smooth bacfkfitting algorithm
# epsilon: convergence criterion for smooth backfitting algorithm
# h: d-dimensional vector of bandwidths
SBF=function(X_test,X_training,Y_training,add=11,max_iteration=50,epsilon=10^-4,h)
{
  n=nrow(X_training)
  d=ncol(X_training)
  N=nrow(X_test)
  
  add=as.integer(add)
  max_iteration=as.integer(max_iteration)
  h=as.double(h)
  
  g=integration.grid(X_test,add)
  ngrid=nrow(g)
  actual_iteration=as.integer(0)
  T=as.integer(1)
  Y_training=matrix(Y_training,n,T)
  yhat=matrix(0,N,T)
  mhat=array(0,dim=c(d,ngrid,T))
  
  result=.Fortran('SBF_continuous_pointwise',n0=N,n=n,d=d,ngrid=ngrid,T=T,X0=X_test,X=X_training,Y=Y_training,g=g,h=h,
                  max_iteration=max_iteration,epsilon=epsilon,actual_iteration=actual_iteration,
                  yhat=yhat,mhat=mhat)
  
  return(list(yhat=result$yhat,mhat=result$mhat,actual_iteration=result$actual_iteration))
  if(any(result$actual_iteration==max_iteration)) print('Some smooth backfitting algorithm did not converge.')
}

# Example
n=200
d=2
X_training=matrix(runif(n*d),nrow=n,ncol=d)
Y_training=sin(X_training[,1]*pi)+cos(X_training[,2]*pi)+rnorm(n,sd=0.1)
h=rep(0.1,d)
X_test=matrix(runif(10*d),nrow=10,ncol=d)
result=SBF(X_test=X_test,X_training=X_training,Y_training=Y_training,h=h)
result$actual_iteration # iteration number r
as.vector(result$yhat) # predicted Y_test
sin(X_test[,1]*pi)+cos(X_test[,2]*pi) # true regression values at X_test
