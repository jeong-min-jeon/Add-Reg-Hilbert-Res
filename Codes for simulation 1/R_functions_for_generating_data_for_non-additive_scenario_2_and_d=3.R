# R functions for generating data for non-additive scenario II and d=3

# Load necessary package
library(pracma)

f123=function(x1,x2,x3,t)
{
  log((x1+x2+x3)*t/2+2)
}

f_error=function(x,t)
{
  exp(-t^4*x)
}

f=function(x1,x2,x3,e,time_vector)
{
  temp=f123(x1,x2,x3,time_vector)*f_error(e,time_vector)
  temp/trapz(time_vector,temp)
}

# Function for data generation
# M: number of monte-carlo simulation
# n: training sample size
# d: dimension of predictor
# N: test sample size
# error_range: error follows U[-error_range,error_range]
# time_vector: vector of evaluation time points
data_make=function(M,n,d,N,error_range,time_vector)
{
  set.seed(0)
  X=array(runif(M*n*d),dim=c(M,n,d))
  e=matrix(runif(M*n,min=-error_range,max=error_range),M,n)
  T=length(time_vector)
  Y=array(,dim=c(M,n,T))
  for(j in 1:M)
  {
    for(i in 1:n)
    {
      Y[j,i,]=f(X[j,i,1],X[j,i,2],X[j,i,3],e[j,i],time_vector)
    }
  }
  set.seed(1)
  test_X=array(runif(M*N*d),dim=c(M,N,d))
  test_e=matrix(runif(M*N,min=-error_range,max=error_range),M,N)
  test_Y=array(,dim=c(M,N,T))
  for(j in 1:M)
  {
    for(i in 1:N)
    {
      test_Y[j,i,]=f(test_X[j,i,1],test_X[j,i,2],test_X[j,i,3],test_e[j,i],time_vector)
    }
  }
  return(list(X=X,Y=Y,test_X=test_X,test_Y=test_Y))
}

data=data_make(100,100,3,100,1,seq(-0.5,0.5,length=101)) # for n=100
data=data_make(100,400,3,100,1,seq(-0.5,0.5,length=101)) # for n=400