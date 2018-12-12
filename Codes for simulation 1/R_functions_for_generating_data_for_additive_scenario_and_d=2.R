# R functions for generating data for additive scenario and d=2

# Load necessary package
library(pracma)

# Defining component maps
f1=function(x,t)
{
  exp(-x*t)
}

f2=function(x,t)
{
  exp(-2*x^2*t^2)
}

# Defining error map
f_error=function(x,t)
{
  exp(-x*t^4)
}

# Defining combined map
f=function(x1,x2,e,time_vector)
{
  temp=f1(x1,time_vector)*f2(x2,time_vector)*f_error(e,time_vector)
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
      Y[j,i,]=f(X[j,i,1],X[j,i,2],e[j,i],time_vector)
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
      test_Y[j,i,]=f(test_X[j,i,1],test_X[j,i,2],test_e[j,i],time_vector)
    }
  }
  return(list(X=X,Y=Y,test_X=test_X,test_Y=test_Y))
}

data=data_make(100,100,2,100,1,seq(-0.5,0.5,length=101)) # for n=100
data=data_make(100,400,2,100,1,seq(-0.5,0.5,length=101)) # for n=400