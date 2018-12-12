# R functions for generating data for Table 2

# Load necessary packages
library(truncdist)
library(pracma)
library(EnvStats)
library(ald)

# Defining component maps
f0=function(t,time_vector)
{
  T=length(time_vector)
  dtrunc(t,spec='cauchy',a=time_vector[1],b=time_vector[T],scale=0.2)
}

f1=function(x,t,time_vector)
{
  T=length(time_vector)
  (dtrunc(t,spec='norm',a=time_vector[1],b=time_vector[T],sd=0.5))^(cos(2*pi*x))
}

f2=function(x,t,time_vector)
{
  T=length(time_vector)
  (dtrunc(t,spec='t',a=time_vector[1],b=time_vector[T],df=0.25))^(sin(2*pi*x))
}

f3=function(x,t,time_vector)
{
  T=length(time_vector)
  normalize=pALD(time_vector[T])-pALD(time_vector[1])
  (dALD(t)/normalize)^(cos(pi*x))
}

f4=function(x,t,time_vector)
{
  T=length(time_vector)
  normalize=pnormMix(time_vector[T], mean1 = 0.3, sd1 = 0.2, mean2 = -0.3, sd2 = 0.2, p.mix = 0.5)-pnormMix(time_vector[1], mean1 = 0.3, sd1 = 0.2, mean2 = -0.3, sd2 = 0.2, p.mix = 0.5)
  (dnormMix(t, mean1 = 0.3, sd1 = 0.2, mean2 = -0.3, sd2 = 0.2, p.mix = 0.5)/normalize)^(2*x-1)
}

# Defining error map
f_error=function(x,t,time_vector)
{
  T=length(time_vector)
  (dtrunc(t,spec='logis',a=time_vector[1],b=time_vector[T]))^x
}

# Defining combined map
f=function(x1,x2,x3,x4,e,time_vector)
{
  temp=f0(time_vector,time_vector)*f1(x1,time_vector,time_vector)*f2(x2,time_vector,time_vector)*f3(x3,time_vector,time_vector)*f4(x4,time_vector,time_vector)*f_error(e,time_vector,time_vector)
  temp/trapz(time_vector,temp)
}

# Function for data generation
# M: number of monte-carlo simulation
# n: sample size
# d: dimension of predictor
# error_range: error follows U[-error_range,error_range]
# time_vector: vector of evaluation time points
data_make=function(M,n,d,error_range,time_vector)
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
      Y[j,i,]=f(X[j,i,1],X[j,i,2],X[j,i,3],X[j,i,4],e[j,i],time_vector)
    }
  }
  return(list(X=X,Y=Y))
}

# Generating data
data=data_make(500,100,4,1,seq(-0.5,0.5,length=101)) # use for n=100
data=data_make(500,400,4,1,seq(-0.5,0.5,length=101)) # use for n=400
