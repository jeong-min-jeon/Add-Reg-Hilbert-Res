# R functions for m(3,4) propsed method for table 2

# Load necessary packages
library(truncdist)
library(pracma)
library(EnvStats)
library(ald)

# Load necessary source file
source('C:/Downloads/R_functions_for_proposed_method_for_density_responses.R')

# Load necessary dll files
dyn.load('C:/Downloads/CBS_continuous_density.dll')
dyn.load('C:/Downloads/SBF_continuous_density.dll')

# Defining true component maps for calculating IMSE, ISB, IV
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

# Function for calculating the norm of density space
density_norm=function(d1,d2,time_vector)
{
  T=length(time_vector)
  first_integral=c()
  for(t in 1:T)
  {
    first_integral[t]=trapz(time_vector,(log(d1[time_vector]/d1[t])-log(d2[time_vector]/d2[t]))^2)
  }
  trapz(time_vector,first_integral)
}

# Function for obtaining IMSE, ISB, IV
# X: array for predictor of dimension (number of monte-carlo simulation(M) by sample size(n) by dimension of predictor(d))
# Y: array for response of dimension (number of monte-carlo(M) by sample size(n) by number of evaluation time points(T))
# time_vector: vector of evaluation time points
# h_add and h_length: d-dimensional vectors 
# seq(min_h[j],min_h[j]+h_add[j],length=h_length[j]) will be the candidate set of j-th bandwidth
# for some small bandwidth min_h[j]
simulation_proposed=function(X,Y,time_vector,h_length,h_add)
{
  M=dim(X)[1]
  n=dim(X)[2]
  d=2
  T=length(time_vector)
  
  # Defining the true first and second component maps
  x_add=101
  grid=seq(0,1,length=x_add)
  grid_matrix=matrix(rep(grid,d),ncol=d)
  m_true=array(,dim=c(d,x_add,T))
  for(k in 1:d)
  {
    for(t in 1:T)
    {
      if(k==1)  m_true[k,,t]=f3(grid,time_vector[t],time_vector)
      if(k==2)  m_true[k,,t]=f4(grid,time_vector[t],time_vector)
    }
  }
  
  # Obtaining oracle information for the 3rd and 4th component maps
  m_true_data=array(,dim=c(d,M,n,T))
  for(k in 1:2)
  {
    for(j in 1:M)
    {
      for(t in 1:T)
      {
        if(k==1)  m_true_data[k,j,,t]=f1(X[j,,1],time_vector[t],time_vector)
        if(k==2)  m_true_data[k,j,,t]=f2(X[j,,2],time_vector[t],time_vector)
      }
    }
  }
  
  for(j in 1:M)
  {
    for(i in 1:n)
    {
      temporary=Y[j,i,]/m_true_data[1,j,i,]/m_true_data[2,j,i,]
      Y[j,i,]=temporary/trapz(time_vector,temporary)
    }
  }
  
  m_hat=array(,dim=c(M,d,x_add,T))
  
  # Estimation
  for(j in 1:M)
  {
    print(j)
    optimal_h=CBS_density(X[j,,3:4],Y[j,,],add=2,Y_grid=time_vector,h_length=h_length,h_add=h_add)$optimal_h
    m_hat[j,,,]=SBF_density(grid_matrix,X[j,,3:4],Y[j,,],h=optimal_h,add=0,Y_grid=time_vector)$m_hat
  }
  
  # Compute average of m_hat
  avg_m_hat=array(,dim=c(d,x_add,T))
  for(k in 1:d)
  {
    for(l in 1:x_add)
    {
      prods=colProds(m_hat[,k,l,]^{1/M})
      avg_m_hat[k,l,]=prods/trapz(time_vector,prods)
    }
  }
  
  bias=c()
  variance=c()
  norm_vector=c()
  norm_temporary=c()
  
  # ISB
  for(k in 1:d)
  {
    for(l in 1:x_add)
    {
      norm_vector[l]=density_norm(m_true[k,l,],avg_m_hat[k,l,],time_vector)
    }
    bias[k]=trapz(grid,norm_vector)
  }
  
  # IV
  for(k in 1:d)
  {
    for(l in 1:x_add)
    {
      for(j in 1:M)
      {
        norm_temporary[j]=density_norm(m_hat[j,k,l,],avg_m_hat[k,l,],time_vector)
      }
      norm_vector[l]=mean(norm_temporary)
    }
    variance[k]=trapz(grid,norm_vector)
  }
  
  return(list(bias=bias,variance=variance))
}

result=simulation_proposed(data$X,data$Y,seq(-0.5,0.5,length=101),h_length=c(21,21),h_add=c(0.2,0.2))

result$bias
result$variance
result$bias+result$variance
