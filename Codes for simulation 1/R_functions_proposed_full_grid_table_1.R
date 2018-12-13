# Proposed method with full grid search for table 1

# Load necessary package
library(pracma)

# Load necessary source file
source('C:/Downloads/R_functions_for_proposed_method_full_grid_for_density_responses.R')

# Load necessary dll files
dyn.load('C:/Downloads/Full_grid_continuous_density.dll')
dyn.load('C:/Downloads/SBF_continuous_density.dll')

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

# Function for results
# X: array for predictor of dimension (number of monte-carlo simulation(M) by training sample size(n) by dimension of predictor(d))
# Y: array for response of dimension (number of monte-carlo(M) by training sample size(n) by number of evaluation time points(T))
# test_X: array for predictor of dimension (number of monte-carlo simulation(M) by test sample size(N) by dimension of predictor(d))
# test_Y: array for response of dimension (number of monte-carlo(M) by test sample size(N) by number of evaluation time points(T))
# time_vector: vector of evaluation time points
# h_add and h_length: d-dimensional vectors 
# seq(min_h[j],min_h[j]+h_add[j],length=h_length[j]) will be the candidate set of j-th bandwidth
# for some small bandwidth min_h[j]
simulation_proposed_full_grid=function(X,Y,test_X,test_Y,time_vector,h_length,h_add)
{
  M=dim(X)[1]
  n=dim(X)[2]
  d=dim(X)[3]
  N=dim(test_X)[2]
  T=length(time_vector)
  
  optimal_h=matrix(,nrow=M,ncol=d)
  Y_hat=array(,dim=c(M,N,T))
  error=matrix(,nrow=M,ncol=N)
  
  t=0
  for(j in 1:M)
  {
    print(j)
    t1=Sys.time()
    result=Full_grid_density(X[j,,],Y[j,,],add=2,Y_grid=time_vector,h_length=h_length,h_add=h_add)
    t2=Sys.time()
    t=t+t2-t1
    optimal_h[j,]=result$optimal_h
    Y_hat[j,,]=SBF_density(test_X[j,,],X[j,,],Y[j,,],h=optimal_h[j,],Y_grid=time_vector)$Y_hat
    for(i in 1:N)
    {
      error[j,i]=density_norm(test_Y[j,i,],Y_hat[j,i,],time_vector)
    }
  }
  
  return(list(optimal_h=optimal_h,Y_hat=Y_hat,error=error,time=t))
}

d=dim(data$X)[3]
result=simulation_proposed_full_grid(data$X,data$Y,data$test_X,data$test_Y,seq(-0.5,0.5,length=101),h_length=rep(21,d),h_add=rep(0.2,d))
