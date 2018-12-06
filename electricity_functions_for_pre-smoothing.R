# R functions for pre-smoothing

# Kernel function for pre-smoothing
K=function(x)
{
  dnorm(x,sd=0.5)
}

# Nadaraya-Watson estimator for pre-smoothing
# t: target time
# T: vector of observed times
# Y: curve observed on T
# h: bandwidth
nw=function(t,T,Y,h)
{
  K_values=K(abs(t-T)/h)
  below=sum(K_values)
  upper=sum(K_values*Y)
  predict=upper/below
  return(predict)
}

# Leave-one-out cross-validatory bandwidth selection for pre-smoothing
# T: vector of equally-spaced observed times re-scailed on [0,1] (min of observed times=0,max of observed times=1)
# Y: curve observed on T
# h_add and h_length: seq(min_h,min_h+h_add,length=h_length) will be the set of candidate bandwidths
# for some small bandwidth min_h which makes kernel smoothing possible
optimal_h_loocv=function(T,Y,h_add,h_length)
{
  n=length(T)
  min_h=1/n+0.001
  h_vector=seq(min_h,min_h+h_add,length=h_length)
  error=matrix(,nrow=h_length,ncol=n)
  error_h=c()
  for(j in 1:h_length)
  {
    for(i in 1:n)
    {
      distance=min(abs(T[i]-T[-i]))
      if(distance<h_vector[j])
      {
        Y_hat=nw(T[i],T[-i],Y[-i],h_vector[j])
        error[j,i]=(Y[i]-Y_hat)^2
      }
    }
    error_h[j]=mean(error[j,is.na(error[j,])==FALSE])
  }
  return(list(h_optimal=h_vector[which.min(error_h)]))
}
