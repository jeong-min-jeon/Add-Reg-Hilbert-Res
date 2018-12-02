# Kernel function (will be used in pre-smoothing)
K=function(x)
{
  dnorm(x,sd=0.5)
}

# Nadaraya-Watson estimator (will be used in pre-smoothing)
# x: target x
# X: observed X
# Y: observed Y
# h: bandwidth
predict_npfda_loocv_univariate=function(x,X,Y,h)
{
  K_values=K(abs(x-X)/h)
  below=sum(K_values)
  upper=sum(K_values*Y)
  predict=upper/below
  return(predict)
}

# Leave-one-out cross-validatory bandwidth selection used in pre-smoothing
# T: a vector of equally-spaced observed times re-scailed on [0,1] 
# Y: an observed curve of length length(T)
# h_add and h_length: seq(min_h,min_h+h_add,length=h_length) will be the set of candidate bandwidths
# for some minimum bandwidth min_h
optimal_h_loocv_univariate=function(T,Y,h_add,h_length)
{
  n=length(T)
  min_h=1/n+0.001 # minimum bandwidth
  h_vector=seq(min_h,min_h+h_add,length=h_length) # the set of candidate bandwidths
  error=matrix(,nrow=h_length,ncol=n)
  error_h=c()
  for(j in 1:h_length)
  {
    for(i in 1:n)
    {
      distance=min(abs(T[i]-T[-i]))
      if(distance<h_vector[j])
      {
        yhat=predict_npfda_loocv_univariate(T[i],T[-i],Y[-i],h_vector[j])
        error[j,i]=(Y[i]-yhat)^2
      }
    }
    error_h[j]=mean(error[j,is.na(error[j,])==FALSE])
  }
  return(list(h_optimal=h_vector[which.min(error_h)],h_initial=min_h))
}