# R code for proposed method for electricity data

# Load necessary packages
library(readxl)
library(pracma)

# Load necessary source file
source('C:/Downloads/R_functions_for_proposed_method_for_functional_responses.R')

# Load necessary dll files
dyn.load('C:/Downloads/CBS_continuous_L2.dll') # path of the CBS_continuous_L2.dll file
dyn.load('C:/Downloads/SBF_continuous_L2.dll') # path of the CBS_continuous_L2.dll file

# Get response
electricity=read_xlsx('C:/Downloads/electricity.xlsx') # path of the electricity.xlsx file
Y_raw=t(electricity[,-1]) # make response as a (sample size,number of observed  times) matrix
n=nrow(Y_raw)
T=ncol(Y_raw)

# Pre-smoothing step
time_vector=seq(0,1,length=T) # equally-spaced observed times re-scailed on [0,1]
eval_points=277 # number of evaluation time points
eval_vector=seq(0,1,length=eval_points) # vector of evaluation time points
h_add=0.1 # (minimum bandwidth,minimum bandwidth+h_add) will be the range for candidate bandwidths 
h_length=101 # number of candidate bandwidths
pre_smoothing_h=c() # its i-th component will be an optimal bandwidth for pre-smoothing the i-th curve
Y=matrix(,n,eval_points) # pre-smoothed response matrix
for(i in 1:n)
{
  pre_smoothing_h[i]=optimal_h_loocv(time_vector,Y_raw[i,],h_add,h_length)$h_optimal
  for(t in 1:eval_points)
  {
    Y[i,t]=nw(eval_vector[t],time_vector,Y_raw[i,],pre_smoothing_h[i])
  }
}

# Get predictor
X1=read_xlsx('C:/Users/u0124228/Downloads/temperature.xlsx')[,2] # path of the temperature.xlsx file
X2=read_xlsx('C:/Users/u0124228/Downloads/cloudiness.xlsx')[,2] # path of the cloudiness.xlsx file
X=cbind(X1,X2)
d=ncol(X)

# Re-scale X for it to take values in [0,1]^d
for(j in 1:d)
{
  X[,j]=(X[,j]-min(X[,j]))/(max(X[,j])-min(X[,j]))
}

# For saving results
optimal_h=matrix(,nrow=n,ncol=d)
Y_hat=array(,dim=c(n,eval_points))
error=c()

# Get ASPE
for(k in 1:n)
{
  print(k)
  X_test=X[k,]
  Y_test=Y[k,]
  X_training=X[-k,]
  Y_training=Y[-k,]
  optimal_h[k,]=CBS_L2(X_training,Y_training)$optimal_h
  Y_hat[k,]=SBF_L2(X_test,X_training,Y_training,h=optimal_h[k,])$Y_hat
  error[k]=trapz(eval_vector,(Y_test-Y_hat[k,])^2)
}

mean(error)
