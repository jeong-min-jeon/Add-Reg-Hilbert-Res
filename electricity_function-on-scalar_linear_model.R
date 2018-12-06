# R code for function-on-scalar additive model by Scheipl et al. (2015) for electricity data

# Load necessary packages
library(readxl)
library(refund)
library(pracma)

# Load necessary source file
source('C:/Downloads/R_functions_for_pre-smoothing.R') # path of the R_functions_for_pre-smoothing.R file

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

Y_hat=matrix(,nrow=n,ncol=eval_points)
error=c()
k=100 # number of spline basis

# Get ASPE
# It takes long time
for(i in 1:n)
{
  print(i)
  data_training=list(Y=Y[-i,],X1=X[-i,1],X2=X[-i,2])
  fit=pffr(Y~X1+X2,data=data_training,yind=eval_vector,bs.yindex=list(bs='ps',k=k,m=c(2,1)),bs.int=list(bs='ps',k=k,m=c(2,1)))
  Y_hat[i,]=predict(fit,newdata=list(X1=X[i,1],X2=X[i,2]))
  error[i]=trapz(eval_vector,(Y[i,]-Y_hat[i,])^2)
}

mean(error)
