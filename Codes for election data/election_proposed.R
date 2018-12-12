# R code for proposed method for election data

# Load necessary R package
library(Compositional)

# Load necessary source file
source('C:/Downloads/R_functions_for_proposed_method_for_simplex_responses.R') # path of the file

# Load necessary dll files
dyn.load('C:/Downloads/CBS_continuous_simplex.dll') # path of the CBS_continuous_simplex.dll file
dyn.load('C:/Downloads/SBF_continuous_simplex.dll') # path of the CBS_continuous_simplex.dll file

# Get data
raw_data <- read.csv("https://raw.githubusercontent.com/OhmyNews/2017-Election/master/170509_Elec_Data-utf8.csv")
selected_data=raw_data[,9:13]
colnames(selected_data)[1]<-'avg_age'
colnames(selected_data)[2]<-'avg_years_of_education'
colnames(selected_data)[3]<-'avg_housing_price'
colnames(selected_data)[4]<-'avg_health_insurance_price'
n=250
d=4
T=5
# Define X
X=as.matrix(selected_data[seq(1,length=n,by=T),1:d])
for(j in 1:d)
{
  X[,j]=(X[,j]-min(X[,j]))/(max(X[,j])-min(X[,j]))
}
# Define Y
Y=matrix(selected_data[,5],nrow=n,ncol=T,byrow='TRUE')
Y=Y[,1:3]
T=3
for(i in 1:n)
{
  Y[i,]=Y[i,]/sum(Y[i,])
}
# Suffle data
set.seed(0)
shuffle=sample(1:n)
X=X[shuffle,]
Y=Y[shuffle,]

# Function for computing square norm of y1 minus y2
# y1 and y2: compositional vectors
comp_distance=function(y1,y2)
{
  D=length(y1)
  dist=matrix(,D,D)
  for(i in 1:D)
  {
    for(j in 1:D)
    {
      dist[i,j]=(log(y1[i]/y1[j])-log(y2[i]/y2[j]))^2
    }
  }
  sum(dist)/2/D
}

# Index of 10 partitions of sample
nfolds=10
folds <- cut(1:n,breaks=nfolds,labels=FALSE)

# For saving results
optimal_h=matrix(,nrow=nfolds,ncol=d)
Y_hat=array(,dim=c(nfolds,25,T))
error=c()

# Get ASPE
for(k in 1:nfolds)
{
  print(k)
  s=which(folds==k)
  X_test=X[s,]
  Y_test=Y[s,]
  X_training=X[-s,]
  Y_training=Y[-s,]
  optimal_h[k,]=CBS_simplex(X_training,Y_training,h_length=rep(41,d),h_add=rep(0.2,d))$optimal_h # Reduced bw grid for fast computation
  Y_hat[k,,]=SBF_simplex(X_test,X_training,Y_training,h=optimal_h[k,])$Y_hat
  temp_dist=c()
  for(i in 1:nrow(Y_test))
  {
    temp_dist[i]=comp_distance(as.vector(Y_test[i,]),as.vector(Y_hat[k,i,]))
  }
  error[k]=sum(temp_dist)/nrow(Y_test)
}

mean(error)
