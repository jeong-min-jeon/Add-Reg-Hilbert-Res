library(Compositional)

# Get data
raw_data <- read.csv("https://raw.githubusercontent.com/OhmyNews/2017-Election/master/170509_Elec_Data-utf8.csv")
selected_data=raw_data[,9:13]
colnames(selected_data)[1]<-'avg_age'
colnames(selected_data)[2]<-'avg_educational_years'
colnames(selected_data)[3]<-'avg_housing_price'
colnames(selected_data)[4]<-'avg_health_insurance_price'
n=250
d=4
T=5
# Define X
X=as.matrix(selected_data[seq(1,length=n,by=T),1:d])
# Define Y
Y=matrix(selected_data[,5],nrow=n,ncol=T,byrow='TRUE')
Y=Y[,1:3]
T=3
for(i in 1:n)
{
  Y[i,]=Y[i,]/sum(Y[i,])
}
# Suffle the sample
set.seed(0)
shuffle=sample(1:n)
X=X[shuffle,]
Y=Y[shuffle,]

# Define norm square for the simplex space
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

# Index number for 10 partitions of sample
nfolds=10
folds <- cut(1:n,breaks=nfolds,labels=FALSE)

# Dirichlet regression
error_diri=c()
for(k in 1:nfolds)
{
  print(k)
  s=which(folds==k)
  X_test=X[s,]
  Y_test=Y[s,]
  X_training=X[-s,]
  Y_training=Y[-s,]
  Y_hat=diri.reg(y=Y_training, x=X_training, xnew = X_test)$est
  temp_dist=c()
  for(i in 1:nrow(Y_hat))
  {
    temp_dist[i]=comp_distance(as.vector(Y_test[i,]),as.vector(Y_hat[i,]))
  }
  error_diri[k]=sum(temp_dist)/nrow(Y_hat)
}
mean(error_diri)

# Nonlinear regression
error_nonlinear=c()
for(k in 1:nfolds)
{
  print(k)
  s=which(folds==k)
  X_test=X[s,]
  Y_test=Y[s,]
  X_training=X[-s,]
  Y_training=Y[-s,]
  Y_hat=ols.compreg(y=Y_training, x=X_training, xnew = X_test)$est
  temp_dist=c()
  for(i in 1:nrow(Y_hat))
  {
    temp_dist[i]=comp_distance(as.vector(Y_test[i,]),as.vector(Y_hat[i,]))
  }
  error_nonlinear[k]=sum(temp_dist)/nrow(Y_hat)
}
mean(error_nonlinear)

# Kullback-Leibler-divergence-based regression
error_kl=c()
for(k in 1:nfolds)
{
  print(k)
  s=which(folds==k)
  X_test=X[s,]
  Y_test=Y[s,]
  X_training=X[-s,]
  Y_training=Y[-s,]
  Y_hat=kl.compreg(y=Y_training, x=X_training, xnew = X_test)$est
  temp_dist=c()
  for(i in 1:nrow(Y_hat))
  {
    temp_dist[i]=comp_distance(as.vector(Y_test[i,]),as.vector(Y_hat[i,]))
  }
  error_kl[k]=sum(temp_dist)/nrow(Y_hat)
}
mean(error_kl)

# Alpha transformation method by Tsagris (2015)
opt_a=c()
error_alpha=c()
for(k in 1:nfolds)
{
  print(k)
  s=which(folds==k)
  X_test=X[s,]
  Y_test=Y[s,]
  X_training=X[-s,]
  Y_training=Y[-s,]
  opt_a[k]=alfareg.tune(y=Y_training, x=X_training, a = seq(-1, 1, by = 0.2), K = 10)$opt
  Y_hat=alfa.reg(y=Y_training, x=X_training, a=opt_a[k], xnew = X_test)$est
  temp_dist=c()
  for(i in 1:nrow(Y_hat))
  {
    temp_dist[i]=comp_distance(as.vector(Y_test[i,]),as.vector(Y_hat[i,]))
  }
  error_alpha[k]=sum(temp_dist)/nrow(Y_hat)
}
mean(error_alpha)