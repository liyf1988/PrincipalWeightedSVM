library(MASS)
library(optimx)

#svm function that return optimized coefficients
#Input: x = covariates (continuous), y = outcome (binary, {-1,1})
#Output: optimized coefficients in svm
svm.f = function(x,y,weight,lambda) {
  n = length(y)
  sigma_x = t(x)%*%x/n;    x_bar = colMeans(x)
  # weight = weight; lambda = lambda
  svm.beta = function(beta, gs=NULL) {
    tmp = sum(sapply(1:n, function(i) ((y[i]==1)-y[i]*weight)*max(1-y[i]*(beta[1]+beta[-1]%*%(x[i]-x_bar)),0)))
    Lambda = t(beta[-1])%*%sigma_x%*%beta[-1] + (lambda/n)*tmp
    return(Lambda[1,1])
  }
  svm.g = function(beta, gs=NULL) {
    tmp = matrix(NA,nrow=n,ncol=dim(x)[2]+1)
    for (i in 1:n) {
      tmp[i,1] = ((y[i]==1)-y[i]*weight)*(abs(1-y[i]*(beta[1]+beta[-1]%*%(x[i]-x_bar)))>0)*y[i]
      tmp[i,2:(dim(x)[2]+1)] = ((y[i]==1)-y[i]*weight)*(abs(1-y[i]*(beta[1]+beta[-1]%*%(x[i]-x_bar)))>0)%*%(x[i]-x_bar)*y[i]
    }
    Lambda_g = c(-(lambda/n)*sum(tmp[,1]),2*sigma_x%*%beta[-1]-(lambda/n)*colSums(tmp[,2:(dim(x)[2]+1)]))
    return(Lambda_g)
  }
  svm.h = function(beta, gs=NULL) {
    return(rbind(rep(0,dim(x)[2]+1),cbind(rep(0,dim(x)[2]),2*sigma_x)))
  }
  svm_opt = optimx(rep(0,dim(x)[2]+1),fn=svm.beta,gr=svm.g,hess=svm.h)
  return(svm_opt)
}

#Perform svm using a vector of weights, and do pca on the coefficients
#Input: x = covariates (continuous), y = outcome (binary, {-1,1})
#output: principal component analysis of the coefficients in svm, which contain the central space of y|x;
#svm coefficients estimated along a vector of weights from 0 to 1
pwsvm = function(x,y){
  beta_vector = matrix(NA,nrow=20,ncol=3)
  w_vector = seq(from=0,to=1,length=20)
  M = matrix(0,dim(x)[2],dim(x)[2])
  for (i in 1:20) {
    coef = svm.f(x,y,w_vector[i],1)[2,1:3]
    beta_vector[i,] = unlist(coef)
    M = M + unlist(coef)[-1]%*%t(unlist(coef)[-1])
  }
  pca_M = prcomp(M,scale.=F,center=F)
  return(list("pc"=pca_M,"M"=M,"beta"=beta_vector))
}


#Demonstration: generate binary y = {-1,1}, continuous x with dimension two
mu = c(0,0)
sigma = diag(2)
x0 = mvrnorm(200,mu,sigma)
x1 = mvrnorm(200,mu,sigma)+matrix(rep(2,400),nrow=200,ncol=2)
x = rbind(x0,x1)
y = c(rep(1,200),rep(-1,200))

options(warn=-1)
wpc = pwsvm(x,y)
options(warn=0)

#First component in summation of svm coefficients explains 97% of variance
first_cmp_svm = apply(x,1,function(x) wpc$pc[[2]][,1]%*%wpc$M%*%x)
#Fit y at weight=0.5
coef = unlist(svm.f(x,y,0.5,1)[2,1:3])
y_fit = (apply(x,1,function(x) coef[1]+coef[-1]%*%x)>0)
#Plot y versus first component given by weighted principal svm; colored by fitted value using weight=0.5
plot(first_cmp_svm,y,col=ifelse(y_fit==1,'red','blue'))

