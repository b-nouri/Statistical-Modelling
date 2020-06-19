DataQC = read.table("D:\\bikestations.txt")
digitsum = function(x) sum(floor(x/10^(0:(nchar(x)-1)))%%10)
studentnumber=772586 # fill in your student number here, this is an example!
mysum = digitsum(studentnumber)
Y = DataQC[,mysum]
X = DataQC[,-mysum]

library(glmnet)
X=as.matrix(X)
#-------------------Ridge Regression-------------------------------------------
fit1=cv.glmnet(X,Y, family = "poisson", nfold=7, alpha=0)
fit1$lambda.min
coef(fit1, s=fit1$lambda.min)
print(fit1)

#-------------------Lasso Regression-------------------------------------------
fit2=cv.glmnet(X,Y, family = "poisson", nfold=7, alpha=1)
fit2$lambda.min
coef(fit2, s=fit2$lambda.min)
print(fit2)

#-------------------elastic net 1 Regression-----------------------------------
fit3=cv.glmnet(X,Y, family = "poisson", nfold=7, alpha=0.6)
fit3$lambda.min
coef(fit3, s=fit3$lambda.min)
print(fit3)
#-------------------elastic net 2 Regression-----------------------------------
fit4=cv.glmnet(X,Y, family = "poisson", nfold=7, alpha=0.1)
fit4$lambda.min
coef(fit4, s=fit4$lambda.min)
print(fit4)
