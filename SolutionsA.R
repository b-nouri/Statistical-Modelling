library(dplyr)
DataQA = read.table("D:\\covidDATA.txt")
modulo11=function(x) {x- floor(x/11)*11}
studentnumber=0772586 # fill in your student number here, this is an example!
mycol=modulo11(studentnumber)
mydataA=DataQA[,c((1+mycol*4):(4+mycol*4))]
mydataA <- mutate(mydataA,day= as.numeric(mydataA$DATE.1))
mydataA <- mutate(mydataA,y= log(mydataA$NEW_IN.1+1))
attach(mydataA)
plot(day,y)

#--------------Fit 1-------------
fit1 <- lm(y~day+pmax(day-19,0))
AIC(fit1)
plot(day,y)
lines(day,fit1$fitted.values)


#--------------Fit 2-------------
library(SemiPar)
fit2 <- spm(y~ f(day,basis="trunc.poly",degree=2),spar.method = "ML")
log.likelihood <- fit2$fit$logLik
total.df <- sum(fit2$aux$df)
AICsemipar <- -2 * (log.likelihood - total.df)
AICsemipar
plot(fit2)
points(day,y)

#--------------Fit 3-------------
fit3 <- lm(y ~ poly(day, 3))   ## polynomial of degree 3
AIC(fit3)
plot(day,y)
lines(day,fit3$fitted.values)

#--------------Fit 4-------------
library(mgcv)
fit4=gam(y~ s(day),family=gaussian(link=identity),method = "ML")
AIC(fit4)
plot(fit4)

