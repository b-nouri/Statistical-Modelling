DataQB = read.table("D:\\BankDefaultData.txt",header=T)
studentnumber = 772586 # fill in your student number here, this is an example!
set.seed(studentnumber)
rownumbers = sample(1:6436,size=1000)
mydataB = DataQB[rownumbers,]
attach(mydataB)

library(gam)

fit1 <- glm(Default~.,data=mydataB, family = binomial(link="logit"))
summary(fit1)

fit2 <- stepAIC(fit1, dirction="both",scope=list(upper=~.,lower=~1))
summary(fit2)

fit3 <- stepAIC(fit1, dirction="both",scope=list(upper=~.^2,lower=~1))
summary(fit3)


CI<- confint(fit2,level = 0.99)
CI
PostAIC(y= mydataB$Default,X=mydataB[,-3],model.set = "allsubsets", quant=0.99,
        common=c(),intercept=T,family="binomial", Pi=NULL)
