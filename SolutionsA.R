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

#--------------Fit 1 - Piecewise linear function-------------
fit1 <- lm(y~day+pmax(day-12,0)+pmax(day-19,0))
AIC(fit1)
plot(day,y)
lines(day,fit1$fitted.values)
summary(fit1)
#--------------Fit 2 - Cubic spline with Truncated Polynomial basis Function -------------
library(SemiPar)
fit2 <- spm(y~ f(day,basis="trunc.poly",degree=3),spar.method = "ML")
log.likelihood <- fit2$fit$logLik
total.df <- sum(fit2$aux$df)
AICspm <- -2 * (log.likelihood - total.df)
AICspm

fit2$fit$coefficients
fit2$info$pen$knots
fit2$fit$coef

summary(fit2)
#--------------Fit 3 - Cubic spline with radial basis Function--------------------
fit3 <- spm(y~ f(day,basis="radial",degree=3),spar.method = "ML")
log.likelihood <- fit3$fit$logLik
total.df <- sum(fit3$aux$df)
AICspm2 <- -2 * (log.likelihood - total.df)
AICspm2

fit3$fit$coefficients
fit3$info$pen$knots
fit3$fit$coef
#------------------Plotting all the fitted means in one plot----------------------------
library(ggplot2)
ggplot(mydataA,aes(day,y)) +
  geom_point() +
  geom_line(aes(y=fit1$fitted.values,colour = "red")) +
  geom_line(aes(y=fit2$fit$fitted,colour = "blue")) +
  geom_line(aes(y=fit3$fit$fitted,colour = "green")) +
  labs(x = "Day") +
  labs(y = "Log of new cases") +
  scale_color_identity(name = "Fitted Models",
                       breaks = c("red", "blue", "green"),
                       labels = c("Piecewise linear function", "Cubic spline with Truncated Polynomial basis Function", "Cubic spline with radial basis Function"),
                       guide = "legend") +
  labs(title = "Estimate of the mean using 3 different estimators") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom")


#-----------------------Plot Fit 1 - Piecewise linear function --------------------------------------------
z1 <- predict(fit1, mydataA, se.fit = TRUE)
alpha <- 0.99
qt1 <- c(-1, 1) * qt((1 - alpha) / 2, 56-fit1$df, lower.tail = FALSE)
CI1 <- z1$fit + outer(z1$se.fit, qt1)
colnames(CI1) <- c("lower", "upper")
CI1 <- as.data.frame(CI1)

ggplot(mydataA,aes(day,y)) +
  geom_point() +
  geom_line(aes(day,fit1$fitted.values),col="green") +
  geom_line(aes(day,CI1$lower)) +
  geom_line(aes(day,CI1$upper)) +
  labs(x = "Day") +
  labs(y = "Log of new cases") +
  labs(title = "Confidence Intervals for Fit1: Piecewise linear function") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


#-----------------------Plot Fit 2 : Cubic spline with Truncated Polynomial basis Function--------------------------------------------
z2 <- predict.spm(fit2, mydataA, se = TRUE)
alpha <- 0.99
qt2 <- c(-1, 1) * qt((1 - alpha) / 2, fit2$aux$df.res, lower.tail = FALSE)
CI2 <- z2$fit + outer(z2$se, qt2)
colnames(CI2) <- c("lower", "upper")
CI2 <- as.data.frame(CI2)


ggplot(mydataA,aes(day,y)) +
  geom_point() +
  geom_line(aes(day,fit2$fit$fitted),col="pink") +
  geom_line(aes(day,CI2$lower)) +
  geom_line(aes(day,CI2$upper)) +
  labs(x = "Day") +
  labs(y = "Log of new cases") +
  labs(title = "Confidence Intervals for Fit2: Cubic spline with Truncated Polynomial basis Function") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#-----------------------Plot Fit 3 - Cubic spline with radial basis Function--------------------------------------------
z3 <- predict.spm(fit3, mydataA, se = TRUE)
alpha <- 0.99
qt3 <- c(-1, 1) * qt((1 - alpha) / 2, fit3$aux$df.res, lower.tail = FALSE)
CI3 <- z3$fit + outer(z3$se, qt3)
colnames(CI3) <- c("lower", "upper")
CI3 <- as.data.frame(CI3)

ggplot(mydataA,aes(day,y)) +
  geom_point() +
  geom_line(aes(day,fit3$fit$fitted),col="blue") +
  geom_line(aes(day,CI3$lower)) +
  geom_line(aes(day,CI3$upper)) +
  labs(x = "Day") +
  labs(y = "Log of new cases") +
  labs(title = "Confidence Intervals for Fit3: Cubic spline with radial basis Function") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))