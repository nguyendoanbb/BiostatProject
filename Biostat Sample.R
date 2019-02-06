############################################################################################
#set up directory
setwd("~/Dropbox/Work Sample") #change this to directory containing dataset
dat <- read.csv("Data for sample_cgd.csv")
library(MASS)
library(faraway)

#rename variables
names(dat)[names(dat)=='Z1'] <- 'treatment'
names(dat)[names(dat)=='Z2'] <- 'pat.inhe'
names(dat)[names(dat)=='Z3'] <- 'age'
names(dat)[names(dat)=='Z4'] <- 'height'
names(dat)[names(dat)=='Z5'] <- 'weight'
names(dat)[names(dat)=='Z6'] <- 'cortico'
names(dat)[names(dat)=='Z7'] <- 'prophy'
names(dat)[names(dat)=='Z8'] <- 'gender'
names(dat)[names(dat)=='Z9'] <- 'hospital'
names(dat)[names(dat)=='T1'] <- 'time'
names(dat)[names(dat)=='D'] <- 'censored'

#change variables' type
dat$treatment <- as.factor(dat$treatment)
dat$pat.inhe <- as.factor(dat$pat.inhe)
dat$cortico <- as.factor(dat$cortico)
dat$prophy <- as.factor(dat$prophy)
dat$gender <- as.factor(dat$gender)
dat$hospital <- as.factor(dat$hospital)
dat$censored <- as.factor(dat$censored)

#final dataset in use
dat <- dat[,-c(1:3,14,16)]

#format censored data
dat$censored <- ifelse(dat$censored == 2, 0, 1)

#descriptive table
summary(dat)

#library setup
library("survminer")
library(survival)

#proportional hazard model
dat$censored <- as.numeric(dat$censored)
cox.fit <- coxph(Surv(time, censored)~., data=dat)
cox.fit4 <- step(cox.fit)
summary(cox.fit4)

plot(survfit(coxph(Surv(time,censored)~strata(treatment), data=dat)))
library(rms)
survplot(npsurv(coxph(Surv(time,censored)~strata(treatment), data=dat)),
         col = c(1,2), lty = c(1,2))

#PH model diagnostic
#checking PH assumption
test.ph <- cox.zph(cox.fit4)
test.ph #no significant p-value, hence PH assumption is satisfied
ggcoxzph(test.ph)

#checking outliers/influential points
ggcoxdiagnostics(cox.fit4, type = 'deviance' , linear.predictions = FALSE,
                 ggtheme = theme_bw(), title = 'Checking Outliers or Influential Observations by Deviance Residuals') #no visible outliers

#checking linearity of continuous variables
plot(dat$age, residuals(cox.fit4, type='martingale'), xlab = 'Age',
     ylab='Martingale residual')
lines(lowess(dat$age, residuals(cox.fit4, type='martingale'), iter=0))
plot(dat$weight, residuals(cox.fit4, type='martingale'), xlab = 'Weight',
     ylab='Martingale residual')
lines(lowess(dat$weight, residuals(cox.fit4, type='martingale'), iter=0))

#testing for interaction with time
coxph(Surv(time, censored)~treatment + cortico + age + weight +
        tt(age) + tt(weight), data=dat)


#parametric model: assume exponential function of time to disease
par.fit <- WeibullReg(Surv(time, censored)~treatment + age + weight + cortico, data=dat)
par.fit
par.fit1 <- survreg(Surv(time, censored)~treatment,data=dat)

#median follow-up time
predict(par.fit1, newdata = list(treatment = c('1','2')), type = 'quantile', p = 0.5)

#plot survival curve using parametric model
plot(predict(par.fit1, newdata=list(treatment='1'),type="quantile",p=seq(.01,.99,by=.01)),
     seq(.99,.01,by=-.01),col="red", type = 'n',
     ylab = 'Survival Probability', xlab = 'Follow-up Time')
lines(predict(par.fit1, newdata=list(treatment='1'),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col=1, lty=1)
lines(predict(par.fit1, newdata=list(treatment='2'),type="quantile",p=seq(.01,.99,by=.01)),seq(.99,.01,by=-.01),col=2, lty=2)
legend(1700,0.8, c('Gamma', 'Placebo'), lty=c(1:2), col=c(1:2),
       cex=0.75)
##########################################################################################################

