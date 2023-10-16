### MATH 4130k Project ###

#heart.data <- read.csv("/Users/Vyene Guo/OneDrive - York University/Desktop/YorkU/MATH 4130K/Project/heart_failure_clinical_records_dataset.csv", header=TRUE)
heart.data <-read.csv("heart_failure_clinical_records_dataset.csv")
# Rename some variables.
names(heart.data)[names(heart.data) == "creatinine_phosphokinase"] <- "CPK"
names(heart.data)[names(heart.data) == "ejection_fraction"] <- "EF"
names(heart.data)[names(heart.data) == "high_blood_pressure"] <- "HBP"
names(heart.data)[names(heart.data) == "serum_creatinine"] <- "SCR"
names(heart.data)[names(heart.data) == "serum_sodium"] <- "SSO"
names(heart.data)[names(heart.data) == "DEATH_EVENT"] <- "status"

# See how many categorical variables
colSums(heart.data==0)

# Turn variable 'EF' into categorical variable
library(dplyr)
heart.data$EF <- cut(heart.data$EF, breaks=c(10, 40, 55, 70, 80), 
                     labels=c('1','2','3','4'))
heart.data$EF <- unclass(heart.data$EF)

# Standardize all the continuous variables
heart.std <- heart.data %>% mutate_at(c('age','CPK','platelets','SCR','SSO'), 
                                      ~(scale(.) %>% as.vector))
head(heart.std)

# The overall survival curve by KM estimators
library(survival)
fit <- survfit(Surv(heart.std$time, heart.std$status)~1)
nj <- fit$n.risk
dj <- fit$n.event
## Kaplan-Meier estimator ##
km <- cumprod(1-dj/nj)
v.km <- km^2*cumsum(dj/(nj*(nj-dj))) # Estimating the variance
km <- c(1, km) # Concatenating a 1 before the first event
# Event times
etimes <- c(0, fit$time)
v.km <- c(0, v.km)
# Plotting Kaplan-Meier curve --> an overall picture
plot(etimes, km, type="s", ylim=c(0,1), 
     xlab="Days", ylab="Survival Probability",
     col="firebrick", lwd=2)
# Adding confidence intervals
lines(etimes, km-1.96*sqrt(v.km), type="s", lty=2, col="darkblue")
lines(etimes, km+1.96*sqrt(v.km), type="s", lty=2, col="darkblue")
abline(h=0.6, col="yellowgreen", lwd=2)
title(main="Kaplan-Meier Survival Curve")

# Cox proportional hazards model (overall)
cph.fit <- coxph(Surv(time,status) ~ 
                   age+CPK+EF+platelets+SCR+SSO+anaemia
                 +diabetes+HBP+sex+smoking, data=heart.std)
cph.fit
#plot(survfit(cph.fit), ylim=c(0,1), xlab="Days", ylab="S(t)")
#abline(h=0.6, col="green")


library(gtsummary)  #install.packages("gtsummary")
tbl_regression(cph.fit, exponentiate = TRUE)

# Applying LASSO regression to select significant variables
library(glmnet)
X <- model.matrix(cph.fit)
y <- Surv(heart.std$time, heart.std$status)
set.seed(20230414)
cv.lasso <- cv.glmnet(X, y, alpha=1, family="cox", nfolds=10)
best.lambda <- cv.lasso$lambda.min
coef <- coef(cv.lasso, s=best.lambda)
coef

# Reduced model: only includes covariates 
# 'age', 'CPK', 'EF', 'SCR', 'SSO', 'anaemia' and 'HBP'
cph.fit2 <- coxph(Surv(time, status) ~ age + CPK + EF + SCR + SSO 
                  + anaemia + HBP, 
                  data = heart.std)
cph.fit2
tbl_regression(cph.fit2, exponentiate = TRUE)

# Create a survival object for HBP = 0 and anaemia = 0
sf00 <- survfit(cph.fit2, newdata = data.frame(HBP = 0, anaemia = 0,
                                              age = mean(heart.std$age),
                                              CPK = mean(heart.std$CPK),
                                              EF = mean(heart.std$EF),
                                              SCR = mean(heart.std$SCR), 
                                              SSO = mean(heart.std$SSO)))
# Create a survival object for HBP = 1 and anaemia = 0
sf10 <- survfit(cph.fit2, newdata = data.frame(HBP = 1, anaemia = 0,
                                              age = mean(heart.std$age),
                                              CPK = mean(heart.std$CPK),
                                              EF = mean(heart.std$EF),
                                              SCR = mean(heart.std$SCR), 
                                              SSO = mean(heart.std$SSO)))
# Create a survival object for HBP = 0 and anaemia = 1
sf01 <- survfit(cph.fit2, newdata = data.frame(HBP = 0, anaemia = 1,
                                               age = mean(heart.std$age),
                                               CPK = mean(heart.std$CPK),
                                               EF = mean(heart.std$EF),
                                               SCR = mean(heart.std$SCR), 
                                               SSO = mean(heart.std$SSO)))
# Create a survival object for HBP = 1 and anaemia = 0
sf11 <- survfit(cph.fit2, newdata = data.frame(HBP = 1, anaemia = 1,
                                               age = mean(heart.std$age),
                                               CPK = mean(heart.std$CPK),
                                               EF = mean(heart.std$EF),
                                               SCR = mean(heart.std$SCR), 
                                               SSO = mean(heart.std$SSO)))
# Plot the survival curves
plot(sf00, col = "yellowgreen", lty = 1, ylim = c(0.3, 1), xlab = "Days", 
     ylab = "Survival Probability", lwd=2, conf.int=FALSE)
lines(sf10, col = "firebrick", lty = 1, lwd=2, conf.int=FALSE)
lines(sf01, col = "navy", lty = 1, lwd=2, conf.int=FALSE)
lines(sf11, col = "black", lty = 1, lwd=2, conf.int=FALSE)
legend("bottomleft", legend = c("HBP = 0, anaemia = 0", "HBP = 1", 
                                "anaemia = 1",
                                "HBP = 1, anaemia = 1"), 
       col = c("yellowgreen", "firebrick", "navy", "black"), 
       lty = 1, cex=0.45, bg="white", 
       lwd=2)
title(main="Comparing Survival Curves With Levels of HBP and Anaemia", 
      line=2.5)



# Create a survival object for EF = 3
sf1 <- survfit(cph.fit2, newdata = data.frame( EF = 3,
                                               HBP = 0, anaemia = 0,
                                               age = mean(heart.std$age),
                                               CPK = mean(heart.std$CPK),
                                               SCR = mean(heart.std$SCR), 
                                               SSO = mean(heart.std$SSO)))
# Create a survival object for EF = 2
sf2 <- survfit(cph.fit2, newdata = data.frame( EF = 2,
                                               HBP = 0, anaemia = 0,
                                               age = mean(heart.std$age),
                                               CPK = mean(heart.std$CPK),
                                               SCR = mean(heart.std$SCR), 
                                               SSO = mean(heart.std$SSO)))
# Create a survival object for EF = 1
sf3 <- survfit(cph.fit2, newdata = data.frame( EF = 1,
                                               HBP = 0, anaemia = 0,
                                               age = mean(heart.std$age),
                                               CPK = mean(heart.std$CPK),
                                               SCR = mean(heart.std$SCR), 
                                               SSO = mean(heart.std$SSO)))
# Plot the survival curves
plot(sf1, col = "yellowgreen", lty = 1, ylim = c(0.3, 1), xlab = "Days", 
     ylab = "Survival Probability", lwd=2, conf.int=FALSE)
lines(sf2, col = "firebrick", lty = 1, lwd=2, conf.int=FALSE)
lines(sf3, col = "black", lty = 1, lwd=2, conf.int=FALSE)
legend("bottomleft", 
       legend = c("55 < EF <=70", "40 < EF <= 55", "EF <= 40"), 
       col = c("yellowgreen", "firebrick", "black"), 
       lty = 1, cex=1, bg="white", 
       lwd=2)
title(main="Comparing Survival Curves With Levels of EF", line=2.5)




# Testing proportional hazards assumption
cz <- cox.zph(cph.fit2)
cz
plot(cz)


# Check assumption for high blood pressure 
# (by create KM curves separately for HBP)

km.hbp <- survfit(Surv(time, status)~HBP, data=heart.std)
plot(km.hbp, fun="cloglog")
title(main="log-log plot for HBP", line=2.5)

# Check assumption for anaemia
km.anaemia <- survfit(Surv(time, status)~anaemia, data=heart.std)
plot(km.anaemia, fun="cloglog")
title(main="log-log plot for Anaemia", line=2.5)

