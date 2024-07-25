hazard_Bspline <- read.csv("hazard_Bspline.csv")

# Get cdf
library(pracma)  
cumulative_hazard <- cumtrapz(hazard_Bspline$time.week, hazard_Bspline$hazard) # numerical integration
survival <- exp(-cumulative_hazard)
cdf <- 1 - survival

plot(hazard_Bspline$time.week, survival, type="s", ylim=c(0,1))

# Data
n <- nrow(hazard_Bspline)
data <- data.frame(time=sample(hazard_Bspline$time.week, n, 
                               replace = TRUE, prob = as.numeric(cdf)), 
                   event=rep(1,n))

# Estimate hazard from data
library(muhaz)
m=muhaz(data$time,data$event)
plot(m)

# Alternative 
data <- data.frame(time=hazard_Bspline$time.week, event=rep(1,n))
m=muhaz(data$time,data$event)
plot(m)





fit<-flexsurvspline(formula=Surv(time,event)~1,data=data,k=3) 
plot(fit)
plot(hazard_Bspline$time.week, survival[,1], col="green")




#################################
##### DATA

data <- data.frame(time=sample(hazard_Bspline$time.week, n, 
                               replace = TRUE, prob = as.numeric(cdf)), 
                   event=rep(1,n))

##### AIC MODEL SELECTION
fit0 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="exp") 
fit1 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="gamma") 
fit2 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="weibull")
fit3 <- flexsurvreg(formula=Surv(time,event)~1,data=data,,dist="lnorm") 
fit4 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="llogis") 
fit5 <- flexsurvreg(formula=Surv(time,event)~1,data=data,,dist="gompertz")
fit8 <- flexsurvspline(formula=Surv(time,event)~1,data=data, k=1) 
fit9 <- flexsurvspline(formula=Surv(time,event)~1,data=data, k=2) 
fit10 <- flexsurvspline(formula=Surv(time,event)~1,data=data, k=3) 

fitstats.flexsurvreg = function(x){
  ll = x$loglik
  aic = x$AIC
  k = length(x$coefficients)
  n = sum(x$data$m["(weights)"])
  aicc = aic + ((2 * k) * (k + 1) / (n - k - 1))
  bic = - 2 * ll + (k * log(n))
  data.frame(
    Df = k,
    "n2ll" = -2 * ll,
    AIC = aic,
    AICc = aicc,
    BIC = bic
  )
}

gof <-rbind(fitstats.flexsurvreg(fit0),
            fitstats.flexsurvreg(fit1),
            fitstats.flexsurvreg(fit2),
            fitstats.flexsurvreg(fit3),
            fitstats.flexsurvreg(fit4),
            fitstats.flexsurvreg(fit5),
            fitstats.flexsurvreg(fit8),
            fitstats.flexsurvreg(fit9),
            fitstats.flexsurvreg(fit10))
gof$mdoel <- c("exp", "gamma", "weibull","lnorm","llogis","gompertz", 
                "rcs1", "rcs2", "rcs3")
gof
xtable(gof)

par(mfrow=c(1,2))
plot(fit9)
plot(fit10)
