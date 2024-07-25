library(flexsurv)

generate_data <- function(n, shape, rate, seed){
  time <- rgamma(n, shape, rate=rate)
  event <- rep(1, n)
  data <- data.frame(time, event)
  return(data)
}

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


set.seed(120798)

# generate data
n=300
shape=2
rate=1
data <- generate_data(n, shape, rate)

# fit null and alternmative models
fit0 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="gamma") 
fit1 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="weibull")
fit2 <- flexsurvreg(formula=Surv(time,event)~1,data=data,,dist="lnorm") 
fit3 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="llogis") 
fit4 <- flexsurvreg(formula=Surv(time,event)~1,data=data,,dist="gompertz")
fit7 <-flexsurvspline(formula=Surv(time,event)~1,data=data, k=1) 
fit8 <-flexsurvspline(formula=Surv(time,event)~1,data=data, k=2) 
fit9 <-flexsurvspline(formula=Surv(time,event)~1,data=data, k=3) 

gof <-rbind(fitstats.flexsurvreg(fit0),
            fitstats.flexsurvreg(fit1),
            fitstats.flexsurvreg(fit2),
            fitstats.flexsurvreg(fit3),
            fitstats.flexsurvreg(fit4),
            fitstats.flexsurvreg(fit7),
            fitstats.flexsurvreg(fit8),
            fitstats.flexsurvreg(fit9))
gof$mdoel <- c("gamma", "weibull","lnorm","llogis","gompertz", 
               "rcs1", "rcs2", "rcs3")
gof

xtable(gof)

rel_AIC_diff <- -(fit0$AIC-fit1$AIC)/fit0$AIC*100
