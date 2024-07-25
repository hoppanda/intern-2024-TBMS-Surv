
generate_data <- function(n, meanlog, sdlog){
  time <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)
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
meanlog=0
sdlog=1
data <- generate_data(n, meanlog, sdlog)

# fit null and alternmative models
fit0 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="lnorm") 
fit1 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="weibull")
fit2 <- flexsurvreg(formula=Surv(time,event)~1,data=data,,dist="gamma") 
fit3 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="llogis") 
fit4 <- flexsurvreg(formula=Surv(time,event)~1,data=data,,dist="gompertz") 
fit5 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="gengamma") 
fit6 <- flexsurvreg(formula=Surv(time,event)~1,data=data,,dist="genf") 
fit8 <-flexsurvspline(formula=Surv(time,event)~1,data=data, k=1) 
fit9 <-flexsurvspline(formula=Surv(time,event)~1,data=data, k=2) 
fit10 <-flexsurvspline(formula=Surv(time,event)~1,data=data, k=3) 


gof <-rbind(fitstats.flexsurvreg(fit0),
            fitstats.flexsurvreg(fit1),
            fitstats.flexsurvreg(fit2),
            fitstats.flexsurvreg(fit3),
            fitstats.flexsurvreg(fit4),
            fitstats.flexsurvreg(fit5),
            fitstats.flexsurvreg(fit6),
            fitstats.flexsurvreg(fit8),
            fitstats.flexsurvreg(fit9),
            fitstats.flexsurvreg(fit10))
gof$mdoel <- c("lnorm", "weibull","gamma","llogis","gompertz", "gengamma",
               "genf", "rcs1", "rcs2", "rcs3")
gof

rel_AIC_diff <- -(fit0$AIC-fit8$AIC)/fit0$AIC*100
