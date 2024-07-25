library(flexsurv)

log_like_contribution <- function(model_density, model_cdf, par, data, i) {
  if(data$event[i]) {
    return(log(model_density(data$time[i], par[1], par[2])))
  } else {
    return(log(1 - model_cdf(data$time[i], par[1], par[2])))
  }
}

log_like <- function(model_density, model_cdf, par, data) {
  # computes log-likelihood of a given model for given param and data
  n <- nrow(data)
  return(sum(sapply(1:n, log_like_contribution, 
                    model_density = model_density, 
                    model_cdf = model_cdf, 
                    par = par, 
                    data = data)))
}
log_like_rcs <- function(coeff, knots, data) {
  ll <- NULL
  for (i in 1:nrow(data)){
    if(data$event[i]) {
      ll[i] <- log(dsurvspline(data$time[i], gamma=coeff, knots=knots))
    } else {
      ll[i] <- log(1 - psurvspline(data$time[i], gamma=coeff, knots=knots))
    }
  }
  return(sum(ll))
}

##### DATA

data <- read.csv("DANUBE_1.csv")
data <- data[,-c(1,4)]


##### AIC MODEL SELECTION
fit0 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="exp") 
fit1 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="gamma") 
fit2 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="weibull")
fit3 <- flexsurvreg(formula=Surv(time,event)~1,data=data,,dist="lnorm") 
fit4 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="llogis") 
fit5 <- flexsurvreg(formula=Surv(time,event)~1,data=data,,dist="gompertz") 
fit6 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="gengamma") 
fit7 <- flexsurvreg(formula=Surv(time,event)~1,data=data,,dist="genf") 
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
    AIC = aic,
    BIC = bic
  )
}

gof <-rbind(fitstats.flexsurvreg(fit0),
            fitstats.flexsurvreg(fit1),
            fitstats.flexsurvreg(fit2),
            fitstats.flexsurvreg(fit3),
            fitstats.flexsurvreg(fit4),
            fitstats.flexsurvreg(fit5),
            fitstats.flexsurvreg(fit6),
            fitstats.flexsurvreg(fit7),
            fitstats.flexsurvreg(fit8),
            fitstats.flexsurvreg(fit9),
            fitstats.flexsurvreg(fit10))
model <- c("exp", "gamma", "weibull","lnorm","llogis","gompertz", 
               "gengamma", "genf", "rcs1", "rcs2", "rcs3")
gof <- cbind(model, gof)
gof

library(xtable)
xtable(gof)



##### PROCEDURE LOG-NORMAL vs FRC3

R <- 1000
fit0 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="lnorm")
fit1 <- flexsurvspline(formula=Surv(time,event)~1,data=data, k=3)
par0 <- c(fit0$coefficients[1], exp(fit0$coefficients[2]))
par1 <-  fit1$coefficients
knots1 <- fit1$knots

par(mfrow=c(1,2))
plot(fit0)
plot(fit1)

sim_data <- simulate(fit0, nsim=R, seed=12345, censtime = max(data$time))

lr<-numeric(R)
boot_par1<-matrix(nrow=R,ncol=length(par1))
boot_knots1<-matrix(nrow=R,ncol=length(knots1))
aic <- NULL
for(i in 1:R) {
  usedata <- data.frame(time=sim_data[[paste0("time_",i)]],
                        event=sim_data[[paste0("event_",i)]])
  boot_fit0 <- flexsurvreg(formula=Surv(time,event)~1,
                           data=usedata, dist="lnorm")
  boot_fit1 <- flexsurvspline(formula=Surv(time,event)~1,
                              data=usedata,k=3)
  lr[i]<-boot_fit0$loglik - boot_fit1$loglik
  boot_par1[i,]<-boot_fit1$coefficients
  boot_knots1[i,]<- boot_fit1$knots
  aic[i] <- boot_fit0$AIC<=boot_fit1$AIC
}

# compute \gamma^*
star_par1 <- apply(boot_par1, 2, mean)
star_knots1 <- apply(boot_knots1, 2, mean)

# compute mean adjustment
mu <- numeric(R)
for(i in 1:R){
  usedat <- data.frame(time=sim_data[[paste0("time_",i)]],
                       event=sim_data[[paste0("event_",i)]])
  mu[i]<- log_like(dlnorm,plnorm,par0,usedat) -
    log_like_rcs(star_par1,star_knots1,usedat) 
}

#check
mean(mu)
mean(lr)

#compute cox's statistic
cox <- (lr - mean(mu)) / sd(mu)
hist(cox)

# test normality for cox's test statistic
#p_shapiro <- shapiro.test(cox)$p.value
#p_ks <- ks.test(cox, "pnorm")$p.value

# compute log-LR and cox's stat for original data
lr_obs <- fit0$loglik - fit1$loglik
cox_obs <- (lr_obs - mean(lr)) / sd(lr)

# compute p-values
p_lr <- mean(lr <= lr_obs)
p_cox <- mean(cox <= cox_obs)

# normal approx
pnorm(cox_obs) #tilde
pnorm((lr_obs-mean(mu))/sd(mu)) #hat



