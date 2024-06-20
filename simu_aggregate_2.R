library(flexsurv)
library(purrr)
library(doFuture)
library(doRNG)
library(tictoc)

generate_data <- function(n, shape, rate, seed){
  time <- rgamma(n, shape, rate=rate)
  event <- rep(1, n)
  data <- data.frame(time, event)
  return(data)
}

log_density <- function(model, par, data) {
  d <- length(par)
  if(d==1) return(log(model(data[,1], par)))
  if(d==2) return(log(model(data[,1], par[1], par[2])))
  if(d==3) return(log(model(data[,1], par[1], par[2], par[3])))
  if(d==4) return(log(model(data[,1], par[1], par[2], par[3], par[4])))
}

log_like <- function(model, par, data) {
  # computes log-likelihood of a given model for given param and data
  return(sum(log_density(model,par,data)))
}

simu <- function(n, shape, rate, seed, R, R_low){
  
  set.seed(seed)
  
  # generate data
  data <- generate_data(n, shape, rate)
  
  # fit null and alternmative models
  fit0 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="gamma") 
  fit1 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="weibull") 
  par0 <- exp(fit0$coefficients)
  par1 <- exp(fit1$coefficients)
  
  ##### TODO: add AIC and BIC to the list returned by simu()
  
  # Simulate bootstrap samples from null fitted model
  sim_data <- simulate(fit0, nsim=R)
  
  # fit null and alternative model for each boostrap dataset
  lr<-numeric(R)
  boot_par1<-matrix(nrow=R,ncol=length(par1))
  for(i in 1:R) {
    usedata <- data.frame(time=sim_data[[paste0("time_",i)]],
                          event=sim_data[[paste0("event_",i)]])
    boot_fit0 <- flexsurvreg(formula=Surv(time,event)~1,
                             data=usedata,dist="gamma")
    boot_fit1<- flexsurvreg(formula=Surv(time,event)~1,
                            data=usedata,dist="weibull")
    lr[i]<-boot_fit0$loglik - boot_fit1$loglik
    boot_par1[i,]<-exp(boot_fit1$coefficients)
  }
  
  # compute \gamma^*
  star_par1 <- apply(boot_par1, 2, mean)
  star_par1_low <- apply(boot_par1[1:R_low, ], 2, mean)
  
  # compute mean adjustment
  mu <- numeric(R)
  for(i in 1:R){
    usedat <- data.frame(time=sim_data[[paste0("time_",i)]],
                         event=sim_data[[paste0("event_",i)]])
    mu[i]<-log_like(dgamma,par0,usedat) - 
      log_like(dweibull,star_par1,usedat)
  }
  mu_low <- numeric(R_low)
  for(i in 1:R_low){
    usedat <- data.frame(time=sim_data[[paste0("time_",i)]],
                         event=sim_data[[paste0("event_",i)]])
    mu_low[i]<-log_like(dgamma,par0,usedat) - 
      log_like(dweibull,star_par1_low,usedat)
  }
  
  # compute cox's statistic
  cox <- numeric(R)
  for(i in 1:R){
    cox[i] <- (lr[i] - mean(mu)) / (sd(mu) * sqrt(R-1))
  }
  cox_low <- numeric(R_low)
  for(i in 1:R_low){
    cox_low[i] <- (lr[i] - mean(mu_low)) / (sd(mu_low) * sqrt(R_low-1))
  }
  
  # compute log-LR and cox's stat for original data
  lr_obs <- fit0$loglik - fit1$loglik 
  mu_obs <- (fit0$loglik - log_like(dweibull, star_par1, data))
  sd_obs <- sd(log_density(dgamma,par0,data) - 
                 log_density(dweibull, star_par1, data)) * sqrt(n-1)
  cox_obs <- (lr_obs -  mu_obs) / sd_obs
  mu_obs_low <- (fit0$loglik - log_like(dweibull, star_par1_low, data))
  sd_obs_low <- sd(log_density(dgamma,par0,data) - 
                     log_density(dweibull, star_par1_low, data)) * sqrt(n-1)
  cox_obs_low <- (lr_obs -  mu_obs_low) / sd_obs_low
  
  
  # compute p-values
  p_lr <- mean(lr >= lr_obs)
  p_cox <- mean(cox >= cox_obs)
  p_cox_low <- mean(cox_low >= cox_obs_low)
  
  return(c(LR_p_value=p_lr,
           Cox_p_value=p_cox,
           Cox_low_p_value=p_cox_low,
           hat_mean=mean(mu),
           hat_variance=var(mu),
           hat_mean_low=mean(mu_low),
           hat_variance_low=var(mu_low),
           LR_obs=lr_obs,
           Cox_obs=cox_obs,
           Cox_obs_low=cox_obs_low))
}



####### multiple repetitions of the procedure

tic()
doFuture::registerDoFuture()
doRNG::registerDoRNG()
future::plan("multisession", workers = 16)
N_sim<-1000
load("seed.Rdata")
res<- foreach(
  i=1:N_sim,
  .combine = "rbind",
  .errorhandling = "pass"
) %dorng% {
  simu(n= 300, shape = 2, rate = 1, 
       seed = seed[i], R = 5000, R_low = 250)
}
toc()

# save results
save(res, file = "simu_result.Rdata")
save(seed, file = "seed.Rdata")
#load("simu_result.Rdata")

results<- as.data.frame(res)

# p-values
par(mfrow=c(1,1))
hist(results$LR_p_value)
hist(results$Cox_p_value)
hist(results$Cox_low_p_value)

#mean
mean(results$LR_obs) #sample variance of observed LR statistic
mean(results$hat_mean) #empirical mean (from bootstrap)
mean(results$hat_mean_low) #empirical mean (from bootstrap)

#variance
var(results$LR_obs) #sample variance of observed LR statistic
mean(results$hat_variance) #empirical variance (from bootstrap)
mean(results$hat_variance_low) #empirical variance (from bootstrap)
