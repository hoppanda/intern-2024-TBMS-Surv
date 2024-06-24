################### GAMMA VS WEIBULL ####################

library(flexsurv)
library(purrr)
library(doFuture)
library(doRNG)
library(tictoc)
library(ggplot2)

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
    lr[i]<-(boot_fit0$loglik - boot_fit1$loglik) / n
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
    mu[i]<-(log_like(dgamma,par0,usedat) - 
              log_like(dweibull,star_par1,usedat)) / n
  }
  mu_low <- numeric(R_low)
  for(i in 1:R_low){
    usedat <- data.frame(time=sim_data[[paste0("time_",i)]],
                         event=sim_data[[paste0("event_",i)]])
    mu_low[i]<- (log_like(dgamma,par0,usedat) - 
                   log_like(dweibull,star_par1_low,usedat)) / n
  }
  
  # compute cox's statistic
  cox <- numeric(R)
  for(i in 1:R){
    cox[i] <- (lr[i] - mean(mu)) / sd(mu)
  }
  cox_low <- numeric(R_low)
  for(i in 1:R_low){
    cox_low[i] <- (lr[i] - mean(mu_low)) / sd(mu_low)
  }
  
  # test normality for cox's test statistic
  p_normal_shapiro <- shapiro.test(cox)$p.value
  p_normal_low_shapiro <- shapiro.test(cox_low)$p.value
  p_normal_ks <- ks.test(cox, "pnorm")$p.value
  p_normal_low_ks <- ks.test(cox_low, "pnorm")$p.value
  
  # compute log-LR and cox's stat for original data
  lr_obs <- (fit0$loglik - fit1$loglik) / n
  mu_obs <- (fit0$loglik - log_like(dweibull, star_par1, data)) / n
  sd_obs <- sd(log_density(dgamma,par0,data) - 
                 log_density(dweibull, star_par1, data))
  cox_obs <- (lr_obs -  mu_obs) / sd_obs
  cox_obs_2 <- (lr_obs - mean(lr)) / sd(lr)
  mu_obs_low <- (fit0$loglik - log_like(dweibull, star_par1_low, data)) / n
  sd_obs_low <- sd(log_density(dgamma,par0,data) - 
                     log_density(dweibull, star_par1_low, data))
  cox_obs_low <- (lr_obs -  mu_obs_low) / sd_obs_low
  cox_obs_low_2 <- (lr_obs -  mean(lr[1:R_low])) / sd(lr[1:R_low])
  
  # compute p-values
  p_lr <- mean(lr >= lr_obs)
  p_lr_low <- mean(lr[1:R_low] >= lr_obs)
  p_cox <- mean(cox >= cox_obs)
  p_cox_2 <- mean(cox >= cox_obs_2)
  p_cox_low <- mean(cox_low >= cox_obs_low)
  p_cox_low_2 <- mean(cox_low >= cox_obs_low_2)
  
  return(c(LR_p_value=p_lr,
           LR_low_p_value=p_lr_low,
           Cox_p_value=p_cox,
           Cox_p_value_2=p_cox_2,
           Cox_low_p_value=p_cox_low,
           Cox_low_p_value_2=p_cox_low_2,
           hat_mean=mean(mu),
           hat_variance=var(mu),
           hat_mean_low=mean(mu_low),
           hat_variance_low=var(mu_low),
           LR_obs=lr_obs,
           Cox_obs=cox_obs,
           Cox_obs_low=cox_obs_low,
           shapiro_p_value=p_normal_shapiro,
           shapiro_p_value_low=p_normal_low_shapiro,
           ks_p_value=p_normal_ks,
           ks_p_value_low=p_normal_low_ks))
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
save(res, file = "simu_result_5.Rdata")
#save(seed, file = "seed.Rdata")
#load("simu_result.Rdata")

results <- as.data.frame(res)

# LR and Cox p-values
par(mfrow=c(1,1))
hist(results$LR_p_value, main="LR p-value", 
     xlab = "p-value", xlim = c(0,1))
hist(results$LR_low_p_value, main="LR p-value R_low", 
     xlab = "p-value", xlim = c(0,1))
hist(results$Cox_p_value, main="Cox p-value R", 
     xlab = "p-value", xlim = c(0,1))
hist(results$Cox_p_value_2, main="Cox p-value", 
     xlab = "p-value", xlim = c(0,1))
hist(results$Cox_low_p_value, main="Cox p-value R_low", 
     xlab = "p-value", xlim = c(0,1))
hist(results$Cox_low_p_value_2, main="Cox p-value R_low", 
     xlab = "p-value", xlim = c(0,1))

#mean
mean(results$LR_obs) #sample mean of observed LR statistic
mean(results$hat_mean) #empirical mean (from bootstrap)
mean(results$hat_mean_low) #empirical mean (from bootstrap)

#variance
var(results$LR_obs) #sample variance of observed LR statistic
mean(results$hat_variance) #empirical variance (from bootstrap)
mean(results$hat_variance_low) #empirical variance (from bootstrap)

# test normality
hist(results$shapiro_p_value, main="Shapiro test p-value", 
     xlab = "p-value", xlim = c(0,1))
hist(results$shapiro_p_value_low, main="Shaprio test p-value R_low", 
     xlab = "p-value", xlim = c(0,1))
hist(results$ks_p_value, main="KS test p-value", 
     xlab = "p-value", xlim = c(0,1))
hist(results$ks_p_value_low, main="KS test p-value R_low", 
     xlab = "p-value", xlim = c(0,1))

# test normality 1 simulation
#qqnorm(cox)
#qqline(cox)
#hist(cox, breaks = 50)
ggplot(as.data.frame(cox), aes(x = cox)) + geom_histogram()
ggplot(as.data.frame(cox), aes(x = cox)) + geom_density()
ggplot(as.data.frame(cox), aes(sample = cox)) + geom_qq() + geom_qq_line()
shapiro.test(cox)
ks.test(cox, "pnorm")

