################### GAMMA VS WEIBULL ####################

library(flexsurv)
library(purrr)
library(doFuture)
library(doRNG)
library(tictoc)
library(ggplot2)
library(gridExtra)

generate_data_gamma <- function(n, shape, rate, lambda){
  surv_time <- rgamma(n, shape = shape, rate = rate)
  cens_time <- rexp(n, lambda)
  time <- pmin(surv_time, cens_time)
  event <- (surv_time <= cens_time)
  data <- data.frame(time, event)
  return(data)
}

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

simu <- function(n, shape, rate, lambda, seed, R){
  
  set.seed(seed)
  
  # generate data
  data <- generate_data_gamma(n, shape, rate, lambda)
  
  # fit null and alternmative models
  fit0 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="gamma") 
  fit1 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="weibull") 
  par0 <- exp(fit0$coefficients)
  par1 <- exp(fit1$coefficients)
  
  # fit censoring model
  fit_cens <- survfit(formula=Surv(time,1-event)~1,data=data) 
  
  # Simulate bootstrap samples from null fitted model and censoring fitted model
  sim_data_surv <- simulate(fit0, nsim=R)
  
  cens_cdf <- unique(1 - c(1,fit_cens$surv))
  cens_times <- sort(data$time[data$event==0])
  cens_pmf <- c(diff(cens_cdf))
  
  # Combine together the bootstrap data
  sim_data <- sim_data_surv #just to initialize
  for(i in 1:R){
    surv_time <- sim_data_surv[[paste0("time_",i)]]
    cens_time <- sample(cens_times, n, replace = TRUE, prob = cens_pmf)
    sim_data[[paste0("time_",i)]] <- pmin(surv_time, cens_time)
    sim_data[[paste0("event_",i)]] <- (surv_time <= cens_time)
  }
  
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
  
  # compute mean adjustment
  mu <- numeric(R)
  for(i in 1:R){
    usedat <- data.frame(time=sim_data[[paste0("time_",i)]],
                         event=sim_data[[paste0("event_",i)]])
    mu[i]<-log_like(dgamma,pgamma,par0,usedat) - 
      log_like(dweibull,pweibull,star_par1,usedat)
  }
  
  # compute cox's statistic
  cox <- (lr - mean(mu)) / sd(mu)
  
  # test normality for cox's test statistic
  p_shapiro <- shapiro.test(cox)$p.value
  p_ks <- ks.test(cox, "pnorm")$p.value
  
  # compute log-LR and cox's stat for original data
  lr_obs <- fit0$loglik - fit1$loglik
  cox_obs <- (lr_obs - mean(lr)) / sd(lr)
  
  # compute p-values
  p_lr <- mean(lr <= lr_obs)
  p_cox <- mean(cox <= cox_obs)
  
  return(c(LR_p_value=p_lr,
              Cox_p_value=p_cox,
              hat_mean=mean(mu),
              hat_variance=var(mu),
              LR_obs=lr_obs,
              Cox_obs=cox_obs,
              shapiro_p_value=p_shapiro,
              ks_p_value=p_ks))
}



####### multiple repetitions of the procedure

load("seed.Rdata")
N_sim<-10000
tic()
doFuture::registerDoFuture()
doRNG::registerDoRNG()
future::plan("multisession", workers = 8)
res<- foreach(
  i=1:N_sim,
  .combine = "rbind",
  .errorhandling = "pass"
) %dorng% {
  simu(n= 300, shape = 2, rate = 1, lambda = 0.41,
       seed = seed[i], R = 100)
}
toc()


# save results
save(res, file = "simu_result_cens_E_3.Rdata")
#load("simu_result_cens_1.Rdata")

results <- as.data.frame(res)

# type I error
mean(as.numeric(results$LR_p_value[-8886])<0.05)
mean(as.numeric(results$Cox_p_value[-8886])<0.05)

# type I error normal approximation
mean(as.numeric(results$Cox_obs[-8886])<qnorm(0.05))
mean((as.numeric(results$LR_obs[-8886])-as.numeric(results$hat_mean[-8886]))/
       sqrt(as.numeric(results$hat_variance[-8886]))<qnorm(0.05))

