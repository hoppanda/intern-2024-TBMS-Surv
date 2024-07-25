################### GAMMA VS WEIBULL ####################

library(flexsurv)
library(purrr)
library(doFuture)
library(doRNG)
library(tictoc)
library(ggplot2)
library(gridExtra)

n=300
shape=2
rate=1
lambda=0.1
seed = 230301
R = 100


generate_data <- function(n, shape, rate, lambda){
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

simu <- function(n, shape, rate, lambda, seed, R, R_a, R_b){
  
  set.seed(seed)
  
  # generate data
  data <- generate_data(n, shape, rate, lambda)
  max_time <- max(data$time)
  
  # fit null and alternmative models
  fit0 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="gamma") 
  fit1 <- flexsurvreg(formula=Surv(time,event)~1,data=data,dist="weibull") 
  par0 <- exp(fit0$coefficients)
  par1 <- exp(fit1$coefficients)
  
  ##### TODO: add AIC and BIC to the list returned by simu()
  
  # Simulate bootstrap samples from null fitted model
  sim_data <- simulate(fit0, nsim=R, censtime = max_time)
  
  # fit null and alternative model for each boostrap dataset
  lr<-numeric(R)
  boot_par0<-matrix(nrow=R,ncol=length(par0))
  boot_par1<-matrix(nrow=R,ncol=length(par1))
  for(i in 1:R) {
    usedata <- data.frame(time=sim_data[[paste0("time_",i)]],
                          event=sim_data[[paste0("event_",i)]])
    boot_fit0 <- flexsurvreg(formula=Surv(time,event)~1,
                             data=usedata,dist="gamma")
    boot_fit1<- flexsurvreg(formula=Surv(time,event)~1,
                            data=usedata,dist="weibull")
    lr[i]<-boot_fit0$loglik - boot_fit1$loglik
    boot_par0[i,]<-exp(boot_fit0$coefficients)
    boot_par1[i,]<-exp(boot_fit1$coefficients)
  }
  
  # compute \gamma^*
  star_par1 <- apply(boot_par1, 2, mean)
  star_par1_a <- apply(boot_par1[1:R_a, ], 2, mean)
  star_par1_b <- apply(boot_par1[1:R_b, ], 2, mean)
  
  # compute mean adjustment
  mu <- numeric(R)
  for(i in 1:R){
    usedat <- data.frame(time=sim_data[[paste0("time_",i)]],
                         event=sim_data[[paste0("event_",i)]])
    mu[i]<-log_like(dgamma,pgamma,par0,usedat) - 
      log_like(dweibull,pweibull,star_par1,usedat)
  }
  mu_a <- numeric(R_a)
  for(i in 1:R_a){
    usedat <- data.frame(time=sim_data[[paste0("time_",i)]],
                         event=sim_data[[paste0("event_",i)]])
    mu_a[i]<- log_like(dgamma,pgamma,par0,usedat) - 
      log_like(dweibull,pweibull,star_par1_a,usedat)
  }
  mu_b <- numeric(R_b)
  for(i in 1:R_b){
    usedat <- data.frame(time=sim_data[[paste0("time_",i)]],
                         event=sim_data[[paste0("event_",i)]])
    mu_b[i]<- log_like(dgamma,pgamma,par0,usedat) - 
      log_like(dweibull,pweibull,star_par1_b,usedat)
  }
  
  # compute cox's statistic
  cox <- (lr - mean(mu)) / sd(mu)
  cox_a <- (lr[1:R_a] - mean(mu_a)) / sd(mu_a)
  cox_b <- (lr[1:R_b] - mean(mu_b)) / sd(mu_b)
  
  # test normality for cox's test statistic
  p_shapiro <- shapiro.test(cox)$p.value
  p_shapiro_a <- shapiro.test(cox_a)$p.value
  p_shapiro_b <- shapiro.test(cox_b)$p.value
  p_ks <- ks.test(cox, "pnorm")$p.value
  p_ks_a <- ks.test(cox_a, "pnorm")$p.value
  p_ks_b <- ks.test(cox_b, "pnorm")$p.value
  
  # compute log-LR and cox's stat for original data
  lr_obs <- fit0$loglik - fit1$loglik
  cox_obs <- (lr_obs - mean(lr)) / sd(lr)
  cox_obs_a <- (lr_obs -  mean(lr[1:R_a])) / sd(lr[1:R_a])
  cox_obs_b <- (lr_obs -  mean(lr[1:R_b])) / sd(lr[1:R_b])
  
  # compute p-values
  p_lr <- mean(lr <= lr_obs)
  p_lr_a <- mean(lr[1:R_a] <= lr_obs)
  p_lr_b <- mean(lr[1:R_b] <= lr_obs)
  p_cox <- mean(cox <= cox_obs)
  p_cox_a <- mean(cox_a <= cox_obs_a)
  p_cox_b <- mean(cox_b <= cox_obs_b)
  
  return(c(LR_p_value=p_lr,
              LR_p_value_a=p_lr_a,
              LR_p_value_b=p_lr_b,
              Cox_p_value=p_cox,
              Cox_p_value_a=p_cox_a,
              Cox_p_value_b=p_cox_b,
              hat_mean=mean(mu),
              hat_variance=var(mu),
              hat_mean_a=mean(mu_a),
              hat_variance_a=var(mu_a),
              hat_mean_b=mean(mu_b),
              hat_variance_b=var(mu_b),
              LR_obs=lr_obs,
              Cox_obs=cox_obs,
              Cox_obs_a=cox_obs_a,
              Cox_obs_b=cox_obs_b,
              shapiro_p_value=p_shapiro,
              shapiro_p_value_a=p_shapiro_a,
              shapiro_p_value_b=p_shapiro_b,
              ks_p_value=p_ks,
              ks_p_value_a=p_ks_a,
              ks_p_value_b=p_ks_b))
}



####### multiple repetitions of the procedure

load("seed.Rdata")
N_sim<-1000

tic()
doFuture::registerDoFuture()
doRNG::registerDoRNG()
future::plan("multisession", workers = 8)
res<- foreach(
  i=1:N_sim,
  .combine = "rbind",
  .errorhandling = "pass"
) %dorng% {
  simu(n= 300, shape = 2, rate = 1, lambda=0.1,
       seed = seed[i], R = 250, R_a = 150, R_b = 100)
}
toc()

# save results
save(res, file = "simu_result_cens_admin_1.Rdata")
#load("simu_result_9.Rdata")

results <- as.data.frame(res)

# mean
mean(results$LR_obs) #sample mean of observed LR statistic
mean(results$hat_mean) #empirical mean (from bootstrap)
mean(results$hat_mean_a) #empirical mean (from bootstrap)
mean(results$hat_mean_b) #empirical mean (from bootstrap)

# variance
var(results$LR_obs) #sample variance of observed LR statistic
mean(results$hat_variance) #empirical variance (from bootstrap)
mean(results$hat_variance_a) #empirical variance (from bootstrap)
mean(results$hat_variance_b) #empirical variance (from bootstrap)

# type I error
mean(results$LR_p_value<0.05)*100
mean(results$LR_p_value_a<0.05)*100
mean(results$LR_p_value_b<0.05)*100
mean(results$Cox_p_value<0.05)*100
mean(results$Cox_p_value_a<0.05)*100
mean(results$Cox_p_value_b<0.05)*100

# type I error normal approximation
mean(results$Cox_obs<qnorm(0.05))*100
mean(results$Cox_obs_a<qnorm(0.05))*100
mean(results$Cox_obs_b<qnorm(0.05))*100
mean((results$LR_obs-results$hat_mean)/
       sqrt(results$hat_variance)<qnorm(0.05))*100
mean((results$LR_obs-results$hat_mean_a)/
       sqrt(results$hat_variance_a)<qnorm(0.05))*100
mean((results$LR_obs-results$hat_mean_b)/
       sqrt(results$hat_variance_b)<qnorm(0.05))*100

# p-value normal approximation
p_norm_approx <- pnorm(results$Cox_obs)
hist(p_norm_approx)
mean(p_norm_approx<0.05)

p_norm_approx_2 <- pnorm((results$LR_obs-results$hat_mean)/
                           sqrt(results$hat_variance))
hist(p_norm_approx_2)
mean(p_norm_approx_2<0.05) 

# LR and Cox p-values
par(mfrow=c(1,1))
hist(results$LR_p_value, main="LR p-value R", 
     xlab = "p-value", xlim = c(0,1))
hist(results$LR_p_value_a, main="LR p-value R_a", 
     xlab = "p-value", xlim = c(0,1))
hist(results$LR_p_value_b, main="LR p-value R_b", 
     xlab = "p-value", xlim = c(0,1))
hist(results$Cox_p_value, main="Cox p-value R", 
     xlab = "p-value", xlim = c(0,1))
hist(results$Cox_p_value_a, main="Cox p-value R_a", 
     xlab = "p-value", xlim = c(0,1))
hist(results$Cox_p_value_b, main="Cox p-value R_b", 
     xlab = "p-value", xlim = c(0,1))

# test normality
hist(results$shapiro_p_value, main="Shapiro test p-value", 
     xlab = "p-value", xlim = c(0,1))
hist(results$shapiro_p_value_a, main="Shaprio test p-value R_a", 
     xlab = "p-value", xlim = c(0,1))
hist(results$shapiro_p_value_b, main="Shaprio test p-value R_b", 
     xlab = "p-value", xlim = c(0,1))
hist(results$ks_p_value, main="KS test p-value", 
     xlab = "p-value", xlim = c(0,1))
hist(results$ks_p_value_a, main="KS test p-value R_a", 
     xlab = "p-value", xlim = c(0,1))
hist(results$ks_p_value_b, main="KS test p-value R_b", 
     xlab = "p-value", xlim = c(0,1))

# parameters
par(mfrow=c(1,2))
obs_par0_shape=as.vector(sapply(results$par0_obs, function(x) x[1]))
boot_par0_shape=as.vector(sapply(results$boot_par0, function(x) x[,1]))
hist(boot_par0_shape, freq = FALSE, main = "Histogram of shape for the N simulations",
     ylim = c(0, max(density(boot_par0_shape)$y, density(obs_par0_shape)$y)),
     col = "lightblue",)
hist(obs_par0_shape, freq = FALSE, add=TRUE, col = rgb(1, 0, 0, alpha = 0.5))
legend("topright", legend = c("boostrap", "observed"), fill = c("lightblue", rgb(1, 0, 0, alpha = 0.5)))

obs_par0_rate=as.vector(sapply(results$par0_obs, function(x) x[2]))
boot_par0_rate=as.vector(sapply(results$boot_par0, function(x) x[,2]))
hist(boot_par0_rate, freq = FALSE, main = "Histogram of rate for the N simulations",
     ylim = c(0, 4), col = "lightblue",)
hist(obs_par0_rate, freq = FALSE, add=TRUE, col = rgb(1, 0, 0, alpha = 0.5))
legend("topright", legend = c("boostrap", "observed"), fill = c("lightblue", rgb(1, 0, 0, alpha = 0.5)))

# test normality 1 simulation
p1=ggplot(as.data.frame(cox), aes(x = cox)) +
  geom_histogram(aes(y = ..density..))+ geom_density()
p2=ggplot(as.data.frame(cox), aes(sample = cox)) + geom_qq() + geom_qq_line()
grid.arrange(p1, p2, ncol = 2)

p1=ggplot(as.data.frame(cox_low), aes(x = cox_low)) +
  geom_histogram(aes(y = ..density..), bins = 15)+ geom_density()
p2=ggplot(as.data.frame(cox_low), aes(sample = cox_low)) + geom_qq() + geom_qq_line()
grid.arrange(p1, p2, ncol = 2)


shapiro.test(cox)
ks.test(cox, "pnorm")


# plots p-values
plot(results$LR_p_value, results$Cox_p_value, xlab="LR p-value", 
     ylab = "Cox p-value", main="Cox vs LR p-values with R")

# hazard shape
curve(dgamma(x,2,1)/(1-pgamma(x,2,1)), xlim=c(0,10), 
      ylab = "h(x)", main = "Hazard function of Gamma(2,1)")
