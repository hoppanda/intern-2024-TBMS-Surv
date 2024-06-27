################### GAMMA VS WEIBULL ####################

library(flexsurv)
library(purrr)
library(doFuture)
library(doRNG)
library(tictoc)
library(ggplot2)
library(gridExtra)

generate_data <- function(n, shape, rate){
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

simu <- function(n, shape, rate, seed, R, R_a, R_b){
  
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
  star_par1_a <- apply(boot_par1[1:R_a, ], 2, mean)
  star_par1_b <- apply(boot_par1[1:R_b, ], 2, mean)
  
  # compute mean adjustment
  mu <- numeric(R)
  for(i in 1:R){
    usedat <- data.frame(time=sim_data[[paste0("time_",i)]],
                         event=sim_data[[paste0("event_",i)]])
    mu[i]<-(log_like(dgamma,par0,usedat) - 
              log_like(dweibull,star_par1,usedat)) / n
  }
  mu_a <- numeric(R_a)
  for(i in 1:R_a){
    usedat <- data.frame(time=sim_data[[paste0("time_",i)]],
                         event=sim_data[[paste0("event_",i)]])
    mu_a[i]<- (log_like(dgamma,par0,usedat) - 
                   log_like(dweibull,star_par1_a,usedat)) / n
  }
  mu_b <- numeric(R_b)
  for(i in 1:R_b){
    usedat <- data.frame(time=sim_data[[paste0("time_",i)]],
                         event=sim_data[[paste0("event_",i)]])
    mu_b[i]<- (log_like(dgamma,par0,usedat) - 
                   log_like(dweibull,star_par1_b,usedat)) / n
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
  lr_obs <- (fit0$loglik - fit1$loglik) / n
  cox_obs <- (lr_obs - mean(lr)) / sd(lr)
  cox_obs_a <- (lr_obs -  mean(lr[1:R_a])) / sd(lr[1:R_a])
  cox_obs_b <- (lr_obs -  mean(lr[1:R_b])) / sd(lr[1:R_b])
  
  # compute p-values
  p_lr <- mean(lr >= lr_obs)
  p_lr_a <- mean(lr[1:R_a] >= lr_obs)
  p_lr_b <- mean(lr[1:R_b] >= lr_obs)
  p_cox <- mean(cox >= cox_obs)
  p_cox_a <- mean(cox_a >= cox_obs_a)
  p_cox_b <- mean(cox_b >= cox_obs_b)
  
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
future::plan("multisession", workers = 16)
res<- foreach(
  i=1:N_sim,
  .combine = "rbind",
  .errorhandling = "pass"
) %dorng% {
  simu(n= 300, shape = 2, rate = 1, 
       seed = seed[i], R = 5000, R_a = 250, R_b = 100)
}
toc()

# save results
save(res, file = "simu_result_7.Rdata")
#load("simu_result_7.Rdata")

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
mean(abs(results$Cox_obs)>qnorm(0.975))*100
mean(abs(results$Cox_obs_a)>qnorm(0.975))*100
mean(abs(results$Cox_obs_b)>qnorm(0.975))*100
mean(abs((results$LR_obs-results$hat_mean)/
           sqrt(results$hat_variance))>qnorm(0.975))*100
mean(abs((results$LR_obs-results$hat_mean_a)/
           sqrt(results$hat_variance_a))>qnorm(0.975))*100
mean(abs((results$LR_obs-results$hat_mean_b)/
           sqrt(results$hat_variance_b))>qnorm(0.975))*100

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
