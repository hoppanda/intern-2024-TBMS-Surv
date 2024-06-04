library(flexsurv)

### derive p-value distribution

fit<-flexsurvspline(formula=Surv(futime,fustat)~1,data=ovarian,k=2) 
fit_a<-flexsurvspline(formula=Surv(futime,fustat)~1,data=ovarian,k=3) 
plot(fit)
nd <- data.frame(id=1:100)

simdat <- simulate(fit,nsim=1000,seed=12345, censtime = max(ovarian$futime), newdata = nd) 
sim_lrs <- numeric(1000)
sim_lrs_fixed <- numeric(1000)
sim_msg <- character(1000)

for(i in 1:1000){
  if(i%%50==0) print(i)
  skip_ind <- FALSE
  usedat <- data.frame(time=simdat[[paste0("time_",i)]],
                       event=simdat[[paste0("event_",i)]])
  tryCatch({
    mod1_default <- flexsurvspline(formula=Surv(time,event)~1,data=usedat,k=2) 
    mod1_fixed_2knot <- flexsurvspline(formula=Surv(time,event)~1,data=usedat,knots = fit$knots[2:3]) 
    mod2_default <- flexsurvspline(formula=Surv(time,event)~1,data=usedat,k=3) 
    mod2_fixed_3knot <- flexsurvspline(formula=Surv(time,event)~1,data=usedat,knots = fit_a$knots[2:4]) 
  },
  error = function(e) {
    skip_ind <<- TRUE
    sim_lrs[i] <<- NA
    sim_lrs_fixed[i] <<- NA
    sim_msg[i] <<- "error"
  })

  if(skip_ind) next
  sim_lrs[i] <- -2*(mod1_default$loglik - mod2_default$loglik) 
  sim_lrs_fixed[i] <- -2*(mod1_fixed_2knot$loglik - mod2_fixed_3knot$loglik) 
}

sum(is.na(sim_lrs))
sum(is.na(sim_lrs_fixed))

sim_pv <- 1-pchisq(na.omit(sim_lrs),df=1)
sim_pv_fixed <- 1-pchisq(na.omit(sim_lrs_fixed),df=1)
 

par(mfrow=c(1,1),mgp=c(1.5,0.5,0),cex.lab = 1, cex.main=0.9, tcl=-0.13)
hist(sim_pv,freq=FALSE,col="dodgerblue",border="white",
     main="Distribution of p-values of LRT via simulation\n (algorithm-determined knots per simulated dataset)",
     xlab="p-value",ylab="Density")
hist(sim_pv_fixed,freq=FALSE,col="dodgerblue",border="white",
     main="Distribution of p-values of LRT via simulation\n (pre-defined knots)")

hist(sim_lrs,freq=FALSE,col="dodgerblue",border="white",
     main="Distribution of Likelihood Ratios Statistics via simulation\n (algorithm-determined knots per simulated dataset)",
     xlab="-2*(LogLik_Null-LogLik_Alternative)",ylab="Density")
lines(x=seq(0.1,15,length=500),y=seq(0.1,15,length=500),col="red")

sim_lrs_fixed2 <- numeric(1000)
sim_msg2 <- character(1000)

for(i in 1:1000){
  if(i%%50==0) print(i)
  skip_ind <- FALSE
  usedat <- data.frame(time=simdat[[paste0("time_",i)]],
                       event=simdat[[paste0("event_",i)]])
  tryCatch({
    mod1_fixed_2knot_new <- flexsurvspline(formula=Surv(time,event)~1,data=usedat,knots = fit$knots[2:3]) 
    mod2_fixed_3knot_new <- flexsurvspline(formula=Surv(time,event)~1,data=usedat,knots = c(fit$knots[2:3], 6.15)) 
  },
  error = function(e) {
    skip_ind <<- TRUE
    sim_lrs_fixed2[i] <<- NA
    sim_msg2[i] <<- "error"
  })
  
  if(skip_ind) next
  sim_lrs_fixed2[i] <- -2*(mod1_fixed_2knot_new$loglik - mod2_fixed_3knot_new$loglik) 
}

sim_pv_fixed2 <- 1-pchisq(na.omit(sim_lrs_fixed2),df=1)

hist(sim_pv_fixed2,freq=FALSE,col="dodgerblue",border="white",
     main="Distribution of p-values of LRT via simulation\n (pre-defined knots)",
     xlab="p-value", ylab="Density")
