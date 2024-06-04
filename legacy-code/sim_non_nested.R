library(flexsurv)

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


fit0<-flexsurvreg(formula=Surv(recyrs,censrec)~1,data=bc[bc$group=="Poor",],dist="exp") 
fit1<-flexsurvreg(formula=Surv(recyrs,censrec)~1,data=bc[bc$group=="Poor",],dist="weibull") 
fit2<-flexsurvreg(formula=Surv(recyrs,censrec)~1,data=bc[bc$group=="Poor",],dist="lnorm") 
fit3<-flexsurvreg(formula=Surv(recyrs,censrec)~1,data=bc[bc$group=="Poor",],dist="llogis") 
fit4<-flexsurvreg(formula=Surv(recyrs,censrec)~1,data=bc[bc$group=="Poor",],dist="gompertz") 
fit5<-flexsurvreg(formula=Surv(recyrs,censrec)~1,data=bc[bc$group=="Poor",],dist="gamma") 
fit6<-flexsurvreg(formula=Surv(recyrs,censrec)~1,data=bc[bc$group=="Poor",],dist="gengamma") 
fit7<-flexsurvreg(formula=Surv(recyrs,censrec)~1,data=bc[bc$group=="Poor",],dist="genf") 

fit_rcs2 <-flexsurvspline(formula=Surv(recyrs,censrec)~1,data=bc[bc$group=="Poor",], k=2) 
fit_rcs3 <-flexsurvspline(formula=Surv(recyrs,censrec)~1,data=bc[bc$group=="Poor",], k=3) 

gof <-rbind(fitstats.flexsurvreg(fit0),
            fitstats.flexsurvreg(fit1),
            fitstats.flexsurvreg(fit2),
            fitstats.flexsurvreg(fit3),
            fitstats.flexsurvreg(fit4),
            fitstats.flexsurvreg(fit5),
            fitstats.flexsurvreg(fit6),
            fitstats.flexsurvreg(fit7))
gof$mdoel <- c("exp","weibull","lnorm","llogis","gompertz","gamma","gengamma","genf")
gof


par(mfrow=c(2,2),tcl=-0.13,mgp=c(1.5,0.5,0),cex.main=1,cex.lab=1,mar=c(4,2,2,1))
plot(fit2,main="Log-Normal")
plot(fit3,main="Log-Logistics")
plot(fit6,main="Generalized-Gamma")
plot(fit7,main="Generalized-F")

par(mfrow=c(1,1))
plot(fit0,main="Exponential")

sim1 <- simulate(fit2, nsim=1000, seed=12345, censtime = max(bc$recyrs[bc$group=="Poor"]))

lr_lnorm_llogis <- numeric(1000)
lr_lnorm_gengamma <- numeric(1000)
lr_lnorm_genf <- numeric(1000)

for(i in 1:1000){
  if(i%%50==0) print(i)
  skip_ind <- FALSE
  usedat <- data.frame(time=sim1[[paste0("time_",i)]],
                       event=sim1[[paste0("event_",i)]])
  
  tryCatch({
    tmpfit_lnorm <- flexsurvreg(formula=Surv(time,event)~1,data=usedat,dist="lnorm") 
    tmpfit_llogis <- flexsurvreg(formula=Surv(time,event)~1,data=usedat,dist="llogis") 
    tmpfit_gengamma <- flexsurvreg(formula=Surv(time,event)~1,data=usedat,dist="gengamma") 
    tmpfit_genf <- flexsurvreg(formula=Surv(time,event)~1,data=usedat,dist="genf") 
  },
  error = function(e) {
    skip_ind <<- TRUE
    lr_lnorm_llogis[i] <<- NA
    lr_lnorm_gengamma[i] <<- NA
    lr_lnorm_genf[i] <<- NA
  })
  
  if(skip_ind) next
  lr_lnorm_llogis[i] <- -2*(tmpfit_lnorm$loglik - tmpfit_llogis$loglik) 
  lr_lnorm_gengamma[i] <- -2*(tmpfit_lnorm$loglik - tmpfit_gengamma$loglik) 
  lr_lnorm_genf[i] <- -2*(tmpfit_lnorm$loglik - tmpfit_genf$loglik) 
}


pv_lnorm_llogis <- sum(lr_lnorm_llogis>-2*(fit2$loglik - fit3$loglik),na.rm=TRUE)/(sum(!is.na((lr_lnorm_llogis))))
pv_lnorm_gengamma <- sum(lr_lnorm_gengamma>-2*(fit2$loglik - fit6$loglik),na.rm=TRUE)/(sum(!is.na((lr_lnorm_gengamma))))
pv_lnorm_genf <- sum(lr_lnorm_genf>-2*(fit2$loglik - fit7$loglik),na.rm=TRUE)/(sum(!is.na((lr_lnorm_genf))))

par(mfrow=c(1,3),tcl=-0.13,mgp=c(1.5,0.5,0),cex.main=1,cex.lab=1)
hist(lr_lnorm_llogis, xlab="-2*(LogLik(Null) - LogLik(Alternative))", ylab="",main="Histrogram of Likelihood Ratio Statistics between Models\n with Log-Normal (Null) and Log-Logistic (Alternative)\n via a parametric Bootstrap procedure",
     col="lightgrey",border="white")
abline(v=-2*(fit2$loglik - fit3$loglik),lty=4, col="deeppink")
legend("topright",bty="n", legend=paste("p-value:",round(pv_lnorm_llogis,3)),cex=1.1,text.col="deeppink",text.font=7)
hist(lr_lnorm_gengamma, xlab="-2*(LogLik(Null) - LogLik(Alternative))", ylab="",main="Histrogram of Likelihood Ratio Statistics between Models\n with Log-Normal (Null) and Generalized-Gamma (Alternative)\n via a parametric Bootstrap procedure",
     col="lightgrey",border="white")
abline(v=-2*(fit2$loglik - fit6$loglik),lty=4, col="deeppink")
legend("topright",bty="n", legend=paste("p-value:",round(pv_lnorm_gengamma,3)),cex=1.1,text.col="deeppink",text.font=7)
hist(lr_lnorm_genf, xlab="-2*(LogLik(Null) - LogLik(Alternative))", ylab="",main="Histrogram of Likelihood Ratio Statistics between Models\n with Log-Normal (Null) and Generalized-F (Alternative)\n via a parametric Bootstrap procedure",
     col="lightgrey",border="white")
abline(v=-2*(fit2$loglik - fit7$loglik),lty=4, col="deeppink")
legend("topright",bty="n", legend=paste("p-value:",round(pv_lnorm_genf,3)),cex=1.1,text.col="deeppink",text.font=7)



sim_exp <- simulate(fit0, nsim=1000, seed=12345, censtime = max(bc$recyrs[bc$group=="Poor"]))
lr_exp_lnorm <- numeric(1000)

for(i in 1:1000){
  if(i%%50==0) print(i)
  skip_ind <- FALSE
  usedat <- data.frame(time=sim_exp[[paste0("time_",i)]],
                       event=sim_exp[[paste0("event_",i)]])

  tryCatch({
    tmpfit_exp <- flexsurvreg(formula=Surv(time,event)~1,data=usedat,dist="exp")
    tmpfit_lnorm <- flexsurvreg(formula=Surv(time,event)~1,data=usedat,dist="lnorm")
  },
  error = function(e) {
    skip_ind <<- TRUE
    lr_exp_lnorm[i] <<- NA
  })

  if(skip_ind) next
  lr_exp_lnorm[i] <- -2*(tmpfit_exp$loglik - tmpfit_lnorm$loglik)
}


sim_lnorm <- simulate(fit2, nsim=1000, seed=12345, censtime = max(bc$recyrs[bc$group=="Poor"]))
lr_lnorm_rcs2 <- numeric(1000)

for(i in 1:1000){
  if(i%%50==0) print(i)
  skip_ind <- FALSE
  usedat <- data.frame(time=sim_lnorm[[paste0("time_",i)]],
                       event=sim_lnorm[[paste0("event_",i)]])
  
  tryCatch({
    tmpfit_lnorm <- flexsurvreg(formula=Surv(time,event)~1,data=usedat,dist="lnorm")
    tmpfit_rcs2 <- flexsurvspline(formula=Surv(time,event)~1,data=usedat,k=2)
  },
  error = function(e) {
    skip_ind <<- TRUE
    lr_lnorm_rcs2[i] <<- NA
  })
  
  if(skip_ind) next
  lr_lnorm_rcs2[i] <- -2*(tmpfit_lnorm$loglik - tmpfit_llogis$loglik)
}


sim_rcs2 <- simulate(fit_rcs2, nsim=1000, seed=12345, censtime = max(bc$recyrs[bc$group=="Poor"]))
lr_rcs2_rcs3 <- numeric(1000)

for(i in 1:1000){
  if(i%%50==0) print(i)
  skip_ind <- FALSE
  usedat <- data.frame(time=sim_rcs2[[paste0("time_",i)]],
                       event=sim_rcs2[[paste0("event_",i)]])
  
  tryCatch({
    tmpfit_rcs2 <- flexsurvspline(formula=Surv(time,event)~1,data=usedat,k=2)
    tmpfit_rcs3 <- flexsurvspline(formula=Surv(time,event)~1,data=usedat,k=3)
  },
  error = function(e) {
    skip_ind <<- TRUE
    lr_rcs2_rcs3[i] <<- NA
  })
  
  if(skip_ind) next
  lr_rcs2_rcs3[i] <- -2*(tmpfit_rcs2$loglik - tmpfit_rcs3$loglik)
}

pv_exp_lnorm <- sum(lr_exp_lnorm>-2*(fit0$loglik - fit2$loglik),na.rm=TRUE)/(sum(!is.na((lr_exp_lnorm))))
pv_lnorm_rcs2 <- sum(lr_lnorm_rcs2>-2*(fit2$loglik - fit_rcs2$loglik),na.rm=TRUE)/(sum(!is.na((lr_lnorm_rcs2))))
pv_rcs2_rcs3 <- sum(lr_rcs2_rcs3>-2*(fit_rcs2$loglik - fit_rcs3$loglik),na.rm=TRUE)/(sum(!is.na((lr_rcs2_rcs3))))




