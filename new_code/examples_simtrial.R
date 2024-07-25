library(simtrial)
library(survival)
library(muhaz)
library(bshazard)

example <- list(ex1_delayed_effect, ex2_delayed_effect, ex3_cure_with_ph,
                ex4_belly, ex5_widening, ex6_crossing)
par(mfrow=c(6,2),mar=c(2, 2, 1, 1))
for(j in 1:6){
  for(i in 0:1){
    data <- example[[j]]
    time <- data$month[data$trt==i]
    event <- data$evntd[data$trt==i]
    fit=muhaz(time, event, min.time = min(time),max.time = max(time))
    plot(fit, main = paste0("Example ",j,", Arm ", i))
  }
}






# EX1
data(ex1_delayed_effect) 
km1 <- with(ex1_delayed_effect,survfit(Surv(month, evntd) ~ trt))
km1 
plot(km1) 
with(subset(ex1_delayed_effect, trt == 1), survfit(Surv(month, evntd) ~ trt)) 
with(subset(ex1_delayed_effect, trt == 0), survfit(Surv(month, evntd) ~ trt))


fit=muhaz(ex1_delayed_effect$month[ex1_delayed_effect$trt==1], 
          ex1_delayed_effect$evntd[ex1_delayed_effect$trt==1],
          min.time = min(ex1_delayed_effect$month[ex1_delayed_effect$trt==1]),
          max.time = max(ex1_delayed_effect$month[ex1_delayed_effect$trt==1]))
plot(fit)


# EX2
data(ex2_delayed_effect)
km1 <- with(ex2_delayed_effect, survfit(Surv(month, evntd) ~ trt)) 
km1 
plot(km1) 
with(subset(ex2_delayed_effect, trt == 1), survfit(Surv(month, evntd) ~ trt)) 
with(subset(ex2_delayed_effect, trt == 0), survfit(Surv(month, evntd) ~ trt))


fit=muhaz(ex2_delayed_effect$month[ex2_delayed_effect$trt==1], 
          ex2_delayed_effect$evntd[ex2_delayed_effect$trt==1])
plot(fit)

# EX3
data(ex3_cure_with_ph)
km1 <- with(ex3_cure_with_ph, survfit(Surv(month, evntd) ~ trt)) 
km1 
plot(km1)
fit=muhaz(ex3_cure_with_ph$month, ex3_cure_with_ph$evntd)
plot(fit)


data(ex4_belly) 
#KM both arms
km1 <- with(ex4_belly, survfit(Surv(month, evntd) ~ trt)) 
km1 
plot(km1)
#KM hazard
km1 <-survfit(Surv(month, evntd) ~ 1, data=ex4_belly[ex4_belly$trt==1,]) 
plot(km1)
time <- km1$time
survival <- km1$surv
cumulative_hazard <- -log(survival)
hazard <- diff(cumulative_hazard) / diff(time)
plot(time[-1], hazard, type="l")
# muhaz
time=ex4_belly$month[ex4_belly$trt==1]
event=ex4_belly$evntd[ex4_belly$trt==1]
fit=muhaz(time, ex4_belly$evntd[ex4_belly$trt==1],
          min.time = min(time),max.time = max(time))
plot(fit, xlab = "Time (months)", main=" Estimated hazard rate")
#bshazard
haz=bshazard(Surv(month, evntd) ~ 1,data=ex4_belly[ex4_belly$trt==1,])
plot(haz, xlab = "Time (months)", main=" Estimated hazard rate", conf.int = FALSE)
#flexsurvreg hazard
data=ex4_belly[ex4_belly$trt==1,]
data=data.frame(time=data$month, event=data$evntd)
fit10 <- flexsurvspline(formula=Surv(time,event)~1,data=data, k=3)
sum_surv <- as.data.frame(summary(fit10))
time <- sum_surv$time
survival <- sum_surv$est
cumulative_hazard <- -log(survival)
hazard <- diff(cumulative_hazard) / diff(time)
plot(time[-1], hazard, type="l")




data(ex5_widening) 
km1 <- with(ex5_widening, survfit(Surv(month, evntd) ~ trt)) 
km1
plot(km1)
fit=muhaz(ex5_widening$month, ex5_widening$evntd)
plot(fit)



data(ex6_crossing) 
km1 <- with(ex6_crossing, survfit(Surv(month, evntd) ~ trt)) 
km1 
plot(km1)
fit=muhaz(ex6_crossing$month[ex6_crossing$trt==1],
          ex6_crossing$evntd[ex6_crossing$trt==1])
plot(fit)

