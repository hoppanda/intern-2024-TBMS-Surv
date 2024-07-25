library(muhaz)
library(bshazard)
library(flexsurvreg)

data <- read.csv("DANUBE_0.csv")
data <- read.csv("DANUBE_1.csv")
data <- read.csv("IMP150_0.csv")
data <- read.csv("IMP150_1.csv")


#muhaz
m=muhaz(data$time,data$event,min.time = min(data$time),max.time = max(data$time))
plot(m)

#bshazard
haz=bshazard(Surv(time,event) ~ 1,data=data)
plot(haz, xlab = "Time", main="Treatment estimated hazard rate", conf.int = FALSE)


km1 <-survfit(Surv(time,event) ~ 1, data=data) 
plot(km1)
