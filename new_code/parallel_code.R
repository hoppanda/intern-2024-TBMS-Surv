rm(list=ls())
# library packages
library(parallel)
library(doParallel)
library(foreach)

detectCores()

# no_cores <- detectCores(logical = TRUE) 
# cl <- makeCluster(no_cores-1, type='PSOCK')  
my.cluster <- parallel::makeCluster(
  16, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()


# foreach parallel - 1 secs
time_res <- system.time(
  foreach(i = 1:24, .combine=cbind) %dopar% {
    Sys.sleep(1)
  }
)
stopCluster(my.cluster)

test <- data.frame(time = time_res[3])
test

#saveRDS(test, file = "/SFS/user/ctc/farinare/code")

# default without using parallel - 10 secs
system.time(
  for(i in 1:23){
    Sys.sleep(1)
  }
)
