library(Matrix)
library(combinat)
library(doMC)

source("MCEM_k.R")

## load the kos or nips dataset
load("Data/sparse_kos.RData")
# load("Data/sparse_nips.RData")


set.seed(0)

## number of random initial values for the MCMC-EM algorithm
nstart = 12

## list of number of topics
k_list = c(5,7,10,15,20)

## specify number of topics for the job
my_index = 1
k = k_list[my_index]


## function for parallel programming on CPUs
registerDoMC(nstart)
em_nstart <- function(j){
  lda_varname = paste0("MCEM_k_", k, "_j_", j)
  MCEM_est = MCEM(sparse_D, k = k, num_it = 200, gibbs_steps = 100, log_step = 5, burn_in = 0)
  MCEM_est$lkh = lkh_cal(sparse_D, MCEM_est$C, MCEM_est$W)
  my_result =  MCEM_est
  assign(lda_varname, my_result)
  save(list = lda_varname, file = paste0(lda_varname,".RData"))
}


## run in parallel
cat("\n\n\n")
print(paste0("$$$$$$$$$$$$ k = ", k, "  $$$$$$$$$$$$$$"))
time = Sys.time()
foreach(j = c(1:nstart))%dopar%{
  em_nstart(j = j)
}
new_time = Sys.time()
print(paste("time:", round(new_time - time, 2)))

