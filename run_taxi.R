library(Matrix)

set.seed(10)
source("MCEM_k.R")


## read data
taxi = read.csv("Data/taxi_19.csv")

## pick and dropoff locations
pickup_name = as.numeric(sapply(2:ncol(taxi), function(x)strsplit(colnames(taxi)[x], "X")[[1]][2]))
dropoff_name = taxi[,1]

D = as(as.matrix(taxi[-which(dropoff_name==264|dropoff_name==265),
                       -c(1,which(pickup_name==264|pickup_name==265))]), "dgCMatrix")

## number of meta-states
k = 9

MCEM_9 = MCEM(D, k = k, num_it = 100, gibbs_steps = 50, log_step = 1, burn_in = 0)

## estimators
em_C = MCEM_9$C
em_W = MCEM_9$W

## estimated disaggregation and aggregation distributions
C_df = data.frame(location_id = dropoff_name[-which(dropoff_name==264|dropoff_name==265)], C = em_C)
W_df = data.frame(location_id = pickup_name[-which(pickup_name==264|pickup_name==265)], W = t(em_W))


