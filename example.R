library(Matrix)
library(combinat)
# setwd("~/Desktop/TopicModel_code/code_submission")
source("MCEM.R")
set.seed(0)


V = 1200
d = 1000
k = 5
n = 1000

C = t(rdirichlet(k, rep(0.1, V)))
W = t(rdirichlet(d, rep(0.1, k)))
U = C%*%W
D = apply(U, 2, FUN = function(x) rmultinom(1, n, x))
sparse_D = as(D, "sparseMatrix") # to sparse matrix
 

time = Sys.time()
MMSB_est = MCEM(sparse_D, k, num_it = 5, gibbs_steps = 20, burn_in = 0, log_step = 1)
new_time = Sys.time()
new_time - time


MMSB_est$C
MMSB_est$W
 
 
 
 
 
 

