library(Matrix)
library(combinat)
library(here)

here::i_am("example.R")
source(here("MCEM.R"))

set.seed(0)

V = 1200 # number of vocab
d = 1000 # number of docs
k = 5 # number of topics 
n = 1000 # number of words in each doc

## generate vocab-doc matrix
C = t(rdirichlet(k, rep(0.1, V)))
W = t(rdirichlet(d, rep(0.1, k)))
U = C%*%W
D = apply(U, 2, FUN = function(x) rmultinom(1, n, x))
sparse_D = as(D, "sparseMatrix") # to sparse matrix
 
## implement MCMC-EM algorithm
time = Sys.time()
MMSB_est = MCEM(sparse_D, k, num_it = 5, gibbs_steps = 20, burn_in = 0, log_step = 1)
new_time = Sys.time()
new_time - time

## estimators for topic matrix and mixing matrix
MMSB_est$C
MMSB_est$W
