library(gtools) 
library(Rcpp)

## source helper functions written in C++
sourceCpp('MCEM_helper.cpp')


## main function to implement MCMC-EM algorithm for topic model estimation
## without log-likelihood reported in the run
MCEM = function(D, k, num_it = 50, gibbs_steps = 100, thin = 1, tol = 10^-9, alpha = rep(1, k), 
                 log_step = 10, burn_in = 50, verbose=TRUE, C_0 = NULL){
  # D: term-doc matrix (V by d), should be in the form of dgCMatrix
  # k: number of topics 
  # num_it: number of iterations for EM
  # gibbs_step: number of iterations in MCMC in each E-step
  # thin:  keeps every thin sample in E-step 
  # tol: if relative change is below tol, then stop
  # alpha: prior of columns of W
  # log_step: number of EM iterations to print relative change
  # burn_in: disgard the first burn_in samples in MCMC
  # verbose: whether or not to print relative change
  # C_0: initial value of C
  
  d = ncol(D) # number of docs
  V = nrow(D) # number of vocab
  
  #### initialize topic matrix C #####
  if (!is.null(C_0)){
    C_e = C_0
  }else{
    C_e = t(rdirichlet(k, rep(1,V)))
  }
  
  time = Sys.time()

  for (s in (1:num_it)){
    ## initialize mixing matrix W
    W_0 = t(rdirichlet(d, alpha))
    
    #### E step #####
    ## Gibbs sampling for Z and W
    E_step_result = E_step(C = C_e, D = D, W = W_0, alpha = alpha,  
                           gibbs_steps =  gibbs_steps, thin = thin, 
                           burn_in = burn_in)
    C_mean = E_step_result$C_mean
    W_mean = E_step_result$W_mean
    ## estimator for mixing matrix W
    W = W_mean 

    #### M step ####
    ## estimator for topic matrix C
    C_e = apply(C_mean,2, function(x) x/(sum(x))) 
    
    ## print number of iterations, relative change of log-likelihood and running time
    if (verbose && s%%log_step == 0){
      new_time = Sys.time()
      print(paste("iteration:", s, ", time:", round(new_time-time,2)))
      time = new_time
    }
  }

  ## return the estimators for the topic matrix and mixing matrix
  return (list(C = C_e, W = W_mean))
}


