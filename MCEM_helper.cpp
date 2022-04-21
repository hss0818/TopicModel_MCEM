#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


//' draw sample from categorical distribution 
//' using built-in uniform distribution sampler
//' by comparing the uniform sample with cumulative sum of p
//' 
//' @param P: categorical distribution parameter
//' @return a sample from categorical distribution
// [[Rcpp::export]]
unsigned int rmultinomial(const arma::vec& P){
  double u = R::runif(0, sum(P));
  double tmp = 0;
  for (int i = 0; i < P.size(); ++i){
    tmp += P(i);
    if (u < tmp){
      return i;
    }
  }
  return P.size() - 1;
}


//' draw sample from multinomial distribution
//' using categorical distribution sampler
//'
//' @param n, ps: multinomial distribution parameter
//' @return a sample from multinomial distribution
// [[Rcpp::export]]
IntegerVector MultinomCalt(unsigned int n, arma::vec & ps) {
  IntegerVector res(ps.size());
  for (unsigned int i = 0; i < n; i++){
    unsigned int index = rmultinomial(ps);
    res(index) += 1;
  }
  return res;
}


//' draw sample from Dirichlet distribution
//' using built-in Gamma distribution sampler
//'
//' @param deltas: Dirichlet distribution parameter
//' @return a sample from Dirichlet distribution
// [[Rcpp::export]]
arma::vec rDirichlet(const arma::vec& deltas){
  unsigned int C = deltas.n_elem;
  arma::vec Xgamma(C);
  
  //generating gamma(deltac,1)
  for(unsigned int c=0;c<C;c++){
    Xgamma(c) = R::rgamma(deltas(c),1.0);
  }
  return Xgamma/sum(Xgamma);
}


//' implement one iteration of Gibbs sampling within the E-step of the MCMC-EM algorithm
//' 
//' @param D: term-doc matrix (V by d)
//' @param C: current estimator of the topic matrix (V by k)
//' @param W: Gibbs sample of W|Z from last iteration
//' @param alpha: prior of columns of mixing matrix W
//' @return C_hat: intermediate variable calculated from Gibbs sample of Z|X,C, 
//' will be used for the estimation of the topic matrix
//' @return W_hat: Gibbs sample of W|Z
// [[Rcpp::export]]
Rcpp::List one_it(const arma::sp_mat & D, const arma::mat& C, const arma::mat W, const arma::vec & alpha){
  unsigned int V = D.n_rows; // number of vocab
  unsigned int d = D.n_cols; // number of docs
  unsigned int k = C.n_cols; // number of topics
  
  // initialize matrices and iterators
  arma::mat C_hat = arma::zeros<arma::mat>(V, k);
  arma::mat W_beta = arma::zeros<arma::mat>(k, d);
  arma::mat W_hat = arma::zeros<arma::mat>(k, d);
  unsigned int i;
  unsigned int j;
  unsigned int x;
  
  // iterate over all elements of D
  for(arma::sp_mat::const_iterator it = D.begin(); it != D.end(); ++it){
    i = it.col(); // Col position: document
    j = it.row(); // Row position: word index
    x =  *it;    // Value: word frequent
    arma::vec p = C.row(j).t()% W.col(i); // element-wise multiplication, multinomial parameter
    
    // sample Z|X,C
    if (x > 4){
      IntegerVector a = MultinomCalt(x, p);
      for (unsigned int l = 0; l < k; l++)
      {
        C_hat(j, l) = C_hat(j,l) + a(l);
        W_beta(l, i) = W_beta(l,i) + a(l);
      }
    }
    else{
      for (unsigned int l = 0; l < x; l++)
      {
        int a = rmultinomial(p) ;
        C_hat(j, a) = C_hat(j,a) + 1;
        W_beta(a, i) = W_beta(a,i) + 1;
      }
    }
  }
  
  // sample W|Z
  for (unsigned f = 0; f < d; f++){
    arma::vec p_w = W_beta.col(f) + alpha; // Dirichlet parameter
    W_hat.col(f) = rDirichlet(p_w);
  }
  
  return Rcpp::List::create(Rcpp::Named("C_hat")=C_hat, 
                            Rcpp::Named("W_hat")=W_hat);
}


//' implement the E-step of the MCMC-EM Algorithm (Algorithm 1 in supplement material)
//' 
//' @param C: current estimator of the topic matrix (V by k)
//' @param D: term-doc matrix (V by d)
//' @param W: randomly initialized mixing matrix W as the start point for Gibbs sampling
//' @param alpha: prior of columns of mixing matrix W
//' @param gibbs_steps: number of iterations in Gibbs sampling
//' @param thin: keeps every thin sample in Gibbs sampling
//' @return C_mean: raw estimator of the topic matrix C before column normalization
//' @return W_mean: estimator of the mixing matrix W calculated from Gibbs samples of Z|W
// [[Rcpp::export]]
Rcpp::List E_step(arma::mat & C, const arma::sp_mat & D, arma::mat W,
                  const arma::vec & alpha, unsigned int gibbs_steps, unsigned int thin,
                  double burn_in){
  unsigned int V = D.n_rows; // number of vocab
  unsigned int d = D.n_cols; // number of docs
  unsigned int k = C.n_cols; // number of topics
  double count = 0.0; // count of number of Gibbs samples
  
  // initialize cumulative sum for matrix from Gibbs sampling
  arma::mat C_cul = arma::zeros<arma::mat>(V, k);
  arma::mat W_cul = arma::zeros<arma::mat>(k, d);
  
  // Gibbs sampling
  for (unsigned int t = 0; t < gibbs_steps; t++){
    List L = one_it(D, C, W, alpha); // call function to implement one iteration
    arma::mat C_hat = L["C_hat"];
    arma::mat W_hat = L["W_hat"];
    W = W_hat;
    
    // cumulate the sample for every thin sample after burn-in
    if (((t - 1)%thin ==0) && (t > burn_in)){
      count = count + 1.0;
      C_cul = C_cul + C_hat;
      W_cul = W_cul + W_hat;
    }
  }
  
  // average the cumulative sum
  arma::mat C_mean = C_cul/count;
  arma::mat W_mean = W_cul/count;
  
  return  Rcpp::List::create(Rcpp::Named("C_mean")=C_mean, 
                             Rcpp::Named("W_mean")=W_mean);
}


//' calculate the log-likelihood 
//' 
//' @param D: term-doc matrix (V by d)
//' @param C: estimator of the topic matrix (V by k)
//' @param W: estimator of the mixing matrix (k by d)
//' @param thr: to avoid numeric issue when take log
//' @return lkh: the log-likelihood 
// [[Rcpp::export]]
double lkh_cal(const arma::sp_mat & D, const arma::mat& C, const arma::mat& W, const double thr=1e-9){
  unsigned int i;
  unsigned int j;
  unsigned int x;
  double lkh;
  lkh = 0.0;
  
  // iterate over all elements of D
  for(arma::sp_mat::const_iterator it = D.begin(); it != D.end(); ++it){
    i = it.row(); // Row position: word index
    j = it.col(); // Col position: document
    x =  *it;    // Value: word frequency
    double u = dot(C.row(i),W.col(j));
    if (u < thr){
      u = thr;
    }
    lkh = lkh + x*log(u);
  }
  
  return lkh;
}


