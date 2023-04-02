#include <random>
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' @title Fisher-Yates sampler
//' @description This function draws a without replacement sample using the Fisher-Yates sampling algorithm
//' @param n_tot the total number of cells
//' @param M the maximum number of cells in a given gRNA group
//' @param B the number of WOR samples to generate
//' @noRd
// [[Rcpp::export]]
SEXP fisher_yates_samlper(int n_tot, int M, int B) {
  // initialize output vector, x vector shared across samples
  std::vector<std::vector<int>>* synth_idx_list = new std::vector<std::vector<int>>(B);
  std::vector<int> x(n_tot);
  
  // initialize the random number generator
  std::mt19937 generator(4);
  std::uniform_real_distribution<double> distribution(0, 1);
  
  // initialize array of i doubles
  std::vector<double> i_doub_array(M);
  for (int i = 0; i < M; i ++) i_doub_array[i] = (double) i;
  
  // initialize remaining pieces
  double n_tot_doub = (double) n_tot, mult_fact, u, temp;
  int pos;

  // loop from 0 to B, generating WOR samples
  for (int j = 0; j < B; j ++) {
    // fill x, initialize v
    for (int i = 0; i < n_tot; i ++) x[i] = i;
    std::vector<int> v(M);
    
    // perform the swap
    for (int i = 0; i < M; i ++) {
      u = distribution(generator);
      pos = floor((n_tot_doub - i_doub_array[i]) * u);
      temp = x[pos];
      x[pos] = x[n_tot - i - 1];
      x[n_tot - i - 1] = temp;
    }
    
    // load the m WOR samples into v
    for (int i = n_tot - M; i < n_tot; i ++) {
      v[i - n_tot + M] = x[i];
    }
    // store v within synth_idx_list
    (*synth_idx_list)[j] = v;
  }
  Rcpp::XPtr<std::vector<std::vector<int>>> ptr(synth_idx_list);
  return ptr;
}


//' @title Hybrid Fisher-Yates/IWOR sampler
//' @description This function draws an inductive without replacement sample using the hybrid Fisher-Yates/IWOR sampling algorithm (developed by Barry et al, to be described in a forthcoming preprint)
//' @param N the number of control cells
//' @param m the minumum number of cells in a given gRNA group
//' @param M the maximum number of cells in a given gRNA group
//' @param B the number of WOR samples to generate
//' @noRd
// [[Rcpp::export]]
SEXP hybrid_fisher_iwor_sampler(int N, int m, int M, int B) {
  // initialize output vector
  std::vector<std::vector<int>>* synth_idx_list = new std::vector<std::vector<int>>(B);

  // initialize x vector shared across samples; also initialize i_doub vector
  std::vector<int> x(N + m);
  std::vector<double> i_doub_array(M);
  for (int i = 0; i <= M; i ++) i_doub_array[i] = (double) i;
  
  // initialize the random number generator
  std::mt19937 generator(4);
  std::uniform_real_distribution<double> distribution(0, 1);
  
  // initialize required variables
  double N_doub = (double) N, m_doub = (double) m, i_doub, p, u;
  int pos, temp;
  
  // loop over B
  for (int j = 0; j < B; j ++) {
    // fill x with 0...N+m-1
    for (int i = 0; i < N + m; i ++) x[i] = i;
    // initialize a new output vector v of size M
    std::vector<int> v(M);
    
    // perform a Fisher-Yates shuffle on x, sampling m elements
    for (int i = 0; i < m; i ++) {
      u = distribution(generator);
      pos = floor((N_doub + m_doub - i_doub_array[i]) * u);
      temp = x[pos];
      x[pos] = x[N + m - i - 1];
      x[N + m - i - 1] = temp; 
    }
    
    // load the m WOR samples into v
    for (int i = 0; i < m; i ++) {
      v[i] = x[i+N];
    }
    
    // perform the inductive WOR step
    x[N] = N + m;
    for (int i = m + 1; i <= M; i++) {
      // sample from iwor distribution
      u = distribution(generator);
      i_doub = i_doub_array[i];
      p = i_doub/(N_doub + i_doub);
      if (u > 1 - p) {
        pos = N;
      } else {
        pos = floor(u * N_doub/(1 - p));
      }
      
      // perform the swap and extraction
      v[i - 1] = x[pos];
      x[pos] = x[N];
      x[N] = N + i;
    }
    
    // store v within synth_idx_list
    (*synth_idx_list)[j] = v;
  }
  // return a pointer to the synthetic idx list
  Rcpp::XPtr<std::vector<std::vector<int>>> ptr(synth_idx_list);
  return ptr;
}
