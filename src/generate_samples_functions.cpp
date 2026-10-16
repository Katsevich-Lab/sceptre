#include <boost/random.hpp>
#include <Rcpp.h>
#include <cmath>
#include "sampling_safety.h"
using namespace Rcpp;

// Native sampler streams are local Boost engines, independent of R draws.
// Generated Rcpp RNGScope behavior is unchanged: an absent R seed may be
// initialized at the wrapper boundary, but an existing R stream is not advanced.


// Reset scratch x to 0..n_tot-1, then perform M partial Fisher-Yates swaps.
// The selected ordered sample occupies x[n_tot-M..n_tot); callers preserve
// that extraction order. Scratch/i_doub_array have n_tot/M capacity, 0<=M<=n_tot.
// Advance the supplied Boost stream by exactly M uniform draws; no R RNG/state.
void draw_wor_sample(boost::random::mt19937& generator, boost::random::uniform_real_distribution<double>& distribution, const std::vector<double>& i_doub_array, std::vector<int>& x, int n_tot, int M) {
  for (int i = 0; i < n_tot; i ++) x[i] = i;
  double n_tot_doub = (double) n_tot, u, temp;
  int pos;
  // perform the swap
  for (int i = 0; i < M; i ++) {
    u = distribution(generator);
    pos = floor((n_tot_doub - i_doub_array[i]) * u);
    temp = x[pos];
    x[pos] = x[n_tot - i - 1];
    x[n_tot - i - 1] = temp;
  }
  return;
}


// Own one B-row bank of M distinct zero-based positions in the current fit
// universe [0,n_tot). Workers reuse prefixes/selected stages of this full bank;
// no source-cell IDs or R one-based positions enter it. Empty B/M is supported.
// Each call restarts Boost mt19937(4), preserving swap/extraction order without
// consuming R random draws. Guard sizes, protect/register a finalizing NULL XPtr,
// then allocate native storage; transfer ownership only after the whole bank is
// ready. R errors at the pointer boundary therefore cannot strand a live bank.
// [[Rcpp::export]]
SEXP fisher_yates_sampler(int n_tot, int M, int B) {
  sceptre::sampling::nonnegative(n_tot, "n_tot");
  sceptre::sampling::nonnegative(M, "M");
  if (M > n_tot) stop("M must not exceed n_tot.");
  sceptre::sampling::bank_extent(B, M);
  Rcpp::XPtr<std::vector<std::vector<int>>> pointer(
    static_cast<std::vector<std::vector<int>>*>(nullptr), true);
  // initialize output vector, x vector shared across samples
  std::unique_ptr<std::vector<std::vector<int>>> synth_idx_list(
    new std::vector<std::vector<int>>(B));
  std::vector<int> x(n_tot);

  // initialize the random number generator
  boost::random::mt19937 generator(4);
  boost::random::uniform_real_distribution<double> distribution(0, 1);

  // initialize array of i doubles
  std::vector<double> i_doub_array(M);
  for (int i = 0; i < M; i ++) i_doub_array[i] = (double) i;

  // loop from 0 to B, generating WOR samples
  for (int j = 0; j < B; j ++) {
    std::vector<int> v(M);
    // perform the swap
    draw_wor_sample(generator, distribution, i_doub_array, x, n_tot, M);
    // load the m WOR samples into v
    for (int i = n_tot - M; i < n_tot; i ++) v[i - n_tot + M] = x[i];
    // store v within synth_idx_list
    (*synth_idx_list)[j] = v;
  }
  return sceptre::sampling::register_bank(pointer, synth_idx_list);
}


// Full B-row bank for NT-control fits with N controls and treatment sizes
// m..M (N>=1, 1<=m<=M). Each stored row has M distinct zero-based positions;
// its size-k PREFIX, m<=k<=M, is valid in the N+k fit universe. A suffix can be
// outside a smaller current fit, so consumers must not validate/read the whole
// row as though it belonged to that smaller view. This prefix coupling is why
// one bank serves different guides instead of redrawing separately per guide.
// Preserve the initial Fisher-Yates extraction, inductive update order and
// fresh Boost mt19937(4) stream. No R draws/input mutation; NULL-XPtr-before-bank
// ownership and size/overflow guards are the same as fisher_yates_sampler.
// [[Rcpp::export]]
SEXP hybrid_fisher_iwor_sampler(int N, int m, int M, int B) {
  if (N < 1 || m < 1 || M < m) stop("Hybrid sampling requires N >= 1 and 1 <= m <= M.");
  sceptre::sampling::sum(N, M, "N + M");
  sceptre::sampling::sum(M, 1, "M + 1");
  sceptre::sampling::bank_extent(B, M);
  Rcpp::XPtr<std::vector<std::vector<int>>> pointer(
    static_cast<std::vector<std::vector<int>>*>(nullptr), true);
  // initialize output vector
  std::unique_ptr<std::vector<std::vector<int>>> synth_idx_list(
    new std::vector<std::vector<int>>(B));

  // initialize x vector shared across samples; also initialize i_doub vector
  std::vector<int> x(N + m);
  std::vector<double> i_doub_array(M + 1);
  for (int i = 0; i < M + 1; i ++) i_doub_array[i] = (double) i;

  // initialize the random number generator
  boost::random::mt19937 generator(4);
  boost::random::uniform_real_distribution<double> distribution(0, 1);

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
  return sceptre::sampling::register_bank(pointer, synth_idx_list);
}


namespace {
// Validate finite probabilities in [0,1] and integer-indexed B/cell capacities
// before allocating any native bank. Endpoints, empty cells and B=0 are valid;
// these checks never shrink a request or replace its statistical mechanism.
void check_crt_request(const NumericVector& probabilities, int B) {
  sceptre::sampling::nonnegative(B, "B");
  if (probabilities.size() > std::numeric_limits<int>::max())
    stop("CRT cell count exceeds native integer capacity.");
  for (double probability : probabilities) {
    if (!std::isfinite(probability) || probability < 0 || probability > 1)
      stop("Fitted probabilities must be finite and in [0, 1].");
  }
  sceptre::sampling::bank_extent(B, static_cast<int>(probabilities.size()));
}
}


// Bernoulli CRT bank: B variable-length rows of zero-based fit-cell positions.
// For each cell in input order, draw its number of inclusions with Boost
// binomial(B,p), then a Fisher-Yates subset of draw rows. Appending cells in
// this order preserves every full-bank draw and the original random stream;
// it is not equivalent to a row-major Bernoulli loop for seeded parity.
// Both CRT workers create the complete B1+B2+B3 bank once per guide and reuse
// it across responses; observed n_trt does not set synthetic row cardinality.
// Fresh Boost mt19937(4), no R random draws or input mutation. Protect the finalizing
// NULL XPtr before native allocation and attach ownership only on success.
// [[Rcpp::export]]
SEXP sample_crt_indices(NumericVector fitted_probabilities, int B) {
  check_crt_request(fitted_probabilities, B);
  Rcpp::XPtr<std::vector<std::vector<int>>> pointer(
    static_cast<std::vector<std::vector<int>>*>(nullptr), true);
  // initialize output vector
  std::unique_ptr<std::vector<std::vector<int>>> synth_idx_list(
    new std::vector<std::vector<int>>(B));

  // initialize the random number generator
  boost::random::mt19937 generator(4);
  boost::random::uniform_real_distribution<double> unif_distribution(0, 1);

  // initialize remaining pieces
  std::vector<double> i_doub_array(B);
  for (int i = 0; i < B; i ++) i_doub_array[i] = (double) i;
  std::vector<int> x(B);
  int M;

  for (int j = 0; j < fitted_probabilities.size(); j ++) {
    boost::random::binomial_distribution<int> binom_distribution(B, fitted_probabilities[j]);
    // draw M
    M = binom_distribution(generator);
    // sample WOR of size M from 1 ... B
    draw_wor_sample(generator, unif_distribution, i_doub_array, x, B, M);
    // place j into the sampled positions
    for (int i = B - M; i < B; i ++) (*synth_idx_list)[x[i]].push_back(j);
  }

  return sceptre::sampling::register_bank(pointer, synth_idx_list);
}
