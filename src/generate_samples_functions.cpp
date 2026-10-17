#include <boost/random.hpp>
#include <Rcpp.h>
#include <cmath>
#include "sampling_safety.h"
using namespace Rcpp;

// These samplers use their own Boost random generators, separate from R's.
// The generated Rcpp wrapper can initialize a missing R seed when entered, but
// these C++ draws do not advance an R random stream that already exists.


// Reset the working array x to cell positions 0 through n_tot-1 and perform M
// partial Fisher-Yates swaps (sampling without replacement). Selected positions
// end up in x[n_tot-M, ..., n_tot-1]; keep that extraction order. The caller
// provides space for n_tot working entries and M output entries, with 0 <= M <=
// n_tot. Advance the supplied Boost generator by exactly M uniform draws.
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


// Generate B draws of M distinct cell positions, starting at 0 among n_tot cells.
// These are positions in the selected fitting data, not original cell IDs.
// Workers can reuse the first k positions of each draw and different stage ranges.
// B or M may be zero. Each call starts Boost mt19937 with seed 4; preserve the
// swap and extraction order to keep the exact random draws.
// Create and protect an empty R external pointer, including its cleanup handler,
// before allocating C++ draws. Attach the finished draws only afterward. An R
// error while creating the pointer can skip C++ destructors, so reversing this
// order could leak all the draws.
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

  // Generate B samples without replacement.
  for (int j = 0; j < B; j ++) {
    std::vector<int> v(M);
    // perform the swap
    draw_wor_sample(generator, distribution, i_doub_array, x, n_tot, M);
    // Copy the selected positions into v in extraction order.
    for (int i = n_tot - M; i < n_tot; i ++) v[i - n_tot + M] = x[i];
    // store v within synth_idx_list
    (*synth_idx_list)[j] = v;
  }
  return sceptre::sampling::register_bank(pointer, synth_idx_list);
}


// Generate B shared draws for NT discovery with N controls and treatment groups
// ranging from m to M cells (N >= 1, 1 <= m <= M). Each draw has M distinct
// positions. For a group of k treated cells, use only the first k entries: those
// positions are valid among its N+k fitting cells. Later entries may exceed
// that smaller fit's size. This construction lets differently sized guide groups
// share draws without generating them again. Keep the original initial shuffle,
// updates and fresh Boost seed 4. Check sizes first and protect the empty R
// pointer before allocating draws, as in fisher_yates_sampler().
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

    // Copy the selected positions into v in extraction order.
    for (int i = 0; i < m; i ++) {
      v[i] = x[i+N];
    }

    // Extend the draw so each allowed treatment size can use its prefix.
    x[N] = N + m;
    for (int i = m + 1; i <= M; i++) {
      // Sample the next position for the larger treatment group.
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
  // Return the R pointer that will free the draws when no longer needed.
  return sceptre::sampling::register_bank(pointer, synth_idx_list);
}


namespace {
// Before allocating draws, check that probabilities are finite and in [0, 1],
// and that B and cell counts fit integer positions. Probabilities 0 and 1,
// no cells, and B = 0 are allowed. Invalid sizes error rather than being reduced.
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


// Generate B CRT draws, each a possibly empty list of cell positions starting at 0.
// Visit cells in input order. For each cell, draw how many of the B samples
// include it using Boost binomial(B, p), then choose those samples by Fisher-Yates
// sampling without replacement. Append the cell to each chosen sample.
// Keep this cell-by-cell order: looping over samples first gives different seeded
// draws. Each CRT worker generates B1+B2+B3 draws once per guide and reuses them
// across genes. A resampled treatment group need not have the observed size.
// Start a fresh Boost mt19937 at seed 4. As above, protect the empty R pointer
// and its cleanup handler before allocating C++ draws, then attach on success.
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
    // Draw how many of the B samples include this cell.
    M = binom_distribution(generator);
    // Choose M distinct draw rows, numbered 0 through B-1.
    draw_wor_sample(generator, unif_distribution, i_doub_array, x, B, M);
    // Include cell j in each chosen draw.
    for (int i = B - M; i < B; i ++) (*synth_idx_list)[x[i]].push_back(j);
  }

  return sceptre::sampling::register_bank(pointer, synth_idx_list);
}
