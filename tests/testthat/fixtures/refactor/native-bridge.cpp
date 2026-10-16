#include <Rcpp.h>
#include <vector>

// Test-only access to production draw banks. This file is not package source.
// Draw entries use the native zero-based convention.

// [[Rcpp::export]]
Rcpp::List refactor_native_materialize(SEXP pointer) {
    Rcpp::XPtr<std::vector<std::vector<int>>> draws(pointer);
    Rcpp::List out(draws->size());
    for (size_t i = 0; i < draws->size(); ++i) out[i] = (*draws)[i];
    return out;
}

// [[Rcpp::export]]
SEXP refactor_native_pointer(Rcpp::List draws, int n_cells) {
    if (n_cells < 1 || draws.size() < 1) {
        Rcpp::stop("A positive cell count and at least one draw are required.");
    }
    std::vector<std::vector<int>> validated(draws.size());
    for (R_xlen_t i = 0; i < draws.size(); ++i) {
        SEXP entry = draws[i];
        if (TYPEOF(entry) != INTSXP) {
            Rcpp::stop("Draw indices must be integer vectors.");
        }
        Rcpp::IntegerVector positions(entry);
        std::vector<bool> seen(n_cells, false);
        for (int position : positions) {
            if (position == NA_INTEGER || position < 0 || position >= n_cells) {
                Rcpp::stop("Draw index outside the native zero-based cell range.");
            }
            if (seen[position]) Rcpp::stop("Draw indices must be unique.");
            seen[position] = true;
            validated[i].push_back(position);
        }
    }
    return Rcpp::XPtr<std::vector<std::vector<int>>>(
        new std::vector<std::vector<int>>(std::move(validated)), true
    );
}
