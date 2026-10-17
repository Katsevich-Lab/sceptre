#include <boost/random.hpp>
#include <boost/uuid/detail/sha1.hpp>
#include <Rcpp.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <vector>

// [[Rcpp::depends(BH)]]

// Immutable characterization oracle: the literal Git blob from
// 21f9ea098b69c444a884a49bc254de711968e3b5:src/generate_samples_functions.cpp.
// Its original includes are already loaded above, outside this namespace.
// Never include candidate package source here or regenerate this oracle from it.
namespace refactor_native_reference_21f9ea0 {
// A macro include prevents sourceCpp's attribute scanner from treating the
// frozen file's historical Rcpp::export comments as new global wrappers.
#define REFACTOR_NATIVE_REFERENCE_FILE "native-sampler-reference.inc"
#include REFACTOR_NATIVE_REFERENCE_FILE
#undef REFACTOR_NATIVE_REFERENCE_FILE
}

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
Rcpp::List refactor_native_reference_sample(std::string name, Rcpp::List args) {
    SEXP pointer = R_NilValue;
    if (name == "fisher_yates_sampler") {
        pointer = refactor_native_reference_21f9ea0::fisher_yates_sampler(
            Rcpp::as<int>(args["n_tot"]), Rcpp::as<int>(args["M"]), Rcpp::as<int>(args["B"])
        );
    } else if (name == "hybrid_fisher_iwor_sampler") {
        pointer = refactor_native_reference_21f9ea0::hybrid_fisher_iwor_sampler(
            Rcpp::as<int>(args["N"]), Rcpp::as<int>(args["m"]),
            Rcpp::as<int>(args["M"]), Rcpp::as<int>(args["B"])
        );
    } else if (name == "crt_index_sampler_fast") {
        pointer = refactor_native_reference_21f9ea0::crt_index_sampler_fast(
            Rcpp::as<Rcpp::NumericVector>(args["fitted_probabilities"]), Rcpp::as<int>(args["B"])
        );
    } else if (name == "crt_index_sampler") {
        pointer = refactor_native_reference_21f9ea0::crt_index_sampler(
            Rcpp::as<Rcpp::NumericVector>(args["fitted_probabilities"]), Rcpp::as<int>(args["B"])
        );
    } else Rcpp::stop("Unknown frozen-reference sampler.");
    Rcpp::Shield<SEXP> protected_pointer(pointer);
    return refactor_native_materialize(protected_pointer);
}

// Hash the actual oracle bytes using Git's blob-header convention. BH keeps
// this portable in package tests without requiring Git or another R package.
// [[Rcpp::export]]
std::string refactor_native_git_blob(std::string path) {
    std::ifstream stream(path.c_str(), std::ios::binary);
    if (!stream) Rcpp::stop("Cannot read the frozen sampler reference.");
    std::string contents((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
    std::string header = "blob " + std::to_string(contents.size());
    header.push_back('\0');
    boost::uuids::detail::sha1 hash;
    hash.process_bytes(header.data(), header.size());
    hash.process_bytes(contents.data(), contents.size());
    boost::uuids::detail::sha1::digest_type digest;
    hash.get_digest(digest);
    std::ostringstream result;
    result << std::hex << std::setfill('0');
    for (size_t i = 0; i < sizeof(digest) / sizeof(digest[0]); ++i) {
        result << std::setw(sizeof(digest[0]) * 2) << static_cast<unsigned long>(digest[i]);
    }
    return result.str();
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
