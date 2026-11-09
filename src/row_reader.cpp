#include "Rcpp.h"
#include "Rtatami.h"

#include <memory>
#include <vector>
using namespace Rcpp;

// Reads rows of a matrix (responses x cells) one at a time through a single
// tatami extractor (a matrix from beachmat::initializeCpp()), so that a
// disk-backed matrix keeps its recently read chunks between reads: rows read
// in order of their position are mostly read from memory. The extractor holds
// state, so each worker creates its own.
namespace {

struct RowExtractor {
  Rtatami::BoundNumericPointer matrix;
  std::unique_ptr<tatami::MyopicSparseExtractor<double, int> > extractor;
  std::vector<double> values;
  std::vector<int> indices;
};

}

// [[Rcpp::export]]
SEXP make_row_extractor(SEXP matrix) {
  RowExtractor* out = new RowExtractor{Rtatami::BoundNumericPointer(matrix), nullptr, {}, {}};
  const auto& mat = out->matrix->ptr;
  out->extractor = mat->sparse(true, tatami::Options());
  out->values.resize(mat->ncol());
  out->indices.resize(mat->ncol());
  return XPtr<RowExtractor>(out, true);
}

// Row `row` (1-based) as a dense numeric vector.
// [[Rcpp::export]]
NumericVector read_extractor_row(SEXP row_extractor, int row) {
  XPtr<RowExtractor> r(row_extractor);
  NumericVector out(r->values.size());
  auto range = r->extractor->fetch(row - 1, r->values.data(), r->indices.data());
  for (int k = 0; k < range.number; k ++) out[range.index[k]] = range.value[k];
  return out;
}
