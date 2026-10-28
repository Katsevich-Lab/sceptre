#include <Rcpp.h>
#include "eminem/eminem.hpp"
#include "tatami/sparse/compress_sparse_triplets.hpp"

#include <algorithm>
#include <climits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
S4 read_mtx(const std::string& path, bool transpose) {
  // Read a Matrix Market file into a dgCMatrix. eminem parses the file and
  // checks the banner, the dimensions and every entry; tatami sorts the
  // entries. With transpose, the rows of the file become columns (Parse writes
  // cells in rows; Cell Ranger writes features in rows).
  std::vector<int> rows, cols;
  std::vector<double> vals;
  int n_rows, n_cols;
  try {
    auto parser = eminem::parse_some_file<int>(path.c_str(), {});
    parser.scan_preamble();
    const auto& banner = parser.get_banner();
    if (banner.object != eminem::Object::MATRIX ||
        banner.format != eminem::Format::COORDINATE ||
        banner.symmetry != eminem::Symmetry::GENERAL ||
        !(banner.field == eminem::Field::INTEGER || banner.field == eminem::Field::REAL ||
          banner.field == eminem::Field::DOUBLE)) {
      throw std::runtime_error("expected a general coordinate matrix of integer or real counts");
    }
    if (parser.get_nlines() > static_cast<unsigned long long>(INT_MAX) ||
        parser.get_nrows() == INT_MAX || parser.get_ncols() == INT_MAX) {
      throw std::runtime_error("too many entries or dimensions for a sparse matrix in R");
    }
    n_rows = transpose ? parser.get_ncols() : parser.get_nrows();
    n_cols = transpose ? parser.get_nrows() : parser.get_ncols();
    rows.reserve(parser.get_nlines());
    cols.reserve(parser.get_nlines());
    vals.reserve(parser.get_nlines());
    parser.scan_real([&](int r, int c, double x) {
      if (transpose) std::swap(r, c);
      rows.push_back(r - 1);
      cols.push_back(c - 1);
      vals.push_back(x);
    });
  } catch (const std::exception& e) {
    stop(path + ": " + e.what());
  }

  // sort the entries by column, then by row, and get the column pointers
  auto p = tatami::compress_sparse_triplets<std::vector<double>, int>(n_cols, vals, cols, rows);
  // a dgCMatrix may not repeat a row within a column: Matrix functions would
  // then disagree about the value, so reject repeated entries
  for (int c = 0; c < n_cols; c ++) {
    if (std::adjacent_find(rows.begin() + p[c], rows.begin() + p[c + 1]) != rows.begin() + p[c + 1]) {
      stop(path + ": the same (row, column) entry appears more than once");
    }
  }

  S4 out("dgCMatrix");
  out.slot("Dim") = IntegerVector::create(n_rows, n_cols);
  out.slot("i") = IntegerVector(rows.begin(), rows.end());
  out.slot("p") = IntegerVector(p.begin(), p.end());
  out.slot("x") = NumericVector(vals.begin(), vals.end());
  return out;
}
