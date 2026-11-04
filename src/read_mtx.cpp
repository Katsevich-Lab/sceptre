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

// Reads the preamble of a Matrix Market file opened with eminem, which checks
// the banner, the dimensions and every entry, and checks that the file holds a
// general coordinate matrix of counts that fits a sparse matrix in R.
template <class Parser>
static void check_mtx_preamble(Parser& parser) {
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
}

// read_mtx() with transpose: the rows of the file become the columns of the
// result. The file is read twice, as ondisc does, so that only the result is
// held in memory: the first pass counts the entries of each column, the second
// places each entry in its column. The entries of a column keep their order in
// the file, which lists them by increasing row when sorted by column then row
// (as Cell Ranger writes them); a column out of order is then sorted.
static S4 read_mtx_transposed(const std::string& path) {
  int n_rows = 0, n_cols = 0;
  std::vector<int> p;
  try {
    auto parser = eminem::parse_some_file<int>(path.c_str(), {});
    check_mtx_preamble(parser);
    n_rows = parser.get_ncols();
    n_cols = parser.get_nrows();
    p.assign(static_cast<size_t>(n_cols) + 1, 0);
    parser.scan_real([&](int r, int, double) { p[r] ++; });
  } catch (const std::exception& e) {
    stop(path + ": " + e.what());
  }
  for (int c = 0; c < n_cols; c ++) p[c + 1] += p[c];
  IntegerVector i(p[n_cols]);
  NumericVector x(p[n_cols]);
  std::vector<int> next(p.begin(), p.end() - 1);
  try {
    auto parser = eminem::parse_some_file<int>(path.c_str(), {});
    check_mtx_preamble(parser);
    if (parser.get_nrows() != n_cols || parser.get_ncols() != n_rows ||
        parser.get_nlines() != static_cast<unsigned long long>(p[n_cols])) {
      throw std::runtime_error("the file changed while it was read");
    }
    parser.scan_real([&](int r, int c, double val) {
      int k = next[r - 1] ++;
      i[k] = c - 1;
      x[k] = val;
    });
  } catch (const std::exception& e) {
    stop(path + ": " + e.what());
  }
  std::vector<std::pair<int, double>> column;
  for (int c = 0; c < n_cols; c ++) {
    if (std::is_sorted(i.begin() + p[c], i.begin() + p[c + 1])) continue;
    column.clear();
    for (int k = p[c]; k < p[c + 1]; k ++) column.emplace_back(i[k], x[k]);
    std::stable_sort(column.begin(), column.end(), [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
      return a.first < b.first;
    });
    for (int k = p[c]; k < p[c + 1]; k ++) {
      i[k] = column[k - p[c]].first;
      x[k] = column[k - p[c]].second;
    }
  }
  for (int c = 0; c < n_cols; c ++) {
    if (std::adjacent_find(i.begin() + p[c], i.begin() + p[c + 1]) != i.begin() + p[c + 1]) {
      stop(path + ": the same (row, column) entry appears more than once");
    }
  }
  S4 out("dgCMatrix");
  out.slot("Dim") = IntegerVector::create(n_rows, n_cols);
  out.slot("i") = i;
  out.slot("p") = IntegerVector(p.begin(), p.end());
  out.slot("x") = x;
  return out;
}

// [[Rcpp::export]]
S4 read_mtx(const std::string& path, bool transpose) {
  // Read a Matrix Market file into a dgCMatrix. eminem parses the file and
  // checks the banner, the dimensions and every entry; tatami sorts the
  // entries. With transpose, the rows of the file become columns (Parse writes
  // cells in rows; Cell Ranger writes features in rows), without tatami
  // (read_mtx_transposed()).
  if (transpose) return read_mtx_transposed(path);
  std::vector<int> rows, cols;
  std::vector<double> vals;
  int n_rows, n_cols;
  try {
    auto parser = eminem::parse_some_file<int>(path.c_str(), {});
    check_mtx_preamble(parser);
    n_rows = parser.get_nrows();
    n_cols = parser.get_ncols();
    rows.reserve(parser.get_nlines());
    cols.reserve(parser.get_nlines());
    vals.reserve(parser.get_nlines());
    parser.scan_real([&](int r, int c, double x) {
      rows.push_back(r - 1);
      cols.push_back(c - 1);
      vals.push_back(x);
    });
  } catch (const std::exception& e) {
    stop(path + ": " + e.what());
  }

  // sort the entries by column, then by row, and get the column pointers
  // the template parameters are deduced: their order differs between the
  // tatami versions of assorthead 1.6 and 1.7
  auto p = tatami::compress_sparse_triplets(n_cols, vals, cols, rows);
  // a dgCMatrix may not repeat a row within a column: Matrix functions would
  // then disagree about the value, so reject repeated entries
  for (int c = 0; c < n_cols; c ++) {
    if (std::adjacent_find(rows.begin() + p[c], rows.begin() + p[c + 1]) != rows.begin() + p[c + 1]) {
      stop(path + ": the same (row, column) entry appears more than once");
    }
  }

  // free each triplet vector once it is copied (or no longer needed), so that
  // at most one is held alongside its copy
  std::vector<int>().swap(cols);
  S4 out("dgCMatrix");
  out.slot("Dim") = IntegerVector::create(n_rows, n_cols);
  out.slot("i") = IntegerVector(rows.begin(), rows.end());
  std::vector<int>().swap(rows);
  out.slot("p") = IntegerVector(p.begin(), p.end());
  out.slot("x") = NumericVector(vals.begin(), vals.end());
  std::vector<double>().swap(vals);
  return out;
}
