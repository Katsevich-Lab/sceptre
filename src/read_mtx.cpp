#include <Rcpp.h>
#include "eminem/eminem.hpp"

#include <algorithm>
#include <climits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
S4 read_mtx(const std::string& path, bool transpose) {
  // Read a Matrix Market file into a dgCMatrix. eminem checks the banner, the
  // dimensions and every entry. With transpose, the rows of the file become
  // columns (Parse writes cells in rows; Cell Ranger writes features in rows).
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

  // place the entries column by column, keeping file order within a column
  IntegerVector p(n_cols + 1), i(rows.size());
  NumericVector x(rows.size());
  for (int c : cols) p[c + 1] ++;
  for (int c = 0; c < n_cols; c ++) p[c + 1] += p[c];
  std::vector<int> next(p.begin(), p.end() - 1);
  for (size_t k = 0; k < rows.size(); k ++) {
    int dest = next[cols[k]]++;
    i[dest] = rows[k];
    x[dest] = vals[k];
  }

  // sort rows within each column (files are usually sorted already)
  std::vector<std::pair<int, double>> buffer;
  for (int c = 0; c < n_cols; c ++) {
    int start = p[c], end = p[c + 1];
    if (!std::is_sorted(i.begin() + start, i.begin() + end)) {
      buffer.clear();
      for (int k = start; k < end; k ++) buffer.emplace_back(i[k], x[k]);
      std::sort(buffer.begin(), buffer.end());
      for (int k = start; k < end; k ++) {
        i[k] = buffer[k - start].first;
        x[k] = buffer[k - start].second;
      }
    }
    if (std::adjacent_find(i.begin() + start, i.begin() + end) != i.begin() + end) {
      stop(path + ": the same (row, column) entry appears more than once");
    }
  }

  S4 out("dgCMatrix");
  out.slot("Dim") = IntegerVector::create(n_rows, n_cols);
  out.slot("i") = i;
  out.slot("p") = p;
  out.slot("x") = x;
  return out;
}
