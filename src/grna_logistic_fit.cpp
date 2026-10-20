// Reused from ZiangNiu6/sceptre, commit 7ea1aefd246f87965a16a1cfc5a4bf5b1901e05a.
// GPL-3; experimental integration retains the original numerical solver.
// Exact unpenalized logistic fitter for SCEPTRE's X | Z nuisance model.
//
// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

using Eigen::MatrixXd;
using Eigen::VectorXd;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

double max_abs(const VectorXd& values) {
  return values.size() == 0 ? 0.0 : values.cwiseAbs().maxCoeff();
}

struct Evaluation {
  bool valid = false;
  double loglik = kNaN;
  double gradient_max = kNaN;
  double max_abs_eta = kNaN;
  VectorXd gradient;
  MatrixXd information;
  VectorXd fitted;
};

template <int FixedP>
bool evaluate_impl(const Rcpp::NumericMatrix& X,
                   const VectorXd& beta,
                   const VectorXd& Xty,
                   Evaluation& out) {
  const int n = X.nrow();
  const int p = FixedP > 0 ? FixedP : X.ncol();
  out.valid = false;
  out.loglik = kNaN;
  out.gradient_max = kNaN;
  out.max_abs_eta = kNaN;
  if (X.ncol() != p || beta.size() != p || Xty.size() != p ||
      !beta.allFinite()) {
    return false;
  }

  const double* x = X.begin();
  out.gradient = Xty;
  out.information.setZero(p, p);
  out.fitted.resize(n);
  std::array<double, 16> stack_row;
  std::vector<double> heap_row;
  double* row = stack_row.data();
  if (p > static_cast<int>(stack_row.size())) {
    heap_row.resize(p);
    row = heap_row.data();
  }
  double log_partition = 0.0;
  double max_eta = 0.0;

  for (int i = 0; i < n; ++i) {
    double eta = 0.0;
    for (int j = 0; j < p; ++j) {
      row[j] = x[i + n * j];
      eta += beta[j] * row[j];
    }
    if (!std::isfinite(eta)) return false;
    const double exponential = std::exp(-std::abs(eta));
    const double probability = eta >= 0.0
                                   ? 1.0 / (1.0 + exponential)
                                   : exponential / (1.0 + exponential);
    const double variance = probability * (1.0 - probability);
    if (!std::isfinite(probability) || !std::isfinite(variance)) return false;
    out.fitted[i] = probability;
    log_partition += eta >= 0.0
                         ? eta + std::log1p(exponential)
                         : std::log1p(exponential);
    max_eta = std::max(max_eta, std::abs(eta));
    for (int j = 0; j < p; ++j) {
      const double xj = row[j];
      out.gradient[j] -= probability * xj;
      for (int k = 0; k <= j; ++k) {
        out.information(j, k) += variance * xj * row[k];
      }
    }
  }
  for (int j = 0; j < p; ++j) {
    for (int k = 0; k < j; ++k) {
      out.information(k, j) = out.information(j, k);
    }
  }
  out.loglik = beta.dot(Xty) - log_partition;
  out.gradient_max = max_abs(out.gradient);
  out.max_abs_eta = max_eta;
  out.valid = std::isfinite(out.loglik) &&
              std::isfinite(out.gradient_max) &&
              out.gradient.allFinite() && out.information.allFinite();
  return out.valid;
}

bool evaluate_into(const Rcpp::NumericMatrix& X,
                   const VectorXd& beta,
                   const VectorXd& Xty,
                   Evaluation& out) {
  return X.ncol() == 6 ? evaluate_impl<6>(X, beta, Xty, out)
                       : evaluate_impl<0>(X, beta, Xty, out);
}

Evaluation prevalence_state(const int n,
                            const int n_treated,
                            const double prevalence,
                            const VectorXd& Xty,
                            const Rcpp::NumericVector& column_sums,
                            const Rcpp::NumericMatrix& crossprod) {
  const int p = Xty.size();
  Evaluation out;
  out.gradient.resize(p);
  out.information.resize(p, p);
  const double variance = prevalence * (1.0 - prevalence);
  for (int j = 0; j < p; ++j) {
    out.gradient[j] = Xty[j] - prevalence * column_sums[j];
    for (int k = 0; k < p; ++k) {
      out.information(j, k) = variance * crossprod(j, k);
    }
  }
  out.loglik = static_cast<double>(n_treated) * std::log(prevalence) +
               static_cast<double>(n - n_treated) * std::log1p(-prevalence);
  out.gradient_max = max_abs(out.gradient);
  out.max_abs_eta = std::abs(std::log(prevalence / (1.0 - prevalence)));
  out.valid = std::isfinite(out.loglik) && out.gradient.allFinite() &&
              out.information.allFinite();
  return out;
}

bool solve_information(const MatrixXd& information,
                       const VectorXd& gradient,
                       const double rank_tolerance,
                       VectorXd& step,
                       double& condition_ratio) {
  if (!information.allFinite() || !gradient.allFinite()) return false;
  Eigen::LDLT<MatrixXd> decomposition(information);
  if (decomposition.info() != Eigen::Success || !decomposition.isPositive()) {
    return false;
  }
  const VectorXd diagonal = decomposition.vectorD().cwiseAbs();
  if (diagonal.size() == 0 || !diagonal.allFinite()) return false;
  const double largest = diagonal.maxCoeff();
  const double smallest = diagonal.minCoeff();
  condition_ratio = largest > 0.0 ? smallest / largest : 0.0;
  if (!(largest > 0.0) || !(smallest > rank_tolerance * largest)) {
    return false;
  }
  step = decomposition.solve(gradient);
  return decomposition.info() == Eigen::Success && step.allFinite();
}

VectorXd sparse_Xty(const Rcpp::NumericMatrix& X,
                    const Rcpp::IntegerVector& treated_indices,
                    std::vector<unsigned char>& treated) {
  const int n = X.nrow();
  const int p = X.ncol();
  const double* x = X.begin();
  VectorXd result = VectorXd::Zero(p);
  treated.assign(n, 0);
  for (R_xlen_t k = 0; k < treated_indices.size(); ++k) {
    const int value = treated_indices[k];
    if (value == NA_INTEGER || value < 1 || value > n) {
      Rcpp::stop("treated_indices must contain valid one-based indices");
    }
    const int i = value - 1;
    if (treated[i]) Rcpp::stop("treated_indices cannot contain duplicates");
    treated[i] = 1;
    for (int j = 0; j < p; ++j) result[j] += x[i + n * j];
  }
  return result;
}

bool completely_separated(const VectorXd& fitted,
                          const std::vector<unsigned char>& treated) {
  if (fitted.size() != static_cast<Eigen::Index>(treated.size())) return false;
  double minimum_treated = 1.0;
  double maximum_control = 0.0;
  for (Eigen::Index i = 0; i < fitted.size(); ++i) {
    if (treated[i]) {
      minimum_treated = std::min(minimum_treated, fitted[i]);
    } else {
      maximum_control = std::max(maximum_control, fitted[i]);
    }
  }
  return minimum_treated > 0.5 && maximum_control < 0.5;
}

VectorXd original_coefficients(const VectorXd& beta,
                               const Rcpp::NumericVector& centers,
                               const Rcpp::NumericVector& scales,
                               const bool conditioned) {
  if (!conditioned) return beta;
  VectorXd out(beta.size());
  out[0] = beta[0];
  for (Eigen::Index j = 1; j < beta.size(); ++j) {
    out[j] = beta[j] / scales[j];
    out[0] -= centers[j] * out[j];
  }
  return out;
}

VectorXd original_score(const VectorXd& score,
                        const Rcpp::NumericVector& centers,
                        const Rcpp::NumericVector& scales,
                        const bool conditioned) {
  if (!conditioned) return score;
  VectorXd out(score.size());
  out[0] = score[0];
  for (Eigen::Index j = 1; j < score.size(); ++j) {
    out[j] = scales[j] * score[j] + centers[j] * score[0];
  }
  return out;
}

Rcpp::List solver_result(const VectorXd& beta,
                         const Evaluation& evaluation,
                         const Rcpp::NumericVector& centers,
                         const Rcpp::NumericVector& scales,
                         const bool conditioned,
                         const bool converged,
                         const std::string& reason,
                         const int iterations,
                         const int full_scans,
                         const int step_halvings,
                         const int rank,
                         const double condition_ratio,
                         const double original_score_tolerance,
                         const double conditioned_score_tolerance,
                         const double relative_deviance_change,
                         const bool separation_detected) {
  const VectorXd coefficients = original_coefficients(
      beta, centers, scales, conditioned);
  const VectorXd score = evaluation.gradient.size() == beta.size()
                             ? original_score(evaluation.gradient, centers,
                                              scales, conditioned)
                             : VectorXd::Constant(beta.size(), kNaN);
  return Rcpp::List::create(
      Rcpp::Named("coefficients") = Rcpp::wrap(coefficients),
      Rcpp::Named("fitted.values") = Rcpp::wrap(evaluation.fitted),
      Rcpp::Named("score") = Rcpp::wrap(score),
      Rcpp::Named("converged") = converged,
      Rcpp::Named("iterations") = iterations,
      Rcpp::Named("reason") = reason,
      Rcpp::Named("loglik") = evaluation.loglik,
      Rcpp::Named("deviance") = -2.0 * evaluation.loglik,
      Rcpp::Named("gradient_max") = max_abs(score),
      Rcpp::Named("score_tolerance") = original_score_tolerance,
      Rcpp::Named("conditioned_gradient_max") = evaluation.gradient_max,
      Rcpp::Named("conditioned_score_tolerance") =
          conditioned_score_tolerance,
      Rcpp::Named("relative_deviance_change") = relative_deviance_change,
      Rcpp::Named("max_abs_eta") = evaluation.max_abs_eta,
      Rcpp::Named("rank") = rank,
      Rcpp::Named("information_condition_ratio") = condition_ratio,
      Rcpp::Named("full_cell_scans") = full_scans,
      Rcpp::Named("step_halvings") = step_halvings,
      Rcpp::Named("conditioned") = conditioned,
      Rcpp::Named("separation_detected") = separation_detected,
      Rcpp::Named("used_analytic_prevalence_state") = true);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List prepare_grna_logistic_design_cpp(
    const Rcpp::NumericMatrix& Z,
    const bool condition = true,
    const double rank_tolerance = 1e-12) {
  const int n = Z.nrow();
  const int p = Z.ncol();
  if (n < 2 || p < 1) {
    Rcpp::stop("Z must have at least two rows and one column");
  }
  if (!std::isfinite(rank_tolerance) || rank_tolerance <= 0.0) {
    Rcpp::stop("rank_tolerance must be positive and finite");
  }

  Rcpp::NumericVector centers(p, 0.0);
  Rcpp::NumericVector scales(p, 1.0);
  bool finite = true;
  bool leading_intercept = true;
  for (int i = 0; i < n; ++i) {
    leading_intercept = leading_intercept &&
                        std::isfinite(Z(i, 0)) && Z(i, 0) == 1.0;
    for (int j = 0; j < p; ++j) finite = finite && std::isfinite(Z(i, j));
  }

  Rcpp::NumericMatrix X = condition ? Rcpp::NumericMatrix(n, p)
                                    : Rcpp::NumericMatrix(Z);
  std::string reason = "ok";
  if (!finite) reason = "nonfinite_design";
  if (finite && !leading_intercept) reason = "missing_leading_intercept";

  if (finite && leading_intercept && condition) {
    for (int i = 0; i < n; ++i) X(i, 0) = 1.0;
    for (int j = 1; j < p; ++j) {
      double sum = 0.0;
      for (int i = 0; i < n; ++i) sum += Z(i, j);
      centers[j] = sum / static_cast<double>(n);
      double square_sum = 0.0;
      double max_value = 0.0;
      for (int i = 0; i < n; ++i) {
        const double centered = Z(i, j) - centers[j];
        square_sum += centered * centered;
        max_value = std::max(max_value, std::abs(Z(i, j)));
      }
      scales[j] = std::sqrt(square_sum / static_cast<double>(n));
      if (!std::isfinite(scales[j]) ||
          scales[j] <= std::sqrt(std::numeric_limits<double>::epsilon()) *
                           std::max(1.0, max_value)) {
        reason = "rank_deficient_design";
        break;
      }
      for (int i = 0; i < n; ++i) {
        X(i, j) = (Z(i, j) - centers[j]) / scales[j];
      }
    }
  }

  Rcpp::NumericVector column_sums(p, kNaN);
  Rcpp::NumericMatrix crossprod(p, p);
  std::fill(crossprod.begin(), crossprod.end(), kNaN);
  int rank = 0;
  if (reason == "ok") {
    const Eigen::Map<const MatrixXd> X_map(X.begin(), n, p);
    const VectorXd sums = X_map.colwise().sum();
    const MatrixXd gram = X_map.transpose() * X_map;
    std::copy(sums.data(), sums.data() + p, column_sums.begin());
    std::copy(gram.data(), gram.data() + p * p, crossprod.begin());
    Eigen::SelfAdjointEigenSolver<MatrixXd> decomposition(gram);
    if (decomposition.info() != Eigen::Success ||
        !decomposition.eigenvalues().allFinite()) {
      reason = "rank_check_failed";
    } else {
      const double largest = decomposition.eigenvalues().maxCoeff();
      const double threshold = rank_tolerance * std::max(1.0, largest);
      rank = static_cast<int>((decomposition.eigenvalues().array() > threshold)
                                  .count());
      if (rank != p) reason = "rank_deficient_design";
    }
  }
  if (!Rf_isNull(Z.attr("dimnames"))) X.attr("dimnames") = Z.attr("dimnames");

  return Rcpp::List::create(
      Rcpp::Named("X") = X,
      Rcpp::Named("centers") = centers,
      Rcpp::Named("scales") = scales,
      Rcpp::Named("column_sums") = column_sums,
      Rcpp::Named("crossprod") = crossprod,
      Rcpp::Named("conditioned") = condition,
      Rcpp::Named("valid") = reason == "ok",
      Rcpp::Named("reason") = reason,
      Rcpp::Named("rank") = rank,
      Rcpp::Named("n") = n,
      Rcpp::Named("p") = p);
}

// [[Rcpp::export]]
Rcpp::List fit_grna_logistic_prepared_cpp(
    const Rcpp::List& prepared,
    const Rcpp::IntegerVector& treated_indices,
    const double tolerance = 1e-8,
    const int max_iterations = 25,
    const int max_step_halvings = 20,
    const double separation_eta_limit = 30.0,
    const double rank_tolerance = 1e-12) {
  if (!std::isfinite(tolerance) || tolerance <= 0.0) {
    Rcpp::stop("tolerance must be positive and finite");
  }
  if (max_iterations < 0 || max_step_halvings < 0) {
    Rcpp::stop("iteration limits must be nonnegative");
  }
  if (!std::isfinite(separation_eta_limit) || separation_eta_limit <= 0.0 ||
      !std::isfinite(rank_tolerance) || rank_tolerance <= 0.0) {
    Rcpp::stop("separation_eta_limit and rank_tolerance must be positive");
  }

  const Rcpp::NumericMatrix X = Rcpp::as<Rcpp::NumericMatrix>(prepared["X"]);
  const Rcpp::NumericVector centers =
      Rcpp::as<Rcpp::NumericVector>(prepared["centers"]);
  const Rcpp::NumericVector scales =
      Rcpp::as<Rcpp::NumericVector>(prepared["scales"]);
  const Rcpp::NumericVector column_sums =
      Rcpp::as<Rcpp::NumericVector>(prepared["column_sums"]);
  const Rcpp::NumericMatrix crossprod =
      Rcpp::as<Rcpp::NumericMatrix>(prepared["crossprod"]);
  const bool conditioned = Rcpp::as<bool>(prepared["conditioned"]);
  const bool design_valid = Rcpp::as<bool>(prepared["valid"]);
  const int rank = Rcpp::as<int>(prepared["rank"]);
  const int n = X.nrow();
  const int p = X.ncol();
  // Reject forged/empty prepared payloads before accessing the intercept.
  if (n < 2 || p < 1) {
    Rcpp::stop("prepared design requires at least two rows and one column");
  }
  if (centers.size() != p || scales.size() != p ||
      column_sums.size() != p || crossprod.nrow() != p ||
      crossprod.ncol() != p) {
    Rcpp::stop("prepared design has incompatible components");
  }

  std::vector<unsigned char> treated;
  const VectorXd Xty = sparse_Xty(X, treated_indices, treated);
  const int n_treated = treated_indices.size();
  VectorXd beta = VectorXd::Zero(p);
  Evaluation current;
  if (!design_valid || rank != p) {
    current.gradient = VectorXd::Constant(p, kNaN);
    current.fitted.resize(0);
    return solver_result(
        beta, current, centers, scales, conditioned, false,
        Rcpp::as<std::string>(prepared["reason"]), 0, 0, 0, rank, kNaN,
        kNaN, kNaN, kNaN, false);
  }
  if (n_treated == 0 || n_treated == n) {
    current.gradient = VectorXd::Constant(p, kNaN);
    current.fitted.resize(0);
    return solver_result(beta, current, centers, scales, conditioned, false,
                         "degenerate_binary_response", 0, 0, 0, rank,
                         kNaN, kNaN, kNaN, kNaN, true);
  }

  const double prevalence = static_cast<double>(n_treated) /
                            static_cast<double>(n);
  beta[0] = std::log(prevalence / (1.0 - prevalence));
  current = prevalence_state(n, n_treated, prevalence, Xty, column_sums,
                             crossprod);
  const double conditioned_score_reference = std::max(1.0, max_abs(Xty));
  const double conditioned_score_tolerance =
      tolerance * conditioned_score_reference;
  const VectorXd original_Xty = original_score(
      Xty, centers, scales, conditioned);
  const double original_score_tolerance =
      tolerance * std::max(1.0, max_abs(original_Xty));
  int iterations = 0;
  int full_scans = 0;
  int step_halvings = 0;
  double condition_ratio = kNaN;
  double relative_deviance_change = kNaN;
  std::string reason = "maximum_iterations_reached";
  bool converged = current.valid && current.gradient_max == 0.0;
  Evaluation candidate;
  VectorXd step(p);
  VectorXd candidate_beta(p);

  for (int iteration = 1; !converged && iteration <= max_iterations;
       ++iteration) {
    if (!solve_information(current.information, current.gradient,
                           rank_tolerance, step, condition_ratio)) {
      reason = "rank_deficient_information";
      break;
    }
    double scale = 1.0;
    bool accepted = false;
    for (int halving = 0; halving <= max_step_halvings; ++halving) {
      candidate_beta = beta + scale * step;
      evaluate_into(X, candidate_beta, Xty, candidate);
      ++full_scans;
      const double allowance = 1e-12 * (1.0 + std::abs(current.loglik));
      if (candidate.valid &&
          candidate.loglik + allowance >= current.loglik) {
        accepted = true;
        break;
      }
      if (halving < max_step_halvings) {
        scale *= 0.5;
        ++step_halvings;
      }
    }
    if (!accepted) {
      reason = "step_halving_failed";
      break;
    }
    const double previous_deviance = -2.0 * current.loglik;
    const double candidate_deviance = -2.0 * candidate.loglik;
    relative_deviance_change =
        std::abs(candidate_deviance - previous_deviance) /
        (0.1 + std::abs(candidate_deviance));
    beta.swap(candidate_beta);
    std::swap(current, candidate);
    iterations = iteration;
    if (max_abs(beta) > separation_eta_limit ||
        current.max_abs_eta > separation_eta_limit) {
      reason = "suspected_separation";
      break;
    }
    converged = relative_deviance_change < tolerance &&
                current.gradient_max <= conditioned_score_tolerance;
  }

  if (converged && current.fitted.size() == 0) {
    current.fitted = VectorXd::Constant(n, prevalence);
  }
  bool separation_detected = false;
  if (current.fitted.size() == n) {
    separation_detected = completely_separated(current.fitted, treated);
    if (separation_detected) {
      converged = false;
      reason = "complete_separation";
    }
  }
  if (converged) reason = "ok";
  return solver_result(beta, current, centers, scales, conditioned,
                       converged, reason, iterations, full_scans,
                       step_halvings, rank, condition_ratio,
                       original_score_tolerance,
                       conditioned_score_tolerance,
                       relative_deviance_change, separation_detected);
}
