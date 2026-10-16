#define SCEPTRE_ORIGINAL_DISPERSION "/private/tmp/sceptre-dispersion-fault-proof-ef0c77b/original-instrumented.inc"
#define SCEPTRE_CANDIDATE_DISPERSION "/private/tmp/sceptre-dispersion-fault-proof-ef0c77b/candidate-instrumented.inc"
// Test-only interposition immediately before the actual MLE/MM calls. Each
// included source keeps its arithmetic, try boundary and state assignments.
// Never linked into the package; no production injection interface is added.
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/policies/policy.hpp>
#include <cerrno>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <system_error>
#include <type_traits>
#include <vector>

namespace dispersion_probe {
int failure_stage = 0;
int failure_type = 0;
std::vector<int> reached;

void before(int stage) {
  reached.push_back(stage);
  if (stage != failure_stage) return;
  switch (failure_type) {
  case 0: return;
  case 1: throw std::domain_error("Controlled numerical domain failure.");
  case 2: throw std::overflow_error("Controlled numerical overflow failure.");
  case 3: throw std::bad_alloc();
  case 4: throw std::system_error(std::error_code(EIO, std::generic_category()), "Controlled system failure.");
  case 5: throw std::length_error("Controlled length failure.");
  case 6: Rcpp::stop("Controlled Rcpp failure.");
  case 7: throw 719;
  default: Rcpp::stop("Unknown controlled failure type.");
  }
}
}

namespace original_dispersion {
#include SCEPTRE_ORIGINAL_DISPERSION
}
namespace candidate_dispersion {
#include SCEPTRE_CANDIDATE_DISPERSION
}

// Catch outside the actual estimator, recording propagation without a real
// allocation failure/OOM and without replacing its native numerical functions.
// [[Rcpp::export]]
Rcpp::List refactor_dispersion_case(bool candidate, Rcpp::NumericVector y,
                                   Rcpp::NumericVector mu, double dfr,
                                   int limit, double eps, int stage, int type) {
  dispersion_probe::failure_stage = stage;
  dispersion_probe::failure_type = type;
  dispersion_probe::reached.clear();
  Rcpp::List value;
  std::string escaped;
  try {
    value = candidate ? candidate_dispersion::estimate_theta(y, mu, dfr, limit, eps)
                      : original_dispersion::estimate_theta(y, mu, dfr, limit, eps);
  } catch (const std::domain_error&) { escaped = "domain";
  } catch (const std::overflow_error&) { escaped = "overflow";
  } catch (const std::bad_alloc&) { escaped = "allocation";
  } catch (const std::system_error&) { escaped = "system";
  } catch (const std::length_error&) { escaped = "length";
  } catch (const Rcpp::exception&) { escaped = "Rcpp";
  } catch (...) { escaped = "unknown"; }
  return Rcpp::List::create(Rcpp::Named("value") = value,
                           Rcpp::Named("escaped") = escaped,
                           Rcpp::Named("reached") = dispersion_probe::reached);
}

// These facts bind the reviewed default binary64 CPU dispatch on each CI lane.
// They describe the linked headers/configuration, not a new production policy.
// [[Rcpp::export]]
Rcpp::List refactor_dispersion_policy_facts() {
  using policy = boost::math::policies::policy<>;
  using evaluation = boost::math::policies::evaluation<double, policy>::type;
  using precision = boost::math::policies::precision<double, policy>::type;
  return Rcpp::List::create(
    Rcpp::Named("double_radix") = std::numeric_limits<double>::radix,
    Rcpp::Named("double_digits") = std::numeric_limits<double>::digits,
    Rcpp::Named("effective_precision") = precision::value,
    Rcpp::Named("evaluation_is_long_double") = std::is_same<evaluation, long double>::value,
    Rcpp::Named("evaluation_digits") = std::numeric_limits<evaluation>::digits,
    Rcpp::Named("domain") = static_cast<int>(policy::domain_error_type::value),
    Rcpp::Named("pole") = static_cast<int>(policy::pole_error_type::value),
    Rcpp::Named("overflow") = static_cast<int>(policy::overflow_error_type::value),
    Rcpp::Named("underflow") = static_cast<int>(policy::underflow_error_type::value),
    Rcpp::Named("denorm") = static_cast<int>(policy::denorm_error_type::value),
    Rcpp::Named("throw_on_error") = static_cast<int>(boost::math::policies::throw_on_error),
    Rcpp::Named("ignore_error") = static_cast<int>(boost::math::policies::ignore_error));
}
