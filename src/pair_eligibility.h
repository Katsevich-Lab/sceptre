#ifndef SCEPTRE_PAIR_ELIGIBILITY_H
#define SCEPTRE_PAIR_ELIGIBILITY_H

// A pair passes QC when treatment and control counts both reach their thresholds
// (>=). Calibration sampling calls this rule directly; discovery/QC uses the R
// vector wrapper in pair_eligibility.cpp. Keep the same rule in both places.
// The name legacy refers to the retained statistical rule, not unused code.
inline bool legacy_pair_eligible(int treatment, int control,
                                 int treatment_threshold, int control_threshold) {
  return treatment >= treatment_threshold && control >= control_threshold;
}

// Handle missing R counts like vector comparisons: a definite failure takes
// precedence over a missing count. Calibration counts themselves are complete.
inline int legacy_pair_eligibility_value(int treatment, int control,
                                         int treatment_threshold, int control_threshold,
                                         int missing) {
  if (treatment != missing && control != missing) {
    return legacy_pair_eligible(treatment, control,
                                treatment_threshold, control_threshold);
  }
  if ((treatment != missing && treatment < treatment_threshold) ||
      (control != missing && control < control_threshold)) return 0;
  return missing;
}

#endif
