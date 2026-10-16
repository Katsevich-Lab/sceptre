#ifndef SCEPTRE_PAIR_ELIGIBILITY_H
#define SCEPTRE_PAIR_ELIGIBILITY_H

// The existing observed-count policy. Sampling uses this in its unchanged
// candidate stream; discovery uses the vector adapter in pair_eligibility.cpp.
// Counts and validated scalar thresholds share count units, with inclusive >=.
// "legacy" names the deliberately fixed statistical policy, not a retired alias.
inline bool legacy_pair_eligible(int treatment, int control,
                                 int treatment_threshold, int control_threshold) {
  return treatment >= treatment_threshold && control >= control_threshold;
}

// Preserve the old R vector predicate's missing-value behavior here too. Native
// calibration counts are complete; the vector boundary supplies R's NA sentinel.
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
