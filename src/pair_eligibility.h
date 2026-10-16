#ifndef SCEPTRE_PAIR_ELIGIBILITY_H
#define SCEPTRE_PAIR_ELIGIBILITY_H

// The existing observed-count policy. Sampling uses this in its unchanged
// candidate stream; discovery uses the vector adapter in pair_eligibility.cpp.
inline bool legacy_pair_eligible(int treatment, int control,
                                 int treatment_threshold, int control_threshold) {
  return treatment >= treatment_threshold && control >= control_threshold;
}

#endif
