#include <R.h>
#include <Rinternals.h>
#include <time.h>

/* Shared benchmark-only clock: the same small DLL is loaded in both builds. */
SEXP sceptre_triplet_monotonic(void) {
    struct timespec time;
    if (clock_gettime(CLOCK_MONOTONIC, &time) != 0) Rf_error("Cannot read monotonic clock.");
    return Rf_ScalarReal((double) time.tv_sec + (double) time.tv_nsec / 1e9);
}
