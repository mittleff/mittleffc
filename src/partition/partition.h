#ifndef PARTITION_H
#define PARTITION_H
#include <flint/arb.h>
#include <flint/acb.h>

/**
 * Check if angle a is clockwise from angle b
 *
 * Uses the cross product (sine of angle difference) to determine orientation.
 * Matches Scheme implementation: sin(a - b) > 0 means clockwise
 *
 * @param a      First angle
 * @param b      Second angle
 * @param closed If true, use >= 0 instead of > 0
 * @return 1 if clockwise, 0 otherwise
/*  *\/ */
/* int ml_clockwise(const arb_t a, const arb_t b, int closed); */

/**
 * Determine which region z belongs to
 *
 * Algorithm:
 * 1. If |z| <= taylor_radius -> Region G0 (for Taylor series)
 * 2. Compute threshold r1
 * 3. If |z| >= r1 -> Check regions G1-G4 (for asymptotic expansions)
 * 4. Else -> Check regions G5-G6 (for integral representations)
 *
 * @param z             Complex argument
 * @param alpha         Parameter α
 * @param taylor_radius Threshold for Region G0
 * @param prec          Precision in bits
 * @return Region number (0-6), or -1 if no region matched (error)
 */
int ml_compute_region(const acb_t z, const acb_t alpha, arb_t taylor_radius,
                      slong prec);

#endif /* PARTITION_H */
