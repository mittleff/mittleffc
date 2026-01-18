#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <flint/arb.h>
#include <flint/acb.h>

/**
 * @file algorithm.h
 * @brief Core algorithms for computing the Mittag-Leffler function
 *
 * This module implements three main computational strategies based on
 * the Scheme reference implementation:
 *
 * 1. Taylor series expansion (for small |z|)
 * 2. Asymptotic expansion (for large |z| in regions G1-G4)
 * 3. Integral representation (for regions G5-G6)
 */

/**
 * Compute Mittag-Leffler using Taylor series expansion
 *
 * Implements the series:
 *   E_{a,b}(z) = Σ_{k=0}^{kmax} z^k / Γ(ak + b)
 *
 * This method is used for small |z| (region G0, typically |z| <= taylor_radius).
 *
 * @param result    Output: computed value
 * @param a         Parameter α (alpha)
 * @param b         Parameter β (beta)
 * @param z         Complex argument
 * @param prec      Precision in bits
 */
void ml_taylor(acb_t result, const acb_t a, const acb_t b, const acb_t z, slong prec);

/**
 * Compute Mittag-Leffler using asymptotic expansion
 *
 * Implements the asymptotic formula:
 *   E_{a,b}(z) = (ρ/a) * z^((1-b)/a) * exp(z^(1/a)) - Σ_{k=1}^{kmax} z^(-k) / Γ(b - ak)
 *
 * where ρ depends on the region (G1-G4) and involves erfc for regions G3-G4.
 *
 * This method is used for large |z| in regions G1-G4.
 *
 * @param result       Output: computed value
 * @param a            Parameter α (alpha)
 * @param b            Parameter β (beta)
 * @param z            Complex argument
 * @param region       Region number (1-4)
 * @param prec         Precision in bits
 */
void ml_asymptotics(acb_t result, const acb_t a, const acb_t b, const acb_t z,
                    int region, slong prec);

/**
 * Compute Mittag-Leffler using integral representation
 *
 * Implements the Hankel contour integral representation:
 *   E_{a,b}(z) = A_i + B_i + C_i
 *
 * where:
 *   - A_i is a direct computation term
 *   - B_i is an integral along a ray
 *   - C_i is an integral along a circular arc
 *
 * This method is used for intermediate |z| in regions G5-G6.
 *
 * @param result       Output: computed value
 * @param a            Parameter α (alpha)
 * @param b            Parameter β (beta)
 * @param z            Complex argument
 * @param region       Region number (5 or 6)
 * @param prec         Precision in bits
 */
void ml_integral_rep(acb_t result, const acb_t a, const acb_t b, const acb_t z,
                     int region, slong prec);

/**
 * Core series computation (used by both taylor and asymptotics)
 *
 * Computes:
 *   Σ_{k=kmin}^{kmax} z^(s*k) / Γ(a*s*k + b)
 *
 * where s = +1 for Taylor series, s = -1 for asymptotic series.
 *
 * @param result       Output: computed series sum
 * @param a            Parameter α (alpha)
 * @param b            Parameter β (beta)
 * @param z            Complex argument
 * @param asymptotic   If true, compute asymptotic series (s=-1, kmin=1)
 *                     If false, compute Taylor series (s=+1, kmin=0)
 * @param prec         Precision in bits
 */
void ml_series(acb_t result, const acb_t a, const acb_t b, const acb_t z,
               int asymptotic, slong prec);

#endif /* ALGORITHM_H */
