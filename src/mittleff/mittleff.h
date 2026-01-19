/**
 * @file mittleff.h
 * @brief Mittleff C Library - Main Mittag-Leffler Function
 *
 * This header provides the main function for computing the Mittag-Leffler function
 * E_{α,β}(z), a generalization of the exponential function that appears frequently
 * in fractional calculus and anomalous diffusion processes.
 *
 * @version 0.1
 * @author Victor P. Santos
 * @copyright MIT License
 */

#ifndef ML_MITTLEFF_H
#define ML_MITTLEFF_H

#include <flint/arb.h>
#include <flint/acb.h>

/* Symbol visibility macros */
#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_MITTLEFF
    #define MITTLEFF_EXPORT __declspec(dllexport)
  #else
    #define MITTLEFF_EXPORT __declspec(dllimport)
  #endif
#else
  #if __GNUC__ >= 4
    #define MITTLEFF_EXPORT __attribute__((visibility("default")))
  #else
    #define MITTLEFF_EXPORT
  #endif
#endif

/* Default constants */
#define ML_DEFAULT_PRECISION 53
#define ML_DEFAULT_TAYLOR_RADIUS 0.95

/**
 * @brief Compute the Mittag-Leffler function E_{α,β}(z)
 *
 * The Mittag-Leffler function is a generalization of the exponential function:
 *   E_{α,β}(z) = Σ_{k=0}^∞ z^k / Γ(αk + β)
 *
 * This implementation uses multiple computational strategies:
 * - Taylor series for small |z|
 * - Asymptotic expansion for large |z|
 * - Integral representation for intermediate |z|
 * - Recursion relation for α > 1
 *
 * @param result       Output: computed Mittag-Leffler function value
 * @param alpha        Parameter α (alpha), typically α > 0
 * @param beta         Parameter β (beta)
 * @param z            Complex argument
 * @param prec         Precision in bits (use ML_DEFAULT_PRECISION for double precision)
 * @param taylor_radius Threshold for Taylor series region (typically ML_DEFAULT_TAYLOR_RADIUS)
 *
 * @example
 * @code
 * #include <mittleff/mittleff.h>
 *
 * int main(void) {
 *     acb_t result, alpha, beta, z;
 *     acb_init(result);
 *     acb_init(alpha);
 *     acb_init(beta);
 *     acb_init(z);
 *
 *     // Compute E_{0.5,1}(1.0)
 *     acb_set_d(alpha, 0.5);
 *     acb_set_d(beta, 1.0);
 *     acb_set_d(z, 1.0);
 *
 *     arb_t radius;
 *     arb_init(radius);
 *     arb_set_d(radius, ML_DEFAULT_TAYLOR_RADIUS);
 *
 *     ml_mittleff(result, alpha, beta, z, ML_DEFAULT_PRECISION, radius);
 *
 *     acb_printn(result, 15, 0);
 *     printf("\n");
 *
 *     acb_clear(result);
 *     acb_clear(alpha);
 *     acb_clear(beta);
 *     acb_clear(z);
 *     arb_clear(radius);
 *     return 0;
 * }
 * @endcode
 *
 * @note This function implements the algorithm from Garrappa (2015)
 * @since 0.1
 */
MITTLEFF_EXPORT void
ml_mittleff(acb_t result, const acb_t alpha, const acb_t beta, const acb_t z,
            slong prec, const arb_t taylor_radius);

/**
 * @brief Auxiliary Mittag-Leffler computation (internal)
 *
 * This is the core computation routine called by ml_mittleff after applying
 * the recursion relation. It assumes α ≤ 1 and dispatches to either:
 * - Asymptotic expansion (for large |z|)
 * - Integral representation (for intermediate |z|)
 *
 * @param result        Output: computed value
 * @param alpha         Parameter α (alpha), should satisfy α ≤ 1
 * @param beta          Parameter β (beta)
 * @param z             Complex argument
 * @param prec          Precision in bits
 * @param taylor_radius Threshold for Taylor series region
 *
 * @note This function is exposed for testing but typically called via ml_mittleff
 * @since 0.1
 */
MITTLEFF_EXPORT void
ml_mittleff_aux(acb_t result, const acb_t alpha, const acb_t beta, const acb_t z,
                slong prec, const arb_t taylor_radius);

#endif /* ML_MITTLEFF_H */
