#ifndef QUAD_H
#define QUAD_H

#include <flint/arb.h>
#include <flint/acb.h>

/**
 * @file quad.h
 * @brief Adaptive Simpson quadrature integration
 *
 * This module implements adaptive Simpson's rule for numerical integration
 * as an alternative to FLINT's acb_calc_integrate. It's based on the
 * implementation from scheme_reference_implementation/mittleff/quad.scm.
 *
 * The adaptive Simpson method recursively subdivides intervals until
 * the tolerance is met or maximum recursion depth is reached.
 */

/**
 * Callback function type for real integration
 *
 * @param out    Output: function value at x
 * @param x      Input: point at which to evaluate
 * @param param  User-provided parameter pointer
 * @param prec   Precision in bits
 * @return 0 on success, non-zero on error
 */
typedef int (*quad_real_func)(arb_t out, const arb_t x, void* param, slong prec);

/**
 * Callback function type for complex integration
 *
 * @param out    Output: function value at x
 * @param x      Input: point at which to evaluate
 * @param param  User-provided parameter pointer
 * @param prec   Precision in bits
 * @return 0 on success, non-zero on error
 */
typedef int (*quad_complex_func)(acb_t out, const acb_t x, void* param, slong prec);

/**
 * Integrate a real function using adaptive Simpson's rule
 *
 * @param result   Output: integral value
 * @param func     Function to integrate
 * @param param    User parameter passed to func
 * @param a        Lower integration bound
 * @param b        Upper integration bound
 * @param tol      Absolute tolerance
 * @param max_depth Maximum recursion depth
 * @param prec     Precision in bits
 * @return 0 on success, non-zero on error
 */
int quad_real(arb_t result, quad_real_func func, void* param,
              const arb_t a, const arb_t b, const arb_t tol,
              slong max_depth, slong prec);

/**
 * Integrate a complex function using adaptive Simpson's rule
 *
 * Integrates real and imaginary parts separately using quad_real.
 *
 * @param result   Output: integral value
 * @param func     Function to integrate
 * @param param    User parameter passed to func
 * @param a        Lower integration bound (real)
 * @param b        Upper integration bound (real)
 * @param tol      Absolute tolerance
 * @param max_depth Maximum recursion depth
 * @param prec     Precision in bits
 * @return 0 on success, non-zero on error
 */
int quad_complex(acb_t result, quad_complex_func func, void* param,
                 const acb_t a, const acb_t b, const arb_t tol,
                 slong max_depth, slong prec);

#endif /* QUAD_H */
