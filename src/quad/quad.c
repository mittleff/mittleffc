/**
 * @file quad.c
 * @brief Adaptive Simpson quadrature implementation
 *
 * Based on the Scheme implementation from:
 * scheme_reference_implementation/mittleff/quad.scm
 * and Rosetta Code:
 * https://rosettacode.org/wiki/Numerical_integration/Adaptive_Simpson%27s_method
 */

#include "quad.h"
#include <stdio.h>

/**
 * Compute Simpson's rule estimate for interval [a, b]
 *
 * Returns: m (midpoint), fm (function at midpoint), and estimate
 */
static int simpson_step(arb_t m, arb_t fm, arb_t estimate,
                        const arb_t a, const arb_t fa,
                        const arb_t b, const arb_t fb,
                        quad_real_func func, void* param, slong prec)
{
    arb_t h, temp;
    int ret;

    arb_init(h);
    arb_init(temp);

    /* m = (a + b) / 2 */
    arb_add(m, a, b, prec);
    arb_mul_2exp_si(m, m, -1);

    /* fm = f(m) */
    ret = func(fm, m, param, prec);
    if (ret != 0) {
        arb_clear(h);
        arb_clear(temp);
        return ret;
    }

    /* h = b - a */
    arb_sub(h, b, a, prec);

    /* estimate = (h/6) * (fa + 4*fm + fb) */
    arb_mul_ui(temp, fm, 4, prec);
    arb_add(temp, temp, fa, prec);
    arb_add(temp, temp, fb, prec);
    arb_mul(temp, temp, h, prec);
    arb_div_ui(estimate, temp, 6, prec);

    arb_clear(h);
    arb_clear(temp);
    return 0;
}

/**
 * Recursive adaptive Simpson quadrature
 *
 * @param result     Output: integral estimate
 * @param a, fa      Left endpoint and function value
 * @param b, fb      Right endpoint and function value
 * @param tol        Current tolerance
 * @param whole      Simpson estimate for whole interval
 * @param m, fm      Midpoint and function value
 * @param depth      Remaining recursion depth
 * @param func       Function to integrate
 * @param param      User parameter
 * @param prec       Precision in bits
 */
static int quad_asr_recursive(arb_t result,
                               const arb_t a, const arb_t fa,
                               const arb_t b, const arb_t fb,
                               const arb_t tol, const arb_t whole,
                               const arb_t m, const arb_t fm,
                               slong depth,
                               quad_real_func func, void* param, slong prec)
{
    arb_t lm, flm, left, rm, frm, right;
    arb_t delta, abs_delta, tol_half, threshold;
    arb_t left_result, right_result, correction;
    int ret;

    /* Check termination conditions */
    if (depth <= 0) {
        /* Max depth reached, return current estimate with correction */
        arb_init(left);
        arb_init(right);
        arb_init(delta);
        arb_init(correction);

        ret = simpson_step(lm, flm, left, a, fa, m, fm, func, param, prec);
        if (ret != 0) goto cleanup_depth;

        ret = simpson_step(rm, frm, right, m, fm, b, fb, func, param, prec);
        if (ret != 0) goto cleanup_depth;

        arb_add(delta, left, right, prec);
        arb_sub(delta, delta, whole, prec);
        arb_div_ui(correction, delta, 15, prec);
        arb_add(result, left, right, prec);
        arb_add(result, result, correction, prec);

cleanup_depth:
        arb_clear(left);
        arb_clear(right);
        arb_clear(delta);
        arb_clear(correction);
        return ret;
    }

    arb_init(lm);
    arb_init(flm);
    arb_init(left);
    arb_init(rm);
    arb_init(frm);
    arb_init(right);
    arb_init(delta);
    arb_init(abs_delta);
    arb_init(tol_half);
    arb_init(threshold);
    arb_init(left_result);
    arb_init(right_result);
    arb_init(correction);

    /* Compute Simpson estimate for left half [a, m] */
    ret = simpson_step(lm, flm, left, a, fa, m, fm, func, param, prec);
    if (ret != 0) goto cleanup;

    /* Compute Simpson estimate for right half [m, b] */
    ret = simpson_step(rm, frm, right, m, fm, b, fb, func, param, prec);
    if (ret != 0) goto cleanup;

    /* delta = (left + right) - whole */
    arb_add(delta, left, right, prec);
    arb_sub(delta, delta, whole, prec);

    /* tol_half = tol / 2 */
    arb_mul_2exp_si(tol_half, tol, -1);

    /* Check convergence: |delta| <= 15 * tol or tol_half == tol (underflow) */
    arb_abs(abs_delta, delta);
    arb_mul_ui(threshold, tol, 15, prec);

    if (arb_eq(tol_half, tol) || arb_le(abs_delta, threshold)) {
        /* Converged: return left + right + delta/15 */
        arb_div_ui(correction, delta, 15, prec);
        arb_add(result, left, right, prec);
        arb_add(result, result, correction, prec);
    } else {
        /* Not converged: recurse on both halves */
        ret = quad_asr_recursive(left_result, a, fa, m, fm, tol_half, left,
                                 lm, flm, depth - 1, func, param, prec);
        if (ret != 0) goto cleanup;

        ret = quad_asr_recursive(right_result, m, fm, b, fb, tol_half, right,
                                 rm, frm, depth - 1, func, param, prec);
        if (ret != 0) goto cleanup;

        arb_add(result, left_result, right_result, prec);
    }

cleanup:
    arb_clear(lm);
    arb_clear(flm);
    arb_clear(left);
    arb_clear(rm);
    arb_clear(frm);
    arb_clear(right);
    arb_clear(delta);
    arb_clear(abs_delta);
    arb_clear(tol_half);
    arb_clear(threshold);
    arb_clear(left_result);
    arb_clear(right_result);
    arb_clear(correction);

    return ret;
}

/**
 * Integrate a real function using adaptive Simpson's rule
 */
int quad_real(arb_t result, quad_real_func func, void* param,
              const arb_t a, const arb_t b, const arb_t tol,
              slong max_depth, slong prec)
{
    arb_t fa, fb, m, fm, whole;
    int ret;

    arb_init(fa);
    arb_init(fb);
    arb_init(m);
    arb_init(fm);
    arb_init(whole);

    /* Evaluate at endpoints */
    ret = func(fa, a, param, prec);
    if (ret != 0) goto cleanup;

    ret = func(fb, b, param, prec);
    if (ret != 0) goto cleanup;

    /* Compute initial Simpson estimate */
    ret = simpson_step(m, fm, whole, a, fa, b, fb, func, param, prec);
    if (ret != 0) goto cleanup;

    /* Start recursive adaptive integration */
    ret = quad_asr_recursive(result, a, fa, b, fb, tol, whole, m, fm,
                             max_depth, func, param, prec);

cleanup:
    arb_clear(fa);
    arb_clear(fb);
    arb_clear(m);
    arb_clear(fm);
    arb_clear(whole);

    return ret;
}

/**
 * Helper structure for complex integration
 */
typedef struct {
    quad_complex_func func;
    void* user_param;
    int extract_real; /* 1 for real part, 0 for imaginary part */
} complex_wrapper_param_t;

/**
 * Wrapper to extract real part of complex function
 */
static int complex_real_wrapper(arb_t out, const arb_t x, void* param, slong prec)
{
    complex_wrapper_param_t* p = (complex_wrapper_param_t*)param;
    acb_t x_complex, result;
    int ret;

    acb_init(x_complex);
    acb_init(result);

    acb_set_arb(x_complex, x);
    ret = p->func(result, x_complex, p->user_param, prec);

    if (ret == 0) {
        if (p->extract_real) {
            acb_get_real(out, result);
        } else {
            acb_get_imag(out, result);
        }
    }

    acb_clear(x_complex);
    acb_clear(result);

    return ret;
}

/**
 * Integrate a complex function using adaptive Simpson's rule
 */
int quad_complex(acb_t result, quad_complex_func func, void* param,
                 const acb_t a, const acb_t b, const arb_t tol,
                 slong max_depth, slong prec)
{
    arb_t a_real, b_real, real_part, imag_part;
    complex_wrapper_param_t wrapper_param;
    int ret;

    arb_init(a_real);
    arb_init(b_real);
    arb_init(real_part);
    arb_init(imag_part);

    /* Extract real parts of bounds */
    acb_get_real(a_real, a);
    acb_get_real(b_real, b);

    /* Set up wrapper parameters */
    wrapper_param.func = func;
    wrapper_param.user_param = param;

    /* Integrate real part */
    wrapper_param.extract_real = 1;
    ret = quad_real(real_part, complex_real_wrapper, &wrapper_param,
                    a_real, b_real, tol, max_depth, prec);
    if (ret != 0) goto cleanup;

    /* Integrate imaginary part */
    wrapper_param.extract_real = 0;
    ret = quad_real(imag_part, complex_real_wrapper, &wrapper_param,
                    a_real, b_real, tol, max_depth, prec);
    if (ret != 0) goto cleanup;

    /* Combine into complex result */
    acb_set_arb_arb(result, real_part, imag_part);

cleanup:
    arb_clear(a_real);
    arb_clear(b_real);
    arb_clear(real_part);
    arb_clear(imag_part);

    return ret;
}
