#include <stdio.h>
#include <math.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include "../src/mittleff/mittleff.h"

static double
compute_relative_error(const acb_t computed, const acb_t expected, slong prec)
{
    arb_t diff_re, diff_im, exp_re, exp_im;
    arb_t rel_err_re, rel_err_im, max_rel_err;
    double result;

    /* Check for non-finite values (NaN/Inf) in the ARF midpoints */
    if (!arf_is_finite(arb_midref(acb_realref(computed))) ||
        !arf_is_finite(arb_midref(acb_imagref(computed))) ||
        !arf_is_finite(arb_midref(acb_realref(expected))) ||
        !arf_is_finite(arb_midref(acb_imagref(expected)))) {
        return INFINITY;
    }

    arb_init(diff_re);
    arb_init(diff_im);
    arb_init(exp_re);
    arb_init(exp_im);
    arb_init(rel_err_re);
    arb_init(rel_err_im);
    arb_init(max_rel_err);

    /* Get components */
    arb_set(diff_re, acb_realref(computed));
    arb_set(diff_im, acb_imagref(computed));
    arb_set(exp_re, acb_realref(expected));
    arb_set(exp_im, acb_imagref(expected));

    printf("DEBUG: Before subtraction:\n");
    printf("  computed real: "); arb_printn(diff_re, 20, 0); printf("\n");
    printf("  computed imag: "); arb_printn(diff_im, 20, 0); printf("\n");
    printf("  expected real: "); arb_printn(exp_re, 20, 0); printf("\n");
    printf("  expected imag: "); arb_printn(exp_im, 20, 0); printf("\n");

    /* Compute absolute differences */
    arb_sub(diff_re, diff_re, exp_re, prec);
    arb_sub(diff_im, diff_im, exp_im, prec);
    arb_abs(diff_re, diff_re);
    arb_abs(diff_im, diff_im);

    printf("DEBUG: After subtraction (absolute differences):\n");
    printf("  diff_re: "); arb_printn(diff_re, 20, 0); printf("\n");
    printf("  diff_im: "); arb_printn(diff_im, 20, 0); printf("\n");

    /* Compute absolute values of expected */
    arb_abs(exp_re, exp_re);
    arb_abs(exp_im, exp_im);

    printf("DEBUG: Absolute expected values:\n");
    printf("  |exp_re|: "); arb_printn(exp_re, 20, 0); printf("\n");
    printf("  |exp_im|: "); arb_printn(exp_im, 20, 0); printf("\n");
    printf("  exp_re is_zero: %d\n", arb_is_zero(exp_re));
    printf("  exp_im is_zero: %d\n", arb_is_zero(exp_im));

    /* Compute relative errors (handle zero expected values properly) */
    if (arb_is_zero(exp_re) && arb_is_zero(exp_im)) {
        /* Both components zero: use absolute error */
        printf("DEBUG: Both components zero, using absolute error\n");
        arb_add(max_rel_err, diff_re, diff_im, prec);
    } else {
        /* At least one component non-zero: compute relative errors */
        if (arb_is_zero(exp_re)) {
            printf("DEBUG: exp_re is zero, using absolute error for real\n");
            arb_set(rel_err_re, diff_re);  /* Use absolute error */
        } else {
            printf("DEBUG: exp_re is non-zero, computing relative error for real\n");
            arb_div(rel_err_re, diff_re, exp_re, prec);
        }

        if (arb_is_zero(exp_im)) {
            printf("DEBUG: exp_im is zero, using absolute error for imag\n");
            arb_set(rel_err_im, diff_im);  /* Use absolute error */
        } else {
            printf("DEBUG: exp_im is non-zero, computing relative error for imag\n");
            arb_div(rel_err_im, diff_im, exp_im, prec);
        }

        printf("DEBUG: Relative errors:\n");
        printf("  rel_err_re: "); arb_printn(rel_err_re, 20, 0); printf("\n");
        printf("  rel_err_im: "); arb_printn(rel_err_im, 20, 0); printf("\n");

        /* Take maximum of the two relative errors */
        arb_max(max_rel_err, rel_err_re, rel_err_im, prec);
    }

    printf("DEBUG: max_rel_err: "); arb_printn(max_rel_err, 20, 0); printf("\n");
    printf("DEBUG: is_special: %d\n", arf_is_special(arb_midref(max_rel_err)));

    /* Check if the result is finite before converting to double */
    if (arf_is_special(arb_midref(max_rel_err))) {
        printf("DEBUG: max_rel_err is special, returning INFINITY\n");
        result = INFINITY;
    } else {
        result = arf_get_d(arb_midref(max_rel_err), ARF_RND_NEAR);
        printf("DEBUG: Converted to double: %.17e\n", result);
    }

    arb_clear(diff_re);
    arb_clear(diff_im);
    arb_clear(exp_re);
    arb_clear(exp_im);
    arb_clear(rel_err_re);
    arb_clear(rel_err_im);
    arb_clear(max_rel_err);

    return result;
}

int main(void) {
    acb_t a, b, z, result, expected;
    arb_t radius;
    slong prec = 53;
    double rel_error;
    double tolerance = 2e-6;

    acb_init(a);
    acb_init(b);
    acb_init(z);
    acb_init(result);
    acb_init(expected);
    arb_init(radius);

    /* Test case: α=0.85, β=1.0, z=0 */
    acb_set_d_d(a, 0.85, 0.0);
    acb_set_d_d(b, 1.0, 0.0);
    acb_set_d_d(z, 0.0, 0.0);
    acb_set_d_d(expected, 1.0, 0.0);
    arb_set_d(radius, 0.95);

    printf("Testing E_{0.85,1.0}(0) - Expected: 1.0\n\n");

    ml_mittleff(result, a, b, z, prec, radius);

    printf("\nResult: ");
    acb_printn(result, 20, 0);
    printf("\n\n");

    rel_error = compute_relative_error(result, expected, prec);

    printf("\n=== FINAL RESULT ===\n");
    printf("Relative error: %.17e\n", rel_error);
    printf("Tolerance: %.17e\n", tolerance);
    printf("Pass: %s\n", (rel_error <= tolerance) ? "YES" : "NO");
    printf("isinf(rel_error): %d\n", isinf(rel_error));

    acb_clear(a);
    acb_clear(b);
    acb_clear(z);
    acb_clear(result);
    acb_clear(expected);
    arb_clear(radius);

    return 0;
}
