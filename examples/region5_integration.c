/**
 * @file region5_integration.c
 * @brief Example demonstrating integral representation (region 5) computation
 *
 * This example shows how to compute the Mittag-Leffler function using the
 * integral representation method, which is used for intermediate |z| values
 * in regions 5-6.
 *
 * Test case: alpha=0.5, beta=1.0, z=4.08373780+2.53485316i
 * Expected result: -15817.8109396066 + 54393.0514682909i
 *
 * This test case was particularly challenging because it involves beta=1.0,
 * which creates 0^0 expressions in the integrand. The implementation uses
 * adaptive Simpson integration to avoid numerical issues.
 */

#include <stdio.h>
#include <stdlib.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include "algorithm.h"

int main(void)
{
    acb_t a, b, z, result;
    slong prec = 53;

    acb_init(a);
    acb_init(b);
    acb_init(z);
    acb_init(result);

    /* Set parameters for region 5 test case */
    double alpha = 0.5;
    double beta = 1.0;
    double re_z = 4.08373780000000e+0;
    double im_z = 2.53485316000000e+0;

    acb_set_d_d(a, alpha, 0.0);
    acb_set_d_d(b, beta, 0.0);
    acb_set_d_d(z, re_z, im_z);

    printf("Computing Mittag-Leffler function using integral representation (region 5)\n");
    printf("Parameters:\n");
    printf("  alpha = %.6f\n", alpha);
    printf("  beta  = %.6f\n", beta);
    printf("  z     = %.6f + %.6fi\n", re_z, im_z);
    printf("\n");

    /* Compute using integral representation (region 5) */
    int region = 5;
    ml_integral_rep(result, a, b, z, region, prec);

    printf("Result:\n");
    printf("  E_%.1f,%.1f(z) = ", alpha, beta);
    acb_fprintd(stdout, result, 8);
    printf("\n\n");

    printf("Expected: -15817.811 + 54393.051i\n");
    printf("\n");

    /* Additional diagnostic information */
    printf("Note: This computation uses adaptive Simpson integration instead of\n");
    printf("FLINT's acb_calc_integrate to avoid special value errors that occur\n");
    printf("when beta=1.0 (which creates 0^0 expressions in the integrand).\n");

    acb_clear(a);
    acb_clear(b);
    acb_clear(z);
    acb_clear(result);

    return 0;
}
