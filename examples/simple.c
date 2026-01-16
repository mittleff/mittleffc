/**
 * @file simple.c
 * @brief Simple example demonstrating mittleff library usage
 *
 * This example shows how to compute the Mittag-Leffler function E_{α,β}(z).
 * Compile with: gcc -o simple simple.c $(pkg-config --cflags --libs mittleff)
 */

#include <stdio.h>
#include <mittleff/mittleff.h>
#include <flint/acb.h>
#include <flint/arb.h>

int main(void)
{
    printf("=== Mittag-Leffler Function Example ===\n\n");

    /* Initialize variables */
    acb_t result, alpha, beta, z;
    arb_t radius;

    acb_init(result);
    acb_init(alpha);
    acb_init(beta);
    acb_init(z);
    arb_init(radius);

    /* Set parameters: E_{0.5,1}(1.0) */
    acb_set_d(alpha, 0.5);
    acb_set_d(beta, 1.0);
    acb_set_d(z, 1.0);
    arb_set_d(radius, ML_DEFAULT_TAYLOR_RADIUS);

    printf("Computing E_{%.1f,%.1f}(%.1f)...\n",
           arf_get_d(arb_midref(acb_realref(alpha)), ARF_RND_NEAR),
           arf_get_d(arb_midref(acb_realref(beta)), ARF_RND_NEAR),
           arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR));

    /* Compute Mittag-Leffler function */
    ml_mittleff(result, alpha, beta, z, ML_DEFAULT_PRECISION, radius);

    printf("Result: ");
    acb_printn(result, 15, 0);
    printf("\n");

    /* Clean up */
    acb_clear(result);
    acb_clear(alpha);
    acb_clear(beta);
    acb_clear(z);
    arb_clear(radius);

    printf("\nExample completed successfully!\n");
    return 0;
}
