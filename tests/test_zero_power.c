#include <flint/arb.h>
#include <flint/acb.h>
#include <stdio.h>

int main() {
    acb_t r, exponent, result;
    slong prec = 48;

    acb_init(r);
    acb_init(exponent);
    acb_init(result);

    // Test case: r=0, exponent = (1-beta)/alpha = (1-1.0)/0.5 = 0
    acb_zero(r);
    acb_zero(exponent);  // 0^0

    printf("Computing 0^0...\n");
    acb_pow(result, r, exponent, prec);

    printf("Result: ");
    acb_printd(result, 10);
    printf("\n");

    // Check if result contains special values
    int has_special = arf_is_special(arb_midref(acb_realref(result))) ||
                      arf_is_special(arb_midref(acb_imagref(result)));
    printf("Has special value: %d\n", has_special);

    // Also test with small positive exponent
    acb_zero(r);
    acb_set_d(exponent, 0.001);

    printf("\nComputing 0^0.001...\n");
    acb_pow(result, r, exponent, prec);

    printf("Result: ");
    acb_printd(result, 10);
    printf("\n");

    acb_clear(r);
    acb_clear(exponent);
    acb_clear(result);

    return 0;
}
