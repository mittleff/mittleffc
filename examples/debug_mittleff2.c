#include <stdio.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <mittleff/mittleff.h>

int main(void) {
    acb_t a, b, z, result;
    arb_t radius;
    slong prec = 53;
    
    acb_init(a);
    acb_init(b);
    acb_init(z);
    acb_init(result);
    arb_init(radius);
    
    // Test case: a=0.85, b=1.0, z=-45.0
    acb_set_d(a, 0.85);
    acb_set_d(b, 1.0);
    acb_set_d(z, -45.0);
    arb_set_d(radius, 0.95);
    
    printf("Testing E_{%.2f,%.2f}(%.2f)\n", 0.85, 1.0, -45.0);
    
    // Compute
    ml_mittleff(result, a, b, z, prec, radius);
    
    printf("Full result: ");
    acb_printd(result, 15);
    printf("\n");
    
    printf("Real part: ");
    arb_printd(acb_realref(result), 15);
    printf("\n");
    
    printf("Imag part: ");
    arb_printd(acb_imagref(result), 15);
    printf("\n");
    
    // Check if parts are special
    printf("Real is_special: %d\n", arf_is_special(arb_midref(acb_realref(result))));
    printf("Imag is_special: %d\n", arf_is_special(arb_midref(acb_imagref(result))));
    
    // Check if parts are zero
    printf("Imag is_zero: %d\n", arb_is_zero(acb_imagref(result)));
    
    acb_clear(a);
    acb_clear(b);
    acb_clear(z);
    acb_clear(result);
    arb_clear(radius);
    
    return 0;
}
