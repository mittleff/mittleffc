#include <stdio.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <mittleff/mittleff.h>

int main(void) {
    acb_t a, b, z, result;
    arb_t radius, z_abs;
    slong prec = 53;
    
    acb_init(a);
    acb_init(b);
    acb_init(z);
    acb_init(result);
    arb_init(radius);
    arb_init(z_abs);
    
    // Test case: a=0.85, b=1.0, z=-45.0
    acb_set_d(a, 0.85);
    acb_set_d(b, 1.0);
    acb_set_d(z, -45.0);
    arb_set_d(radius, 0.95);
    
    acb_abs(z_abs, z, prec);
    
    printf("Testing E_{%.2f,%.2f}(%.2f)\n", 0.85, 1.0, -45.0);
    printf("|z| = ");
    arb_printn(z_abs, 10, 0);
    printf("\n");
    printf("Taylor radius = ");
    arb_printn(radius, 10, 0);
    printf("\n");
    
    // Compute
    ml_mittleff(result, a, b, z, prec, radius);
    
    printf("Result: ");
    acb_printn(result, 15, 0);
    printf("\n");
    
    // Check if result is special
    if (arf_is_special(arb_midref(acb_realref(result)))) {
        printf("Real part is NaN/Inf\n");
    }
    if (arf_is_special(arb_midref(acb_imagref(result)))) {
        printf("Imag part is NaN/Inf\n");
    }
    
    acb_clear(a);
    acb_clear(b);
    acb_clear(z);
    acb_clear(result);
    arb_clear(radius);
    arb_clear(z_abs);
    
    return 0;
}
