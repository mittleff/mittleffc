#include <stdio.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <mittleff/mittleff.h>

int main(void) {
    acb_t a, b, z, result;
    arb_t radius, alpha_real, one;
    slong prec = 53;
    
    acb_init(a);
    acb_init(b);
    acb_init(z);
    acb_init(result);
    arb_init(radius);
    arb_init(alpha_real);
    arb_init(one);
    arb_set_d(radius, 0.95);
    arb_one(one);
    
    // Test α > 1 case (should trigger recursion)
    printf("Testing recursion with α > 1:\n\n");
    
    // Test 1: α=1.25, β=1, z=0
    printf("Test 1: E_{1.25,1}(0) [Expected: 1.0]\n");
    acb_set_d(a, 1.25);
    acb_set_d(b, 1.0);
    acb_zero(z);
    
    acb_get_real(alpha_real, a);
    printf("  Alpha = ");
    arb_printn(alpha_real, 10, 0);
    printf("\n  Alpha > 1? %d\n", arb_gt(alpha_real, one));
    
    ml_mittleff(result, a, b, z, prec, radius);
    printf("  Result: ");
    acb_printn(result, 15, 0);
    printf("\n\n");
    
    // Test 2: α=1.25, β=1, z=1
    printf("Test 2: E_{1.25,1}(1) [Should NOT be 0.0138...]\n");
    acb_set_d(z, 1.0);
    ml_mittleff(result, a, b, z, prec, radius);
    printf("  Result: ");
    acb_printn(result, 15, 0);
    printf("\n\n");
    
    // Test 3: α=2.3, β=0.7, z=-200
    printf("Test 3: E_{2.3,0.7}(-200) [Expected: -9.399...]\n");
    acb_set_d(a, 2.3);
    acb_set_d(b, 0.7);
    acb_set_d(z, -200.0);
    ml_mittleff(result, a, b, z, prec, radius);
    printf("  Result: ");
    acb_printn(result, 15, 0);
    printf("\n");
    
    acb_clear(a);
    acb_clear(b);
    acb_clear(z);
    acb_clear(result);
    arb_clear(radius);
    arb_clear(alpha_real);
    arb_clear(one);
    
    return 0;
}
