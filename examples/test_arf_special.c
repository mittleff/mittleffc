#include <stdio.h>
#include <flint/arb.h>
#include <flint/arf.h>

int main(void) {
    arb_t zero, nan, inf, normal;
    arb_init(zero);
    arb_init(nan);
    arb_init(inf);
    arb_init(normal);
    
    arb_zero(zero);
    arb_set_d(normal, 3.14);
    
    // Create NaN and Inf
    arf_t nan_arf, inf_arf;
    arf_init(nan_arf);
    arf_init(inf_arf);
    arf_nan(nan_arf);
    arf_pos_inf(inf_arf);
    arb_set_arf(nan, nan_arf);
    arb_set_arf(inf, inf_arf);
    
    printf("Zero:\n");
    printf("  is_special: %d\n", arf_is_special(arb_midref(zero)));
    printf("  is_finite:  %d\n", arf_is_finite(arb_midref(zero)));
    printf("  is_nan:     %d\n", arf_is_nan(arb_midref(zero)));
    printf("  is_inf:     %d\n", arf_is_inf(arb_midref(zero)));
    printf("  is_zero:    %d\n", arb_is_zero(zero));
    
    printf("\nNormal (3.14):\n");
    printf("  is_special: %d\n", arf_is_special(arb_midref(normal)));
    printf("  is_finite:  %d\n", arf_is_finite(arb_midref(normal)));
    printf("  is_nan:     %d\n", arf_is_nan(arb_midref(normal)));
    printf("  is_inf:     %d\n", arf_is_inf(arb_midref(normal)));
    
    printf("\nNaN:\n");
    printf("  is_special: %d\n", arf_is_special(arb_midref(nan)));
    printf("  is_finite:  %d\n", arf_is_finite(arb_midref(nan)));
    printf("  is_nan:     %d\n", arf_is_nan(arb_midref(nan)));
    printf("  is_inf:     %d\n", arf_is_inf(arb_midref(nan)));
    
    printf("\nInf:\n");
    printf("  is_special: %d\n", arf_is_special(arb_midref(inf)));
    printf("  is_finite:  %d\n", arf_is_finite(arb_midref(inf)));
    printf("  is_nan:     %d\n", arf_is_nan(arb_midref(inf)));
    printf("  is_inf:     %d\n", arf_is_inf(arb_midref(inf)));
    
    arb_clear(zero);
    arb_clear(nan);
    arb_clear(inf);
    arb_clear(normal);
    arf_clear(nan_arf);
    arf_clear(inf_arf);
    
    return 0;
}
