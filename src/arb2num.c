#include "arb2num.h"

#include <arf.h>
#include "new.h"

// Converts an arb_t number to double.
double
arbtod (const arb_t x)
{
    return arf_get_d(arb_midref(x), ARF_RND_NEAR);
}

num_t
num_from_acb (const acb_t x)
{
    arb_t re, im;
    
    arb_init(re); arb_init(im);
    acb_get_real(re, x);
    acb_get_imag(im, x);
    return new(num, arbtod(re), arbtod(im));
}

num_t
num_from_arb (const arb_t x)
{
    return new(num, arbtod(x), 0.0);
}
