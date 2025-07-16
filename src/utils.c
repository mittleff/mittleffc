#include <flint/arb.h>

#include "include/utils.h"

static void
accuracy_from_prec_flint (arb_t result, slong prec)
{
    arb_t tmp, ten_const;

    /* Initialize temporary variables */
    arb_init(tmp);
    arb_init(ten_const);

    arb_one(tmp);
    arb_set_si(ten_const, 10);

    /* Compute 2^(-prec) */
    arb_mul_2exp_si(tmp, tmp, -prec); 
    //arb_printd(tmp, 10); printf("\n"); 
    
    /* Compute log10(2^(-prec)) */
    arb_log_base_ui(tmp, tmp, 10, prec);
    //arb_printd(tmp, 10);
    
    /* Compute abs(log10(...)) */
    arb_abs(tmp, tmp);
    
    /* Compute ceiling of the absolute value */
    arb_ceil(tmp, tmp, prec);
    
    /* Compute 10^(- ceiling(...)) */
    arb_neg(tmp, tmp);
    arb_pow(result, ten_const, tmp, prec);
    
    /* Clear temporary variables */
    arb_clear(tmp);
    arb_clear(ten_const);
}

/* Wrapper for shared library */
double
accuracy_from_prec (unsigned long int prec)
{
    arb_t result_arb;
    double result;

    /* Initialize temporary variables */
    arb_init(result_arb);

    accuracy_from_prec_flint(result_arb, prec);
    result = arf_get_d(arb_midref(result_arb), ARF_RND_NEAR);
    
    /* Clear temporary variables */
    arb_clear(result_arb);

    return result;    
}

