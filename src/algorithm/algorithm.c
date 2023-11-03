#include <math.h>

#include "algorithm.h"
#include "integrate.h"
#include "log.h"
#include "new.h"

static num_t A (const num_t z, const num_t alpha, const num_t beta, const num_t x);

num_t
mittleff0 (const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
    const double _alpha = num_to_double(alpha);
    const double _beta = num_to_double(beta);
    const double _acc = num_to_double(acc);
    
     log_trace("[%s] (alpha, beta, z, acc) = (%g, %g, %g%+gj, %g)",
              __func__,
              _alpha,
              _beta,
              num_real_d(z),
              num_imag_d(z),
              _acc);
     num_t zero = new(num, 0.0, 0.0);
     num_t one = new(num, 1.0, 0.0);
     num_t two = new(num, 2.0, 0.0);

     const double abs_z = num_to_double(num_abs(z));
     const int k1 = (int) (ceil((2 - _beta)/_alpha) + 1);
     const int k2 = (int) (ceil(log(_acc*(1-abs_z))/log(abs_z)) + 1);
     const int kmax = (k1 > k2) ? k1 : k2;

     num_t sum = zero;
     for (int k = 0; k <= kmax; k++)
     {
         num_t k_num = new(num, (double)k, 0.0);
         sum = num_add(
             sum,
             num_mul(num_pow(z, k_num),
                     num_rgamma(num_add(num_mul(alpha, k_num), beta))));
     }

     delete(zero); delete(one); delete(two);
     
    return sum;
}

num_t
mittleff5 (const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
    log_trace("[%s] (alpha=%g, beta=%g, z=(%+.5e, %+.5e), tol=%g)",
              __func__,              
              num_to_double(alpha),
              num_to_double(beta),
              num_real_d(z),
              num_imag_d(z),
              num_to_double(acc));

    const double _alpha = num_to_double(alpha);
    const double _beta = num_to_double(beta);
    const double _acc = num_to_double(acc);
    const double eps = 0.5;

    const num_t zero = new(num, 0.0, 0.0);
    const num_t one = new(num, 1.0, 0.0);
    const num_t two = new(num, 2.0, 0.0);
    const num_t half = new(num, 0.5, 0.0);
    
    num_t res = zero;
    
    /* Compute r_max, equation (4.53) */
    num_t rmax = zero;
    if (num_ge(beta, zero))
    {
        rmax = num_max(
            num_max(one,
                    num_mul(two, num_abs(z))),
            num_pow(num_negative(num_log(new(num, M_PI*eps/6.0, 0.0))), alpha));
    }
    else
        rmax = num_max(
            num_max(num_pow(num_add(num_abs(beta), one), alpha),
                    num_mul(two, num_abs(z))),
            num_pow(
                num_negative(
                    num_mul(two,
                            num_log(new(num, M_PI*eps/(6.0*(fabs(_beta)+2.0)*pow(2*fabs(_beta), _beta)), 0.0)))),
                alpha));

    num_t a_value = A(z, alpha, beta, zero);
    num_t integ_b = integrate_B();
    num_t integ_c = zero;
    
    if (num_le(beta, one)) /* Equation (4.25) */
    {
        res = num_add(a_value, integ_b);        
    }
    else /* Equation (4.26) */
    {       
        res = num_add(a_value, num_add(integ_b, integ_c));         
    }
    
    return res;
}

num_t
mittleff6 (const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
    const num_t zero = new(num, 0.0, 0.0);
    /* Compute r_max, equation (4.54) */
    num_t rmax = zero;
    return new(num, 0.5, 0.0);
}

/****************************************************************/

static num_t
A (const num_t z, const num_t alpha, const num_t beta, const num_t x)
{
    const double _alpha = num_to_double(alpha);
    const double _beta = num_to_double(beta);
    
    return num_mul(
        num_inverse(alpha),
        num_mul(
            num_pow(z, new(num, (1 - _beta)/_alpha, 0.0)),
            num_exp(num_mul(num_pow(z, num_inverse(alpha)), num_cos(num_div(x, alpha))))));
}
