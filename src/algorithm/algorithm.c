#include <math.h>

#include "algorithm.h"
#include "integrate.h"
#include "log.h"
#include "new.h"

num_t
mittleff0 (const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
    const double _alpha = num_to_double(alpha);
    const double _beta = num_to_double(beta);
    const double _acc = num_to_double(acc);
    
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
    const num_t pi = new(num, M_PI, 0.0);
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
    num_t pi_alpha = num_mul(pi, alpha);
    num_t integ_b;
    num_t integ_c;
    
    if (num_le(beta, one)) /* Equation (4.25) */
    {
        integ_b = integrate_B(alpha, beta, z, pi_alpha, zero, rmax);
        res = num_add(a_value, integ_b);        
    }
    else /* Equation (4.26) */
    {
        integ_b = integrate_B(alpha, beta, z, pi_alpha, half, rmax);
        integ_c = integrate_C(alpha, beta, z, half, num_negative(pi_alpha), pi_alpha);
        res = num_add(a_value, num_add(integ_b, integ_c));         
    }

    delete(pi_alpha); delete(a_value); delete(integ_b); delete(integ_c);

    return res;
}

num_t
mittleff6 (const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
    const double eps = 0.5;
    const double _beta = num_to_double(beta);
    const double _alpha = num_to_double(alpha);

    const num_t zero = new(num, 0.0, 0.0);
    const num_t pi = new(num, M_PI, 0.0);
    const num_t one = new(num, 1.0, 0.0);
    const num_t two = new(num, 2.0, 0.0);
    const num_t half = new(num, 0.5, 0.0);
    
    /* Compute r_max, equation (4.54) */
    num_t rmax = zero;
    if (num_ge(beta, zero))
    {
        rmax = num_max(
            num_max(num_pow(two, alpha),
                    num_mul(two, num_abs(z))),
            num_pow(num_negative(num_log(new(num, M_PI*eps*pow(2.0, _beta)/12.0, 0.0))), alpha));
    }
    else // TODO Check with Hansjoerg whether the brackets mean "ceil" here
        rmax = num_max(
            num_max(
                num_pow(num_ceil(num_mul(two, num_add(num_abs(beta), 1))), alpha),
                num_mul(two, num_abs(z))),
            num_pow(
                num_ceil(
                    num_negative(
                        new(num, log((M_PI*pow(2.0,_beta)*eps)/(12*(fabs(_beta)+2)*pow(4*fabs(_beta), fabs(_beta)))), 0.0))), alpha));

    num_t res = zero;

    num_t integ_b;
    num_t integ_c;
    num_t v = new(num, 2.0 * M_PI * _alpha/3.0, 0.0);
        
    if (num_le(beta, one)) /* Equation (4.31) */
    {
        res = integrate_B(alpha, beta, z, v, zero, rmax);
    }
    else /* Equation (4.32) */
    {
        integ_b = integrate_B(alpha, beta, z, v, half, rmax);
        integ_c = integrate_C(alpha, beta, z, half, num_negative(v), v);
        res = num_add(integ_b, integ_c);         
    }

    delete(integ_b); delete(integ_c); delete(v);

    return res;
}


