#include <math.h>

#include "algorithm.h"
#include "integrate.h"
#include "log.h"
#include "new.h"

static num_t asymptotic_series (const num_t z, const num_t alpha, const num_t beta);

void
mittleff0 (double* res,
           const double a, const double b,
           const double x, const double y,
           const double acc)
{
    int k, k1, k2, kmax;
    double abs_z;
    num_t z, sum, tmp, fac1, fac2;

    /* Compute the maximum number of terms kmax to be taken into account for the
     * Taylor series, eq. (4.5) */
    abs_z = sqrt(x*x + y*y);
    k1 = (int) (ceil((2 - b)/a) + 1);
    k2 = (int) (ceil(log(acc * (1 - abs_z))/log(abs_z)) + 1);
    kmax = (k1 > k2) ? k1 : k2;

    /* Sum Taylor series */
    z = new(num), sum = new(num);
    tmp = new(num), fac1 = new(num), fac2 = new(num);
    num_set_d_d(z, x, y), num_set_d(sum, 0.0);
    for (k = 0; k <= kmax; k++)
    {
        /* fac1 <- z**k */
        num_pow_d(fac1, z, (double) k); 

        /* fac2 <- rgamma(alpha * k + beta) */
        num_set_d(fac2, a * k + b); 
        num_rgamma(fac2, fac2);

        /* partial sum = fac1 * fac2 */
        num_mul(tmp, fac1, fac2);
        
        num_add(sum, sum, tmp);
    }
    num_to_d_d(res, sum);
    delete(sum), delete(z), delete(tmp), delete(fac1), delete(fac2);
}

/* num_t */
/* mittleff1 (const num_t alpha, */
/*            const num_t beta, */
/*            const num_t z, */
/*            const num_t acc) */
/* { */
/*     log_debug("[%s] called", __func__); */
/*     num_t fac1, fac2; */

/*     fac1 = num_mul( */
/*         num_mul( */
/*             new(num, 1.0/num_to_d(alpha), 0.0), */
/*             num_pow(z, num_div(num_sub(new(num, 1.0, 0.0), beta), alpha))), */
/*         num_exp(num_pow(z, num_inv(alpha))));	    */
/*     fac2 = asymptotic_series(z, alpha, beta); */
		
/*     return num_add(fac1, fac2); */
/* } */

/* num_t */
/* mittleff2 (const num_t alpha, */
/*            const num_t beta, */
/*            const num_t z, */
/*            const num_t acc) */
/* { */
/*     return asymptotic_series(z, alpha, beta); */
/* } */

/* num_t */
/* mittleff3 (const num_t alpha, */
/*            const num_t beta, */
/*            const num_t z, */
/*            const num_t acc) */
/* { */
/*     /\* num_t a = PARAMS_GET_ALPHA(p); *\/ */
/*     /\* num_t b = PARAMS_GET_BETA(p); *\/ */

/*     /\* th = mp.arg(z**(1/alpha)) - np.pi *\/ */
/*     const double th = num_to_d( */
/*         num_sub( */
/*             num_arg(num_pow(z, num_inv(alpha))), */
/*             new(num, M_PI, 0.0))); */
		
/*     /\* c = th + 1j*th**2/6 - th**3/36 *\/ */
/*     const num_t c = new(num, th - th*th*th/36, th*th/6, 0.0); */
		
/*     /\* fac = (1/(2*alpha))*z**((1 - beta)/alpha)*mp.exp(z**(1/alpha))*mp.erfc(c*mp.sqrt(0.5*mp.fabs(z)**(1/alpha))) *\/ */
/*     num_t fac; */
/*     fac = num_inv(num_mul(new(num, 2.0, 0.0), alpha)); */
/*     fac = num_mul(fac, */
/*                   num_pow( */
/*                       z, */
/*                       num_div(num_sub(new(num, 1.0, 0.0), beta), alpha))); */
/*     fac = num_mul(fac, */
/*                   num_exp( */
/*                       num_pow( */
/*                           z, */
/*                           num_inv(alpha)))); */
/*     fac = num_mul(fac, */
/*                   num_erfc( */
/*                       num_mul( */
/*                           c, */
/*                           num_sqrt( */
/*                               num_mul( */
/*                                   new(num, 0.5, 0.0), */
/*                                   num_pow( */
/*                                       num_abs(z), */
/*                                       num_inv(alpha))))))); */
		
/*     return num_add(fac, asymptotic_series(z, alpha, beta)); */
/* } */

/* num_t */
/* mittleff4 (const num_t alpha, */
/*            const num_t beta, */
/*            const num_t z, */
/*            const num_t acc) */
/* { */
/*       /\* numeric_t a = PARAMS_GET_ALPHA(p); *\/ */
/*     /\* numeric_t b = PARAMS_GET_BETA(p); *\/ */

/*     /\* th = mp.arg(z**(1/alpha)) - np.pi *\/ */
/*     const double th = num_to_d( */
/*         num_sub( */
/*             num_arg( */
/*                 num_pow( */
/*                     z, */
/*                     num_inv(alpha))), */
/*             new(num, M_PI, 0.0))); */
		
/*     /\* c = th + 1j*th**2/6 - th**3/36 *\/ */
/*     const num_t c = new(num, -th + th*th*th/36.0, -th*th/6.0); */
		
/*     /\* fac = (1/(2*alpha))*z**((1 - beta)/alpha)*mp.exp(z**(1/alpha))*mp.erfc(c*mp.sqrt(0.5*mp.fabs(z)**(1/alpha))) *\/ */
/*     num_t fac; */
/*     fac = num_inv(num_mul(new(num, 2.0, 0.0), alpha)); */
/*     fac = num_mul( */
/*         fac, */
/*         num_pow( */
/*             z, */
/*             num_div( */
/*                 num_sub(new(num, 1.0, 0.0), beta), */
/*                 alpha))); */
/*     fac = num_mul( */
/*         fac, */
/*         num_exp( */
/*             num_pow( */
/*                 z, */
/*                 num_inv(alpha)))); */
/*     fac = num_mul( */
/*         fac, */
/*         num_erfc(num_mul( */
/*                      c, */
/*                      num_sqrt( */
/*                          num_mul( */
/*                              new(num, 0.5, 0.0), */
/*                              num_pow( */
/*                                  num_abs(z), num_inv(alpha))))))); */
		
/*     return num_add(fac, asymptotic_series(z, alpha, beta));	 */
/* } */

/* num_t */
/* mittleff5 (const num_t alpha, */
/*            const num_t beta, */
/*            const num_t z, */
/*            const num_t acc) */
/* { */
/*     /\* log_trace("[%s] (alpha=%g, beta=%g, z=(%+.5e, %+.5e), tol=%g)", *\/ */
/*     /\*           __func__,               *\/ */
/*     /\*           num_to_d(alpha), *\/ */
/*     /\*           num_to_d(beta), *\/ */
/*     /\*           num_real_d(z), *\/ */
/*     /\*           num_imag_d(z), *\/ */
/*     /\*           num_to_d(acc)); *\/ */

/*     const double _alpha = num_to_d(alpha); */
/*     const double _beta = num_to_d(beta); */
/*     const double _acc = num_to_d(acc); */
/*     const double eps = 0.5; */

/*     const num_t zero = new(num, 0.0, 0.0); */
/*     const num_t pi = new(num, M_PI, 0.0); */
/*     const num_t one = new(num, 1.0, 0.0); */
/*     const num_t two = new(num, 2.0, 0.0); */
/*     const num_t half = new(num, 0.5, 0.0); */

/*     log_trace("Before compute rmax"); */
    
/*     num_t res = zero; */
    
/*     /\* Compute r_max, equation (4.53) *\/ */
/*     num_t rmax = zero; */
/*     if (num_ge(beta, zero)) */
/*     { */
/*         rmax = num_max( */
/*             3, */
/*             one, */
/*             num_mul(two, num_abs(z)), */
/*             num_pow(num_neg(num_log(new(num, M_PI*eps/6.0, 0.0))), alpha)); */
/*     } */
/*     else */
/*         rmax = num_max( */
/*             3, */
/*             num_pow(num_add(num_abs(beta), one), alpha), */
/*             num_mul(two, num_abs(z)), */
/*             num_pow( */
/*                 num_neg( */
/*                     num_mul(two, */
/*                             num_log(new(num, M_PI*eps/(6.0*(fabs(_beta)+2.0)*pow(2*fabs(_beta), _beta)), 0.0)))), */
/*                 alpha)); */
/*     //log_trace("After compute rmax: %g %g", num_real_d(rmax), num_imag_d(rmax)); */
/*     num_t a_value = A(z, alpha, beta, zero); */
/*     num_t pi_alpha = num_mul(pi, alpha); */
/*     num_t integ_b = zero; */
/*     num_t integ_c = zero; */

/*     log_trace("Before compute integrals"); */
/*     if (num_le(beta, one)) /\* Equation (4.25) *\/ */
/*     { */
/*         integ_b = integrate_B(alpha, beta, z, pi_alpha, zero, rmax); */
/*         res = num_add(a_value, integ_b);         */
/*     } */
/*     else /\* Equation (4.26) *\/ */
/*     { */
/*         integ_b = integrate_B(alpha, beta, z, pi_alpha, half, rmax); */
/*         integ_c = integrate_C(alpha, beta, z, half, num_neg(pi_alpha), pi_alpha); */
/*         res = num_add(a_value, num_add(integ_b, integ_c));          */
/*     } */
/*     log_trace("After compute integrals"); */

/*     delete(pi_alpha); */
/*     delete(a_value); */
/*     delete(integ_b); */
/*     delete(integ_c); */

/*     log_trace("Before return function"); */

/*     return res; */
/* } */

/* num_t */
/* mittleff6 (const num_t alpha, */
/*            const num_t beta, */
/*            const num_t z, */
/*            const num_t acc) */
/* { */
/*     const double eps = 0.5; */
/*     const double _beta = num_to_d(beta); */
/*     const double _alpha = num_to_d(alpha); */

/*     const num_t zero = new(num, 0.0, 0.0); */
/*     const num_t pi = new(num, M_PI, 0.0); */
/*     const num_t one = new(num, 1.0, 0.0); */
/*     const num_t two = new(num, 2.0, 0.0); */
/*     const num_t half = new(num, 0.5, 0.0); */
    
/*     /\* Compute r_max, equation (4.54) *\/ */
/*     num_t rmax = zero; */
/*     if (num_ge(beta, zero)) */
/*     { */
/*         rmax = num_max( */
/*             3, */
/*             num_pow(two, alpha), */
/*             num_mul(two, num_abs(z)), */
/*             num_pow(num_neg(num_log(new(num, M_PI*eps*pow(2.0, _beta)/12.0, 0.0))), alpha)); */
/*     } */
/*     else // TODO Check with Hansjoerg whether the brackets mean "ceil" here */
/*         rmax = num_max( */
/*             3, */
/*             num_pow(num_ceil(num_mul(two, num_add(num_abs(beta), one))), alpha), */
/*             num_mul(two, num_abs(z)), */
/*             num_pow( */
/*                 num_ceil( */
/*                     num_neg( */
/*                         new(num, log((M_PI*pow(2.0,_beta)*eps)/(12*(fabs(_beta)+2)*pow(4*fabs(_beta), fabs(_beta)))), 0.0))), alpha)); */

/*     num_t res = num_from_d(0.0); */

/*     num_t integ_b = num_from_d(0.0); */
/*     num_t integ_c = num_from_d(0.0); */
/*     num_t v = new(num, 2.0 * M_PI * _alpha/3.0, 0.0); */
        
/*     if (num_le(beta, one)) /\* Equation (4.31) *\/ */
/*     { */
/*         res = integrate_B(alpha, beta, z, v, zero, rmax); */
/*     } */
/*     else /\* Equation (4.32) *\/ */
/*     { */
/*         integ_b = integrate_B(alpha, beta, z, v, half, rmax); */
/*         integ_c = integrate_C(alpha, beta, z, half, num_neg(v), v); */
/*         res = num_add(integ_b, integ_c);          */
/*     } */

/*     delete(zero); */
/*     delete(one); */
/*     delete(two); */
/*     delete(half); */
    
/*     delete(integ_b); */
/*     delete(integ_c); */
/*     delete(v); */

/*     return res; */
/* } */


/* /\* Return the rhs of equation (2.3) *\/ */
/* num_t */
/* asymptotic_series (const num_t z, const num_t alpha, const num_t beta) */
/* { */
/*     log_trace("[%s] called", __func__); */
/*     /\*kmax = int(mp.ceil((1/alpha)*mp.fabs(z)**(1/alpha)) + 1)*\/ */
/*     const int kmax = (int) num_to_d( */
/*         num_add( */
/*             num_ceil(num_mul(num_inv(alpha), */
/*                              num_pow(num_abs(z), */
/*                                      num_inv(alpha)))), */
/*             new(num, 1.0, 0.0))); */
/*     log_trace("[%s] summing %d terms", __func__, kmax); */

/*     /\* -sum([z**(-k) * mp.rgamma(beta - alpha*k) for k in range(1, kmax + 1)]) *\/ */
/*     num_t sum = new(num, 0.0, 0.0); */
/*     for (int k = 1; k <= kmax; k++) */
/*     { */
/*         sum = num_add(sum, */
/*                       num_mul( */
/*                           /\* z**(-k) *\/  */
/*                           num_pow(z, new(num, (double) -k, 0.0)), */
/*                           /\* mp.rgamma(beta - alpha*k) *\/ */
/*                           num_rgamma( */
/*                               num_sub(beta, */
/*                                       num_mul(alpha, */
/*                                               new(num, (double) k, 0.0)))))); */
/*     } */
/*     log_trace("[%s] Finished", __func__); */
/*     return num_mul(sum, new(num, -1.0, 0.0)); */
/* } */
