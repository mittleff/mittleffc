#include "algorithm.h"

#include "integrate.h"
#include "new.h"

#include "utils.h"

#include <math.h>
#include <stdbool.h>

#ifdef DEBUG
#include "log.h"
#endif

/* #define MAX2(a,b) ((a>b)?(a):(b)) */
/* #define MAX3(a,b,c) MAX2(a,MAX2(b,c)) */

static void
asymptotic_series (num_t res, const num_t z, const num_t alpha, const num_t beta);

void
mittleff0 (num_t res,
           const num_t alpha, const num_t beta,
           const num_t z,
           const num_t acc)
{
    /* log_trace("[%s] alpha=%g, beta=%g, z=%g%+g, acc=%g", __func__, */
    /*           num_to_d(alpha), num_to_d(beta), num_real_d(z), num_imag_d(z), num_to_d(acc)); */
    
    int k, kmax;
    num_t absz, sum, tmp, fac1, fac2, k1, k2;

    absz = new(num);
    k1 = new(num), k2 = new(num);
    fac1 = new(num), fac2 = new(num);
    tmp = new(num);
    sum = new(num);

    /* Compute the maximum number of terms kmax to be taken into account for the
     * Taylor series, eq. (4.5) */
    num_abs(absz, z);
    
    num_set_d(k1, 1.0);
    num_set_d(tmp, 2.0);
    num_sub(tmp, tmp, beta);
    num_div(tmp, tmp, alpha);
    num_ceil(tmp, tmp);
    num_add(k1, tmp, k1);
    
    num_set_d(k2, 1.0);
    num_set_d(fac1, 1.0);
    num_sub(fac1, fac1, absz);
    num_mul(fac1, fac1, acc);
    num_log(fac1, fac1);
    num_log(fac2, absz);
    num_div(tmp, fac1, fac2);
    num_ceil(tmp, tmp);
    num_add(k2, tmp, k2);
    
    //k1 = (int) (ceil((2 - b)/a) + 1);
    //k2 = (int) (ceil(log(acc * (1 - abs_z))/log(abs_z)) + 1);
    kmax = num_gt(k1, k2) ? (int) num_to_d(k1) : (int) num_to_d(k2);
    delete(k1), delete(k2);

    /* Sum Taylor series */
    for (k = 0; k <= kmax; k++)
    {
        /* fac1 <- z**k */
        num_pow_d(fac1, z, (double) k); 

        /* fac2 <- rgamma(alpha * k + beta) */
        num_set_d(fac2, (double) k);
        num_mul(fac2, fac2, alpha);
        num_add(fac2, fac2, beta);
        num_rgamma(fac2, fac2);

        /* partial sum = fac1 * fac2 */
        num_mul(tmp, fac1, fac2);
        
        num_add(sum, sum, tmp);
    }
    num_set_num(res, sum);
    delete(sum), delete(z), delete(tmp), delete(fac1), delete(fac2);

    //num_print(res, true);
    //log_trace("[%s] Done.", __func__);
    
}

/* compute eq. (2.4) */
void
mittleff1 (num_t res,
           const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m] (Apply asymptotic series (2.4)) Called with parameters:\n\t    \033[1;32malpha\033[0m = %g\n\t    \033[1;32mbeta\033[0m  = %g\n\t    \033[1;32mz\033[0m = %+.14e%+.14e*I\n ",
                 __func__, num_to_d(alpha), num_to_d(beta), num_real_d(z), num_imag_d(z));
#endif     
    UNUSED(acc);

    num_t fac1, fac2, aux1, aux2;

    /*
     * Python code:
     *     fac1 = (1.0/alpha) * mp.power(z, (1.0 - beta)/alpha) * mp.exp(mp.power(z, 1.0/alpha))
     *
     * This equals to
     *     fac1 = (1.0/alpha) * aux1 * aux2
     *     aux1 = mp.power(z, (1.0 - beta)/alpha)
     *     aux2 = mp.exp(mp.power(z, 1.0/alpha))
     */
    fac1 = new(num), aux1 = new(num), aux2 = new(num);
    /* aux1 = mp.power(z, (1.0 - beta)/alpha) */
    num_set_d(aux1, 1.0);
    num_sub(aux1, aux1, beta);
    num_div(aux1, aux1, alpha);
    num_pow(aux1, z, aux1);
    /* aux2 = mp.exp(mp.power(z, 1.0/alpha)) */
    num_inv(aux2, alpha);
    num_pow(aux2, z, aux2);
    num_exp(aux2, aux2);
    /* compute fac1 */
    num_inv(fac1, alpha);
    num_mul(fac1, fac1, aux1);
    num_mul(fac1, fac1, aux2);

    /*
     * Python code:
     *     fac2 = __asymptotic(alpha, beta, z, acc)
     */
    fac2 = new(num);
    asymptotic_series(fac2, z, alpha, beta);

    /*
     * Python code:
     *     res = mp.fadd(fac1, fac2)
     */
    num_add(res, fac1, fac2);

    delete(fac1), delete(fac2);
    delete(aux1), delete(aux2);

    
    
    /* num_t fac1, fac2, _res; */
    /* fac1 = new(num), fac2 = new(num), _res = new(num); */
    /* num_set_d(fac1, 1.0); */
    /* num_sub(fac1, fac1, beta); */
    /* num_div(fac1, fac1, alpha); */
    /* num_pow(fac1, z, fac1); */
    /* num_inv(fac2, alpha); */
    /* num_pow(fac2, z, fac2); */
    /* num_exp(fac2, fac2); */
    /* num_mul(fac1, fac1, fac2); */
    /* num_inv(_res, alpha); */
    /* num_mul(fac1, fac1, _res); */
    /* asymptotic_series(fac2, z, alpha, beta); */
    /* num_add(_res, fac1, fac2); */
    /* num_set_num(res, _res); */
    /* delete(_res), delete(fac1), delete(fac2); */
}

void
mittleff2 (num_t res,
           const num_t alpha, const num_t beta,
           const num_t z,
           const num_t acc)
{
    UNUSED(acc);
    /* log_trace("[%s] alpha=%g, beta=%g, z=%g%+g, acc=%g", __func__, */
    /*           num_to_d(alpha), num_to_d(beta), num_real_d(z), num_imag_d(z), num_to_d(acc)); */
    asymptotic_series(res, z, alpha, beta);
}

/* Apply eq. (2.6) */
void
mittleff3_4 (num_t res,
             const num_t alpha,
             const num_t beta,
             const num_t z,
             const num_t acc,
             const int flag)
{
    UNUSED(acc);
    /* log_trace("[%s] alpha=%g, beta=%g, z=%g%+g, acc=%g", __func__, */
    /*           num_to_d(alpha), num_to_d(beta), num_real_d(z), num_imag_d(z), num_to_d(acc)); */
    num_t c, one, J, pi, aux, th, fac, fac1, fac2, fac3, fac4;

    c = new(num), one = new(num);
    pi = new(num), J = new(num), aux = new(num);
    fac = new(num), fac3 = new(num), fac4 = new(num);
    th = new(num), fac1 = new(num), fac2 = new(num);

    num_set_d(one, 1.0);
    num_set_d(pi, M_PI);
    num_set_d_d(J, 0.0, 1.0);

    /* compute c(theta) */
    num_inv(th, alpha);
    num_pow(th, z, th);
    num_arg(th, th);
    num_sub(th, th, pi);
    num_mul(aux, J, th);
    num_exp(c, aux);
    num_sub(c, aux, c);
    num_add(c, c, one);
    num_mul_d(c, c, 2.0);
    num_sqrt(c, c);

    if (flag == 4)
        num_neg(c, c);

    num_sub(aux, one, beta);
    num_div(aux, aux, alpha);
    num_pow(fac1, z, aux);

    num_inv(aux, alpha);
    num_pow(aux, z, aux);
    num_exp(fac2, aux);

    num_inv(aux, alpha);
    num_abs(fac4, z);
    num_pow(aux, fac4, aux);
    num_mul_d(fac4, aux, 0.5);

    num_sqrt(aux, fac4);
    num_mul(aux, c, aux);
    num_erfc(fac3, aux);

    num_mul(fac, fac1, fac2);
    num_mul(fac, fac, fac3);
    num_div(fac, fac, alpha);
    num_mul_d(fac, fac, 0.5);

    asymptotic_series(aux, z, alpha, beta);

    num_add(res, fac, aux);

    delete(c), delete(one);
    delete(pi), delete(J), delete(aux);
    delete(fac), delete(fac3), delete(fac4);
    delete(th), delete(fac1), delete(fac2);
}

void
mittleff3 (num_t res,
           const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
    /* log_trace("[%s] alpha=%g, beta=%g, z=%g%+g, acc=%g", __func__, */
    /*           num_to_d(alpha), num_to_d(beta), num_real_d(z), num_imag_d(z), num_to_d(acc)); */
    mittleff3_4(res, alpha, beta, z, acc, 3);
}

void
mittleff4 (num_t res,
           const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
    /* log_trace("[%s] alpha=%g, beta=%g, z=%g%+g, acc=%g", __func__, */
    /*           num_to_d(alpha), num_to_d(beta), num_real_d(z), num_imag_d(z), num_to_d(acc)); */
    mittleff3_4(res, alpha, beta, z, acc, 4);
}

static void
compute_rmax (num_t res,
              const num_t alpha,
              const num_t beta,
              const num_t z,
              const num_t eps)
{
    num_t r1, r2, r3, aux, den;

    r1 = new(num), r2 = new(num), r3 = new(num);
    aux = new(num), den = new(num);
    
    if (num_le_d(beta, 1.0))
    {
        if (num_ge_d(beta, 0.0))
        {
            /*
             * r1 = 2.0 * abs(z)
             * r2 = 2.0**alpha
             * r3 = (-2.0 * log(pi * eps * (2.0**beta)/12.0))**alpha
             */
            num_abs(r1, z);
            num_mul_d(r1, r1, 2.0);
            
            num_set_d(r2, 2.0);
            num_pow(r2, r2, alpha);
            
            num_set_d(r3, 2.0);
            num_pow(r3, r3, beta);
            num_mul_d(r3, r3, M_PI*num_to_d(eps)/12.0);
            num_log(r3, r3);
            num_mul_d(r3, r3, -2.0);
            num_pow(r3, r3, alpha);
        }
        else
        {
            /*
             * r1 = (2.0 * (abs(beta) + 1.0))**alpha
             * r2 = 2.0 * abs(z)
             * den = 12.0 * (abs(beta) + 2.0) * (4.0 * abs(beta))**abs(beta)
             * r3 = (-4.0 * log(pi * eps * (2.0**beta)/den))**alpha
             */
            num_abs(r1, beta);
            num_set_d(aux, 1.0);
            num_add(r1, r1, aux);
            num_mul_d(r1, r1, 2.0);
            num_pow(r1, r1, alpha);

            num_abs(r2, z);
            num_mul_d(r2, r2, 2.0);

            num_abs(aux, beta);
            num_mul_d(den, aux, 4.0);
            num_pow(den, den, aux);
            num_mul_d(den, den, 12.0);
            num_add_d(aux, aux, 2.0);
            num_mul(den, den, aux);

            num_set_d(r3, 2.0);
            num_pow(r3, r3, beta);
            num_div(r3, r3, den);
            num_mul_d(r3, r3, M_PI * num_to_d(eps));
            num_log(r3, r3);
            num_mul_d(r3, r3, -4.0);
            num_pow(r3, r3, alpha);
        }
    }
    else
    {
        if (num_ge_d(beta, 0.0))
        {
            /*
             * r1 = 1.0
             * r2 = 2.0 * abs(z)
             * r3 = (-mp.log(pi*eps/6))**alpha
             */
            num_set_d(r1, 1.0);

            num_abs(r2, z);
            num_mul_d(r2, r2, 2.0);

            num_set_d(r3, M_PI * num_to_d(eps) * 1.0/6.0);
            num_log(r3, r3);
            num_neg(r3, r3);
            num_pow(r3, r3, alpha);
        }
        else
        {
            /*
             * r1 = (abs(beta) + 1.0)**alpha
             * r2 = 2.0 * abs(z)
             * den = 6.0 * (abs(beta) + 2.0) * (2.0 * abs(beta))**abs(beta)
             * r3 = (-2.0 * log(pi * eps/den))**alpha
             */
            num_abs(r1, beta);
            num_add_d(r1, r1, 1.0);
            num_pow(r1, r1, alpha);

            num_abs(aux, z);
            num_set_d(r2, 2.0);
            num_pow(r2, r2, aux);

            num_abs(aux, beta);
            num_mul_d(den, aux, 2.0);
            num_pow(den, den, aux);
            num_mul_d(den, den, 6.0);
            num_add_d(aux, aux, 2.0);
            num_mul(den, den, aux);

            num_inv(r3, den);
            num_mul_d(r3, r3, M_PI * num_to_d(eps));
            num_log(r3, r3);
            num_mul_d(r3, r3, -2.0);
            num_pow(r3, r3,alpha);
        }
    }
    num_max3(res, r1, r2, r3);
    delete(r1), delete(r2), delete(r3);
    delete(aux), delete(den);
}

void
mittleff5_6 (num_t res,
             const num_t alpha,
             const num_t beta,
             const num_t z,
             const num_t acc,
             const num_t phi,
             const num_t c2)
{
#ifdef DEBUG
    log_info("\n[\033[1;33m%s\033[0m] Calling with parameters:\n\t    \033[1;32malpha\033[0m = %g\n\t    \033[1;32mbeta\033[0m  = %g\n\t    \033[1;32mz\033[0m = %+.14e%+.14e*I\n\t    \033[1;32mtol\033[0m = %g\n",
             __func__,
             num_to_d(alpha), num_to_d(beta), num_real_d(z), num_imag_d(z), num_to_d(acc));
#endif

    num_t rmax, from, to, aux, int1, int2, integ;

    rmax = new(num);
    from = new(num), to = new(num);
    aux = new(num);
    int1 = new(num), int2 = new(num);
    integ = new(num);

    compute_rmax(rmax, alpha, beta, z, acc);

     num_set_d(aux, 0.0);
     A(aux, z, alpha, beta, aux);        
     num_mul(aux, aux, c2);

    if (num_le_d(beta, 1.0))
    {
        num_set_d(from, 0.0);
        num_set_num(to, rmax);
        integrate_B(integ, alpha, beta, z, phi, from, to, acc);
    }
    else
    {
        num_set_d(from, 0.5);
        num_set_d(to, 1.0 + num_to_d(c2));
        num_mul(to, to, rmax);
        integrate_B(int1, alpha, beta, z, phi, from, to, acc);

        num_neg(from, phi);
        num_set_num(to, phi);
        integrate_C(int2, alpha, beta, z, phi, from, to/* , acc */);

        num_add(integ, int1, int2);
    }

    num_add(res, aux, integ);
    
    delete(rmax);
    delete(from), delete(to);
    delete(aux);
    delete(int1), delete(int2);
    delete(integ);





    
    /* log_trace("[%s] alpha=%g, beta=%g, z=%g%+g, acc=%g", __func__, */
    /*           num_to_d(alpha), num_to_d(beta), num_real_d(z), num_imag_d(z), num_to_d(acc)); */
    /* num_t zero, rmax, aux, fac1, fac2, fac3, two, one, d, phi; */

    /* one = new(num), two = new(num); */
    /* rmax = new(num), aux = new(num); */
    /* fac1 = new(num), fac2 = new(num), fac3 = new(num); */
    /* d = new(num); */
    /* zero = new(num); */
    /* num_set_d(zero, 0.0); */
    /* phi = new(num); */
    
    /* num_set_d(two, 2.0); */
    /* num_set_d(one, 1.0); */
    
    /* if (num_le_d(beta, 1.0)) */
    /* { */
    /*     /\* Compute rmax *\/ */
    /*     num_set_d(rmax, 0.0); */
    /*     if (num_ge_d(beta, 0.0)) */
    /*     { */
    /*         num_abs(fac1, z); */
    /*         num_mul(fac1, fac1, two); */
            
    /*         num_pow(fac2, two, alpha); */

    /*         num_pow(aux, two, beta); */
    /*         num_mul(aux, aux, acc); */
    /*         num_mul_d(aux, aux, M_PI/12.0); */
    /*         num_log(aux, aux); */
    /*         num_mul_d(aux, aux, -2.0); */
    /*         num_pow(fac3, aux, alpha); */

    /*         num_max3(rmax, fac1, fac2, fac3); */
    /*     } */
    /*     else */
    /*     { */
    /*         num_abs(fac1, beta); */
    /*         num_add(fac1, fac1, one); */
    /*         num_mul(fac1, fac1, two); */
    /*         num_pow(fac1, fac1, alpha); */

    /*         num_abs(fac2, z); */
    /*         num_mul(fac2, fac2, two); */

    /*         num_pow(aux, two, beta); */
    /*         num_mul(aux, aux, acc); */

    /*         num_abs(aux, beta); */
    /*         num_mul_d(d, aux, 4.0); */
    /*         num_pow(d, d, aux); */
    /*         num_add(aux, aux, two); */
    /*         num_mul(d, d, aux); */
    /*         num_mul_d(d, d, 12.0); */
            
    /*         num_div(aux, aux, d); */
    /*         num_log(aux, aux); */
    /*         num_mul_d(aux, aux, -2.0); */
    /*         num_pow(fac3, aux, alpha); */

    /*         num_max3(rmax, fac1, fac2, fac3); */
    /*     } */
    /*     /\* log_trace("[%s] rmax=%+.15e", __func__, *\/ */
    /*     /\*       num_to_d(rmax)); *\/ */

    /*     /\* int1 = numerical_integral(lambda r: B(r, alpha, beta, z, c1), 0.0, rmax, acc) *\/ */
    /*     /\* num_t int1, from, to; *\/ */
    /*     /\* int1 = new(num), from = new(num), to = new(num); *\/ */
    /*     /\* num_set_d(int1, 0.0); *\/ */
    /*     /\* num_set_d(phi, c1); *\/ */
    /*     /\* num_set_d(from, 0.0); *\/ */
    /*     /\* num_set_num(to, rmax); *\/ */
    /*     /\* integrate_B(int1, alpha, beta, z, phi, from, to); *\/ */

    /*     /\* log_trace("[%s] int1=%+.15e", __func__, *\/ */
    /*     /\*       num_to_d(int1)); *\/ */

    /*     /\* delete(int1), delete(from), delete(to); *\/ */

    /*     /\* num_set_d(res, 0.0); *\/ */

    /*     num_t from, to; */
    /*     from = new(num), to = new(num); */
    /*     num_set_d(from, 0.0); */
    /*     num_set_num(to, rmax); */
    /*     A(fac1, z, alpha, beta, zero); */
    /*     num_mul_d(fac1, fac1, c2); */
    /*     num_set_d(phi, c1); */
    /*     integrate_B(fac2, alpha, beta, z, phi, from, to, acc); */
    /*     /\* log_trace("[%s] int1=%+.15e", __func__, *\/ */
    /*     /\*       num_to_d(fac2)); *\/ */
    /*     num_add(res, fac1, fac2); */
    /*     delete(from), delete(to); */
    /* } */
    /* else */
    /* { */
    /*     if (num_ge_d(beta, 0.0)) */
    /*     { */
    /*         /\* fac1 = 2.0 * abs(z) *\/ */
    /*         num_abs(fac1, z); */
    /*         num_mul(fac1, fac1, two); */

    /*         /\* fac2 = (-log(pi * eps/6.0))**alpha *\/ */
    /*         num_mul_d(fac2, acc, M_PI/6.0); */
    /*         num_log(fac2, fac2); */
    /*         num_neg(fac2, fac2); */
    /*         num_pow(fac2, fac2, alpha); */

    /*         num_max3(rmax, one, fac1, fac2); */
    /*     } */
    /*     else */
    /*     { */
    /*         /\* fac1 = (abs(beta) + 1.0)**alpha *\/ */
    /*         num_abs(fac1, beta); */
    /*         num_add(fac1, fac1, one); */
    /*         num_pow(fac1, fac1, alpha); */

    /*         /\* fac2 = 2.0 * abs(z) *\/ */
    /*         num_abs(fac2, z); */
    /*         num_mul(fac2, fac2, two); */

    /*         /\* d = (6.0 * (abs(beta) + 2.0) * (2.0 * abs(beta))**abs(beta)) *\/ */
    /*         num_abs(aux, beta); */
    /*         num_add(aux, aux, two); */
    /*         num_abs(d, beta); */
    /*         num_mul(d, d, aux); */
    /*         num_mul_d(d, d, 12.0); */
    /*         num_abs(aux, beta); */
    /*         num_pow(d, d, aux); */

    /*         /\* fac3 = (-2.0 * log(pi * eps/d))**alpha *\/ */
    /*         num_div(fac3, acc, d); */
    /*         num_mul_d(fac3, fac3, M_PI); */
    /*         num_log(fac3, fac3); */
    /*         num_mul_d(fac3, fac3, -2.0); */
    /*         num_pow(fac3, fac3, alpha); */

    /*         num_max3(rmax, fac1, fac2, fac3); */
    /*     } */

    /*     num_t int1, int2; */
    /*     int1 = new(num), int2 = new(num); */
    /*     num_set_d(int1, 0.0), num_set_d(int2, 0.0); */

    /*     num_t from, to; */
    /*     from = new(num), to = new(num); */
        
    /*     num_set_d(from, 0.5); */
    /*     num_set_d(to, (1+c2) * num_to_d(rmax)); */
    /*     num_mul_d(phi, alpha, M_PI); */
    /*     integrate_B(int1, alpha, beta, z, phi, from, to, acc); */


    /*     num_set_d(from, -c1); */
    /*     num_set_d(to, c1); */
    /*     num_set_d(phi, 0.5); */
    /*     integrate_C(int2, alpha, beta, z, phi, from, to); */

    /*     A(aux, z, alpha, beta, zero); */
    /*     num_mul_d(aux, aux, c2); */
        
    /*     num_add(res, int1, int2); */
    /*     num_add(res, res, aux); */
        
    /*     delete(from), delete(to); */
    /*     delete(int1), delete(int2);         */
    /* } */
    /* delete(zero); */
    /* delete(aux); */
    /* delete(one), delete(two), delete(d); */
    /* delete(rmax); */
    /* delete(fac1), delete(fac2), delete(fac3); */
    /* delete(phi); */
}

/* apply eqs. (4.25) and (4.26) */
void
mittleff5 (num_t res,
           const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
    num_t phi, c2;
    phi = new(num), c2 = new(num);
    num_set_d(phi, M_PI*num_to_d(alpha));
    num_set_d(c2, 1.0);
    mittleff5_6(res, alpha, beta, z, acc, phi,c2);
    delete(phi), delete(c2);
}

void
mittleff6 (num_t res,
           const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
    num_t phi, c2;
    phi = new(num), c2 = new(num);
    num_set_d(phi, 2.0*M_PI*num_to_d(alpha)/3.0);
    num_set_d(c2, 0.0);
    mittleff5_6(res, alpha, beta, z, acc, phi, c2);
    delete(phi), delete(c2);
}

/* rhs of equation (2.3) */
void
asymptotic_series (num_t res, const num_t z, const num_t alpha, const num_t beta)
{
    int kmax;
    /*
     * Python code
     *     kmax = int(ceil((1.0/alpha)*abs(z)**(1.0/alpha)) + 1.0)
     * which can be written as
     *     kmax = int( ceil(fac1) + 1.0 )
     *     fac1 = aux1 * aux2**aux1
     *     aux1 = (1.0/alpha)
     *     aux2 = abs(z)
     */
    num_t fac1, aux1, aux2;

    fac1 = new(num), aux1 = new(num), aux2 = new(num);

    num_inv(aux1, alpha);
    num_abs(aux2, z);
    num_pow(fac1, aux2, aux1);
    num_ceil(fac1, fac1);
    kmax = ((int) num_to_d(fac1)) + 1;
    
#ifdef DEBUG
    log_info("\n[\033[1;33m%s\033[0m]\n\t    Summing %d terms of the asymptotic series\n ",
             __func__, kmax);
#endif

    /*
     * Python code:
     *     res = 0.0
     *     for k in range(1, kmax + 1):
     *         res = mp.fadd(res, mp.fmul(-mp.power(z, -k), mp.rgamma(beta - alpha*k)))
     * the loop statement can be rewritten as
     *     res = mp.fadd(res, fac1)
     *     fac1 = mp.fmul(aux1, aux2)
     *     aux1 = -mp.power(z, -k)
     *     aux2 = mp.rgamma(beta - alpha * k)
     */
    int k;
    
    num_set_d(res, 0.0);
    for (k = 1; k <= kmax; k++)
    {
        num_pow_d(aux1, z, ((double) -1.0 * k));
        num_neg(aux1, aux1);
        num_mul_d(aux2, alpha, ((double) -1.0 * k));
        num_add(aux2, aux2, beta);
        num_rgamma(aux2, aux2);
        num_mul(fac1, aux1, aux2);
        num_add(res, res, fac1);
    }

    delete(fac1), delete(aux1), delete(aux2);    
}
