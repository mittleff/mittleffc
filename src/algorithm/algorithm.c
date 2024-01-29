/*
 * This file is part of mittleffc (https://github.com/mittleff/mittleffc).
 *
 * mittleffc is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * mittleffc is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * mittleffc. If not, see <https://www.gnu.org/licenses/>.
 */

/** 
 * @file algorithm.c
 * @brief Implementation for the routines regarding the algorithm for each region of the complex plane.
 */
#include "algorithm.h"

#include "integrate.h"
/* #include "new.h" */

/* #include "utils.h" */

#include <math.h>
#include <stdbool.h>

#ifdef DEBUG
#include "log.h"
#endif

/* static void */
/* asymptotic_series (acb_t res, */
/*                    const acb_t z, */
/*                    const arb_t alpha, */
/*                    const arb_t beta); */

void
mittleff0 (acb_t res,
           const arb_t alpha,
           const arb_t beta,
           const acb_t z,
           const acb_t acc) 
{
    int k, kmax;
    arb_t k1, k2;
    acb_t absz, sum, tmp, fac1, fac2;

    acb_zero(res);

    /* acb_init(absz); */
    /* arb_init(k1); */
    /* arb_init(k2); */
    /* acb_init(fac1); */
    /* acb_init(fac2); */
    /* acb_init(tmp); */
    /* acb_init(sum); */

    /* mag_set_ui_2exp_si(tol, 1, -prec); */

    /* /\*  */
    /*  * Compute the maximum number of terms kmax to be taken into account for the */
    /*  * Taylor series, eq. (4.5)  */
    /*  *\/ */
    /* acb_abs(absz, z); */
    
    /* arb_set_d(k1, 1.0); */
    /* acb_set_d(tmp, 2.0); */
    /* acb_sub(tmp, tmp, beta); */
    /* acb_div(tmp, tmp, alpha); */
    /* acb_ceil(tmp, tmp); */
    /* acb_add(k1, tmp, k1); */
    
    /* arb_set_d(k2, 1.0); */
    /* acb_set_d(fac1, 1.0); */
    /* acb_sub(fac1, fac1, absz); */
    /* acb_mul(fac1, fac1, acc); */
    /* acb_log(fac1, fac1); */
    /* acb_log(fac2, absz); */
    /* acb_div(tmp, fac1, fac2); */
    /* acb_ceil(tmp, tmp); */
    /* acb_add(k2, tmp, k2); */
    
    /* // k1 = (int) (ceil((2 - b)/a) + 1); */
    /* // k2 = (int) (ceil(log(acc * (1 - abs_z))/log(abs_z)) + 1); */
    /* kmax = arb_gt(k1, k2) ? (int) arbtod(k1) : (int) arbtod(k2); */
    /* arb_clear(k1); */
    /* arb_clear(k2); */

    /* /\* Sum Taylor series *\/ */
    /* for (k = 0; k <= kmax; k++) */
    /* { */
    /*     /\* fac1 <- z**k *\/ */
    /*     num_pow_d(fac1, z, (double) k);  */

    /*     /\* fac2 <- rgamma(alpha * k + beta) *\/ */
    /*     num_set_d(fac2, (double) k); */
    /*     num_mul(fac2, fac2, alpha); */
    /*     num_add(fac2, fac2, beta); */
    /*     num_rgamma(fac2, fac2); */

    /*     /\* partial sum = fac1 * fac2 *\/ */
    /*     num_mul(tmp, fac1, fac2); */
        
    /*     num_add(sum, sum, tmp); */
    /* } */
    /* num_set_num(res, sum); */

    /* acb_clear(absz); */
    /* acb_clear(k1),; */
    /* acb_clear(k2); */
    /* acb_clear(fac1); */
    /* acb_clear(fac2); */
    /* acb_clear(tmp); */
    /* acb_clear(sum); */
}

/* compute eq. (2.4) */
void
mittleff1 (acb_t res,
           const arb_t alpha,
           const arb_t beta,
           const acb_t z,
           const acb_t acc) 
{

    acb_zero(res);
/* #ifdef DEBUG */
/*         log_info("\n[\033[1;33m%s\033[0m] (Apply asymptotic series (2.4)) Called with parameters:\n\t    \033[1;32malpha\033[0m = %g\n\t    \033[1;32mbeta\033[0m  = %g\n\t    \033[1;32mz\033[0m = %+.14e%+.14e*I\n ", */
/*                  __func__, num_to_d(alpha), num_to_d(beta), num_real_d(z), num_imag_d(z)); */
/* #endif      */
/*     UNUSED(acc); */

/*     num_t fac1, fac2, aux1, aux2; */

/*     /\* */
/*      * Python code: */
/*      *     fac1 = (1.0/alpha) * mp.power(z, (1.0 - beta)/alpha) * mp.exp(mp.power(z, 1.0/alpha)) */
/*      * */
/*      * This equals to */
/*      *     fac1 = (1.0/alpha) * aux1 * aux2 */
/*      *     aux1 = mp.power(z, (1.0 - beta)/alpha) */
/*      *     aux2 = mp.exp(mp.power(z, 1.0/alpha)) */
/*      *\/ */
/*     fac1 = new(num), aux1 = new(num), aux2 = new(num); */
/*     /\* aux1 = mp.power(z, (1.0 - beta)/alpha) *\/ */
/*     num_set_d(aux1, 1.0); */
/*     num_sub(aux1, aux1, beta); */
/*     num_div(aux1, aux1, alpha); */
/*     num_pow(aux1, z, aux1); */
/*     /\* aux2 = mp.exp(mp.power(z, 1.0/alpha)) *\/ */
/*     num_inv(aux2, alpha); */
/*     num_pow(aux2, z, aux2); */
/*     num_exp(aux2, aux2); */
/*     /\* compute fac1 *\/ */
/*     num_inv(fac1, alpha); */
/*     num_mul(fac1, fac1, aux1); */
/*     num_mul(fac1, fac1, aux2); */

/*     /\* */
/*      * Python code: */
/*      *     fac2 = __asymptotic(alpha, beta, z, acc) */
/*      *\/ */
/*     fac2 = new(num); */
/*     asymptotic_series(fac2, z, alpha, beta); */

/*     /\* */
/*      * Python code: */
/*      *     res = mp.fadd(fac1, fac2) */
/*      *\/ */
/*     num_add(res, fac1, fac2); */

/*     delete(fac1), delete(fac2); */
/*     delete(aux1), delete(aux2); */
}

void
mittleff2 (acb_t res,
           const arb_t alpha,
           const arb_t beta,
           const acb_t z,
           const acb_t acc) 
{
    acb_zero(res);
    /* UNUSED(acc); */
    /* asymptotic_series(res, z, alpha, beta); */
}

/* Apply eq. (2.6) */
void
mittleff3_4 (acb_t res,
             const arb_t alpha,
             const arb_t beta,
             const acb_t z,
             const acb_t acc,
             const int flag) 
{
    acb_zero(res);
    /* UNUSED(acc); */
    /* num_t c, one, J, pi, aux, th, fac, fac1, fac2, fac3, fac4; */

    /* c = new(num), one = new(num); */
    /* pi = new(num), J = new(num), aux = new(num); */
    /* fac = new(num), fac3 = new(num), fac4 = new(num); */
    /* th = new(num), fac1 = new(num), fac2 = new(num); */

    /* num_set_d(one, 1.0); */
    /* num_set_d(pi, M_PI); */
    /* num_set_d_d(J, 0.0, 1.0); */

    /* /\* compute c(theta) *\/ */
    /* num_inv(th, alpha); */
    /* num_pow(th, z, th); */
    /* num_arg(th, th); */
    /* num_sub(th, th, pi); */
    /* num_mul(aux, J, th); */
    /* num_exp(c, aux); */
    /* num_sub(c, aux, c); */
    /* num_add(c, c, one); */
    /* num_mul_d(c, c, 2.0); */
    /* num_sqrt(c, c); */

    /* if (flag == 4) */
    /*     num_neg(c, c); */

    /* num_sub(aux, one, beta); */
    /* num_div(aux, aux, alpha); */
    /* num_pow(fac1, z, aux); */

    /* num_inv(aux, alpha); */
    /* num_pow(aux, z, aux); */
    /* num_exp(fac2, aux); */

    /* num_inv(aux, alpha); */
    /* num_abs(fac4, z); */
    /* num_pow(aux, fac4, aux); */
    /* num_mul_d(fac4, aux, 0.5); */

    /* num_sqrt(aux, fac4); */
    /* num_mul(aux, c, aux); */
    /* num_erfc(fac3, aux); */

    /* num_mul(fac, fac1, fac2); */
    /* num_mul(fac, fac, fac3); */
    /* num_div(fac, fac, alpha); */
    /* num_mul_d(fac, fac, 0.5); */

    /* asymptotic_series(aux, z, alpha, beta); */

    /* num_add(res, fac, aux); */

    /* delete(c), delete(one); */
    /* delete(pi), delete(J), delete(aux); */
    /* delete(fac), delete(fac3), delete(fac4); */
    /* delete(th), delete(fac1), delete(fac2); */
}

void
mittleff3 (acb_t res,
           const arb_t alpha,
           const arb_t beta,
           const acb_t z,
           const acb_t acc) 

{
    acb_zero(res);
    //mittleff3_4(res, alpha, beta, z, acc, 3);
}

void
mittleff4 (acb_t res,
           const arb_t alpha,
           const arb_t beta,
           const acb_t z,
           const acb_t acc)
{
     acb_zero(res);
     //mittleff3_4 (res, alpha, beta, z, acc, 4);
}

/* static void */
/* compute_rmax (num_t res, */
/*               const num_t alpha, */
/*               const num_t beta, */
/*               const num_t z, */
/*               const num_t eps) */
/* { */
/*     num_t r1, r2, r3, aux, den; */

/*     r1 = new(num), r2 = new(num), r3 = new(num); */
/*     aux = new(num), den = new(num); */
    
/*     if (num_le_d(beta, 1.0)) */
/*     { */
/*         if (num_ge_d(beta, 0.0)) */
/*         { */
/*             /\* */
/*              * r1 = 2.0 * abs(z) */
/*              * r2 = 2.0**alpha */
/*              * r3 = (-2.0 * log(pi * eps * (2.0**beta)/12.0))**alpha */
/*              *\/ */
/*             num_abs(r1, z); */
/*             num_mul_d(r1, r1, 2.0); */
            
/*             num_set_d(r2, 2.0); */
/*             num_pow(r2, r2, alpha); */
            
/*             num_set_d(r3, 2.0); */
/*             num_pow(r3, r3, beta); */
/*             num_mul_d(r3, r3, M_PI*num_to_d(eps)/12.0); */
/*             num_log(r3, r3); */
/*             num_mul_d(r3, r3, -2.0); */
/*             num_pow(r3, r3, alpha); */
/*         } */
/*         else */
/*         { */
/*             /\* */
/*              * r1 = (2.0 * (abs(beta) + 1.0))**alpha */
/*              * r2 = 2.0 * abs(z) */
/*              * den = 12.0 * (abs(beta) + 2.0) * (4.0 * abs(beta))**abs(beta) */
/*              * r3 = (-4.0 * log(pi * eps * (2.0**beta)/den))**alpha */
/*              *\/ */
/*             num_abs(r1, beta); */
/*             num_set_d(aux, 1.0); */
/*             num_add(r1, r1, aux); */
/*             num_mul_d(r1, r1, 2.0); */
/*             num_pow(r1, r1, alpha); */

/*             num_abs(r2, z); */
/*             num_mul_d(r2, r2, 2.0); */

/*             num_abs(aux, beta); */
/*             num_mul_d(den, aux, 4.0); */
/*             num_pow(den, den, aux); */
/*             num_mul_d(den, den, 12.0); */
/*             num_add_d(aux, aux, 2.0); */
/*             num_mul(den, den, aux); */

/*             num_set_d(r3, 2.0); */
/*             num_pow(r3, r3, beta); */
/*             num_div(r3, r3, den); */
/*             num_mul_d(r3, r3, M_PI * num_to_d(eps)); */
/*             num_log(r3, r3); */
/*             num_mul_d(r3, r3, -4.0); */
/*             num_pow(r3, r3, alpha); */
/*         } */
/*     } */
/*     else */
/*     { */
/*         if (num_ge_d(beta, 0.0)) */
/*         { */
/*             /\* */
/*              * r1 = 1.0 */
/*              * r2 = 2.0 * abs(z) */
/*              * r3 = (-mp.log(pi*eps/6))**alpha */
/*              *\/ */
/*             num_set_d(r1, 1.0); */

/*             num_abs(r2, z); */
/*             num_mul_d(r2, r2, 2.0); */

/*             num_set_d(r3, M_PI * num_to_d(eps) * 1.0/6.0); */
/*             num_log(r3, r3); */
/*             num_neg(r3, r3); */
/*             num_pow(r3, r3, alpha); */
/*         } */
/*         else */
/*         { */
/*             /\* */
/*              * r1 = (abs(beta) + 1.0)**alpha */
/*              * r2 = 2.0 * abs(z) */
/*              * den = 6.0 * (abs(beta) + 2.0) * (2.0 * abs(beta))**abs(beta) */
/*              * r3 = (-2.0 * log(pi * eps/den))**alpha */
/*              *\/ */
/*             num_abs(r1, beta); */
/*             num_add_d(r1, r1, 1.0); */
/*             num_pow(r1, r1, alpha); */

/*             num_abs(aux, z); */
/*             num_set_d(r2, 2.0); */
/*             num_pow(r2, r2, aux); */

/*             num_abs(aux, beta); */
/*             num_mul_d(den, aux, 2.0); */
/*             num_pow(den, den, aux); */
/*             num_mul_d(den, den, 6.0); */
/*             num_add_d(aux, aux, 2.0); */
/*             num_mul(den, den, aux); */

/*             num_inv(r3, den); */
/*             num_mul_d(r3, r3, M_PI * num_to_d(eps)); */
/*             num_log(r3, r3); */
/*             num_mul_d(r3, r3, -2.0); */
/*             num_pow(r3, r3,alpha); */
/*         } */
/*     } */
/*     num_max3(res, r1, r2, r3); */
/*     delete(r1), delete(r2), delete(r3); */
/*     delete(aux), delete(den); */
/* } */

void
mittleff5_6 (acb_t res,
             const arb_t alpha,
             const arb_t beta,
             const acb_t z,
             const acb_t acc,
             const arb_t phi,
             const arb_t c2)    
{
    acb_zero(res);
/* #ifdef DEBUG */
/*     log_info("\n[\033[1;33m%s\033[0m] Calling with parameters:\n\t    \033[1;32malpha\033[0m = %g\n\t    \033[1;32mbeta\033[0m  = %g\n\t    \033[1;32mz\033[0m = %+.14e%+.14e*I\n\t    \033[1;32mtol\033[0m = %g\n", */
/*              __func__, */
/*              num_to_d(alpha), num_to_d(beta), num_real_d(z), num_imag_d(z), num_to_d(acc)); */
/* #endif */

/*     num_t rmax, from, to, aux, int1, int2, integ; */

/*     rmax = new(num); */
/*     from = new(num), to = new(num); */
/*     aux = new(num); */
/*     int1 = new(num), int2 = new(num); */
/*     integ = new(num); */

/*     compute_rmax(rmax, alpha, beta, z, acc); */

/*      num_set_d(aux, 0.0); */
/*      A(aux, z, alpha, beta, aux);         */
/*      num_mul(aux, aux, c2); */

/*     if (num_le_d(beta, 1.0)) */
/*     { */
/*         num_set_d(from, 0.0); */
/*         num_set_num(to, rmax); */
/*         integrate_B(integ, alpha, beta, z, phi, from, to, acc); */
/*     } */
/*     else */
/*     { */
/*         num_set_d(from, 0.5); */
/*         num_set_d(to, 1.0 + num_to_d(c2)); */
/*         num_mul(to, to, rmax); */
/*         integrate_B(int1, alpha, beta, z, phi, from, to, acc); */

/*         num_neg(from, phi); */
/*         num_set_num(to, phi); */
/*         integrate_C(int2, alpha, beta, z, phi, from, to/\* , acc *\/); */

/*         num_add(integ, int1, int2); */
/*     } */

/*     num_add(res, aux, integ); */
    
/*     delete(rmax); */
/*     delete(from), delete(to); */
/*     delete(aux); */
/*     delete(int1), delete(int2); */
/*     delete(integ); */
}

/* apply eqs. (4.25) and (4.26) */
void
mittleff5 (acb_t res,
           const arb_t alpha,
           const arb_t beta,
           const acb_t z,
           const acb_t acc)
{
    acb_zero(res);
    /* num_t phi, c2; */
    /* phi = new(num), c2 = new(num); */
    /* num_set_d(phi, M_PI*num_to_d(alpha)); */
    /* num_set_d(c2, 1.0); */
    /* mittleff5_6(res, alpha, beta, z, acc, phi,c2); */
    /* delete(phi), delete(c2); */
}

void
mittleff6 (acb_t res,
           const arb_t alpha,
           const arb_t beta,
           const acb_t z,
           const acb_t acc)
{
    acb_zero(res);
    /* num_t phi, c2; */
    /* phi = new(num), c2 = new(num); */
    /* num_set_d(phi, 2.0*M_PI*num_to_d(alpha)/3.0); */
    /* num_set_d(c2, 0.0); */
    /* mittleff5_6(res, alpha, beta, z, acc, phi, c2); */
    /* delete(phi), delete(c2); */
}

/* rhs of equation (2.3) */
/* void */
/* asymptotic_series (acb_t res, */
/*                    const acb_t z, */
/*                    const arb_t alpha, */
/*                    const arb_t beta)     */
/* { */
/*     int kmax; */
/*     /\* */
/*      * Python code */
/*      *     kmax = int(ceil((1.0/alpha)*abs(z)**(1.0/alpha)) + 1.0) */
/*      * which can be written as */
/*      *     kmax = int( ceil(fac1) + 1.0 ) */
/*      *     fac1 = aux1 * aux2**aux1 */
/*      *     aux1 = (1.0/alpha) */
/*      *     aux2 = abs(z) */
/*      *\/ */
/*     num_t fac1, aux1, aux2; */

/*     fac1 = new(num), aux1 = new(num), aux2 = new(num); */

/*     num_inv(aux1, alpha); */
/*     num_abs(aux2, z); */
/*     num_pow(fac1, aux2, aux1); */
/*     num_ceil(fac1, fac1); */
/*     kmax = ((int) num_to_d(fac1)) + 1; */
    
/* #ifdef DEBUG */
/*     log_info("\n[\033[1;33m%s\033[0m]\n\t    Summing %d terms of the asymptotic series\n ", */
/*              __func__, kmax); */
/* #endif */

/*     /\* */
/*      * Python code: */
/*      *     res = 0.0 */
/*      *     for k in range(1, kmax + 1): */
/*      *         res = mp.fadd(res, mp.fmul(-mp.power(z, -k), mp.rgamma(beta - alpha*k))) */
/*      * the loop statement can be rewritten as */
/*      *     res = mp.fadd(res, fac1) */
/*      *     fac1 = mp.fmul(aux1, aux2) */
/*      *     aux1 = -mp.power(z, -k) */
/*      *     aux2 = mp.rgamma(beta - alpha * k) */
/*      *\/ */
/*     int k; */
    
/*     num_set_d(res, 0.0); */
/*     for (k = 1; k <= kmax; k++) */
/*     { */
/*         num_pow_d(aux1, z, ((double) -1.0 * k)); */
/*         num_neg(aux1, aux1); */
/*         num_mul_d(aux2, alpha, ((double) -1.0 * k)); */
/*         num_add(aux2, aux2, beta); */
/*         num_rgamma(aux2, aux2); */
/*         num_mul(fac1, aux1, aux2); */
/*         num_add(res, res, fac1); */
/*     } */

/*     delete(fac1), delete(aux1), delete(aux2);     */
/* } */
