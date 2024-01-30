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

#include "flintutils.h"
#include "types.h"

#ifdef DEBUG
#include "log.h"
#endif

/* static void */
/* asymptotic_series (acb_t res, */
/*                    const acb_t z, */
/*                    const arb_t alpha, */
/*                    const arb_t beta); */

void
mittleff0 (acb_t res, const acb_t z, void * ctx)
{
    ctx_t* p = (ctx_t*) ctx;
    //double a, b;
    int k, kmax;
    arb_t one, aux, absz, k1, k2, eps;
    acb_t  sum, tmp, fac1, fac2;

#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m] Called with parameters:\n\t    \033[1;32malpha\033[0m = %g\n\t    \033[1;32mbeta\033[0m  = %g\n\t    \033[1;32mz\033[0m = %+.14e%+.14e*I\n ",
                 __func__, arbtod(p->alpha), arbtod(p->beta), acb_real_d(z), acb_imag_d(z));
#endif    

    arb_init(one);
    arb_init(absz);
    arb_init(k1);
    arb_init(k2);
    acb_init(fac1);
    acb_init(fac2);
    acb_init(tmp);
    acb_init(sum);
    arb_init(aux);

    arb_one(one);

    arb_mul_2exp_si(eps, one, -p->prec);

    /*
     * Compute the maximum number of terms kmax to be taken into account for the
     * Taylor series, eq. (4.5)
     */
    acb_abs(absz, z, p->prec);
    
    arb_set_d(aux, 2.0);
    arb_sub(aux, aux, p->beta, p->prec);
    arb_div(aux, aux, p->alpha, p->prec);
    arb_ceil(aux, aux, p->prec);
    arb_add(k1, aux, one, p->prec);
    
    arb_sub(fac1, one, absz, p->prec);
    arb_mul(fac1, fac1, eps, p->prec);
    arb_log(fac1, fac1, p->prec);
    arb_log(fac2, absz, p->prec);
    arb_div(aux, fac1, fac2, p->prec);
    arb_ceil(aux, aux, p->prec);
    arb_add(k2, aux, one, p->prec);
    
    // k1 = (int) (ceil((2 - b)/a) + 1);
    // k2 = (int) (ceil(log(acc * (1 - abs_z))/log(abs_z)) + 1);
    kmax = arb_gt(k1, k2) ? (int) arbtod(k1) : (int) arbtod(k2);

    acb_t kk;
    acb_zero(res);
    acb_init(kk);
    /* Sum Taylor series */
    for (k = 0; k <= kmax; k++)
    {
        acb_set_d(kk, (double)k);

        /* fac1 <- z**k */
        acb_pow(fac1, z, kk, p->prec);

        /* fac2 <- rgamma(alpha * k + beta) */
        acb_mul_arb(fac2, kk, p->alpha, p->prec);
        acb_add_arb(fac2, fac2, p->beta, p->prec);
        acb_rgamma(fac2, fac2, p->prec);

        /* partial sum = fac1 * fac2 */
        acb_mul(tmp, fac1, fac2, p->prec);
        
        acb_add(res, res, tmp, p->prec);
    }
    acb_clear(kk);
    arb_clear(aux);
    acb_clear(fac1);
    acb_clear(fac2);
    acb_clear(tmp);
    acb_clear(sum);
    arb_clear(one);
    arb_clear(k1);
    arb_clear(k2);
}

/* compute eq. (2.4) */
void
mittleff1 (acb_t res, const acb_t z, void * ctx)
{
    ctx_t* p = (ctx_t*) ctx;
    
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m] Called with parameters:\n\t    \033[1;32malpha\033[0m = %g\n\t    \033[1;32mbeta\033[0m  = %g\n\t    \033[1;32mz\033[0m = %+.14e%+.14e*I\n ",
                 __func__, arbtod(p->alpha), arbtod(p->beta), acb_real_d(z), acb_imag_d(z));
#endif     

    acb_zero(res);
/* #ifdef DEBUG */
/*         log_info("\n[\033[1;33m%s\033[0m] (Apply asymptotic series (2.4)) Called with parameters:\n\t    \033[1;32malpha\033[0m = %g\n\t    \033[1;32mbeta\033[0m  = %g\n\t    \033[1;32mz\033[0m = %+.14e%+.14e*I\n ", */
/*                  __func__, num_to_d(alpha), num_to_d(beta), num_real_d(z), num_imag_d(z)); */
/* #endif      */
/*     UNUSED(acc); */

/*     num_t fac1, fac2, aux1, aux2; */

	acb_t fac1, fac2, aux1, aux2;

	acb_init(fac1);
    acb_init(fac2);
    acb_init(aux1);
    acb_init(aux2);

    /*
     * Python code:
     *     fac1 = (1.0/alpha) * mp.power(z, (1.0 - beta)/alpha) * mp.exp(mp.power(z, 1.0/alpha))
     *
     * This equals to
     *     fac1 = (1.0/alpha) * aux1 * aux2
     *     aux1 = mp.power(z, (1.0 - beta)/alpha)
     *     aux2 = mp.exp(mp.power(z, 1.0/alpha))
     */
    /* aux1 = mp.power(z, (1.0 - beta)/alpha) */
    acb_set_d(aux1, 1.0);
    acb_sub_arb(aux1, aux1, p->beta, p->prec);
    acb_div_arb(aux1, aux1, p->alpha, p->prec);
    acb_pow(aux1, z, aux1, p->prec);
    /* aux2 = mp.exp(mp.power(z, 1.0/alpha)) */
	acb_set_arb(aux2, p->alpha);
    acb_inv(aux2, aux2, p->prec);
    acb_pow(aux2, z, aux2, p->prec);
    acb_exp(aux2, aux2, p->prec);
    /* compute fac1 */
	acb_set_arb(fac1, p->alpha);
    acb_inv(fac1, fac1, p->prec);
    acb_mul(fac1, fac1, aux1, p->prec);
    acb_mul(fac1, fac1, aux2, p->prec);

    /*
     * Python code:
     *     fac2 = __asymptotic(alpha, beta, z, acc)
     */
    acb_init(fac2);
    asymptotic_series(fac2, z, ctx);

    /*
     * Python code:
     *     res = mp.fadd(fac1, fac2)
     */
    acb_add(res, fac1, fac2, p->prec);

    acb_clear(fac1);
    acb_clear(fac2);
    acb_clear(aux1);
    acb_clear(aux2);
}

void
mittleff2 (acb_t res, const acb_t z, void * ctx) 
{
		ctx_t* p = (ctx_t*) ctx;
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m] Called with parameters:\n\t    \033[1;32malpha\033[0m = %g\n\t    \033[1;32mbeta\033[0m  = %g\n\t    \033[1;32mz\033[0m = %+.14e%+.14e*I\n ",
                 __func__, arbtod(p->alpha), arbtod(p->beta), acb_real_d(z), acb_imag_d(z));
#endif 		
		//acb_zero(res);
		/* UNUSED(acc); */
		asymptotic_series(res, z, ctx);
}

/* Apply eq. (2.6) */
void
mittleff3_4 (acb_t res,
             const acb_t z,
			 void * ctx,
             const int flag) 
{
		ctx_t* p = (ctx_t*) ctx;
		acb_t c, one, J, pi, aux, th, fac, fac1, fac2, fac3, fac4;

		acb_init(c);
		acb_init(one);
		acb_init(pi);
		acb_init(J);
		acb_init(aux);
		acb_init(fac);
		acb_init(fac3);
		acb_init(fac4);
		acb_init(th);
		acb_init(fac1);
		acb_init(fac2);

		acb_one(one);
		acb_const_pi(pi, p->prec);
		acb_onei(J);

		/* compute c(theta) */
		acb_set_arb(th, p->alpha);
		acb_inv(th, th, p->prec);
		acb_pow(th, z, th, p->prec);
		acb_arg(th, th, p->prec);
		acb_sub(th, th, pi, p->prec);
		acb_mul(aux, J, th, p->prec);
		acb_exp(c, aux, p->prec);
		acb_sub(c, aux, c, p->prec);
		acb_add(c, c, one, p->prec);
		acb_set_d(aux, 2.0);
		acb_mul(c, c, aux, p->prec);
		acb_sqrt(c, c, p->prec);

		if (flag == 4)
				acb_neg(c, c);

		acb_sub_arb(aux, one, p->beta, p->prec);
		acb_div_arb(aux, aux, p->alpha, p->prec);
		acb_pow(fac1, z, aux, p->prec);

		acb_set_arb(aux, p->alpha);
		acb_inv(aux, aux, p->prec);
		acb_pow(aux, z, aux, p->prec);
		acb_exp(fac2, aux, p->prec);

		acb_set_arb(aux, p->alpha);
		acb_inv(aux, aux, p->prec);
		acb_abs(fac4, z, p->prec);
		acb_pow(aux, fac4, aux, p->prec);
		acb_set_d(fac4, 0.5);
		acb_mul(fac4, aux, fac4, p->prec);

		acb_sqrt(aux, fac4, p->prec);
		acb_mul(aux, c, aux, p->prec);
		acb_hypgeom_erfc(fac3, aux, p->prec);

		acb_mul(fac, fac1, fac2, p->prec);
		acb_mul(fac, fac, fac3, p->prec);
		acb_div_arb(fac, fac, p->alpha, p->prec);
		acb_set_d(aux, 0.5);
		acb_mul(fac, fac, aux, p->prec);

		asymptotic_series(aux, z, ctx);

		acb_add(res, fac, aux, p->prec);

		acb_clear(c);
		acb_clear(one);
		acb_clear(pi);
		acb_clear(J);
		acb_clear(aux);
		acb_clear(fac);
		acb_clear(fac3);
		acb_clear(fac4);
		acb_clear(th);
		acb_clear(fac1);
		acb_clear(fac2);
}

void
mittleff3 (acb_t res, const acb_t z, void * ctx) 

{
		ctx_t* p = (ctx_t*) ctx;
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m] Called with parameters:\n\t    \033[1;32malpha\033[0m = %g\n\t    \033[1;32mbeta\033[0m  = %g\n\t    \033[1;32mz\033[0m = %+.14e%+.14e*I\n ",
                 __func__, arbtod(p->alpha), arbtod(p->beta), acb_real_d(z), acb_imag_d(z));
#endif 		
		mittleff3_4(res, z, ctx, 3);
}

void
mittleff4 (acb_t res, const acb_t z, void * ctx)
{
		ctx_t* p = (ctx_t*) ctx;
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m] Called with parameters:\n\t    \033[1;32malpha\033[0m = %g\n\t    \033[1;32mbeta\033[0m  = %g\n\t    \033[1;32mz\033[0m = %+.14e%+.14e*I\n ",
                 __func__, arbtod(p->alpha), arbtod(p->beta), acb_real_d(z), acb_imag_d(z));
#endif 		
	 mittleff3_4(res, z, ctx, 4);
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
mittleff5 (acb_t res, const acb_t z, void * ctx)
{
		ctx_t* p = (ctx_t*) ctx;
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m] Called with parameters:\n\t    \033[1;32malpha\033[0m = %g\n\t    \033[1;32mbeta\033[0m  = %g\n\t    \033[1;32mz\033[0m = %+.14e%+.14e*I\n ",
                 __func__, arbtod(p->alpha), arbtod(p->beta), acb_real_d(z), acb_imag_d(z));
#endif 	
    acb_zero(res);
    /* num_t phi, c2; */
    /* phi = new(num), c2 = new(num); */
    /* num_set_d(phi, M_PI*num_to_d(alpha)); */
    /* num_set_d(c2, 1.0); */
    /* mittleff5_6(res, alpha, beta, z, acc, phi,c2); */
    /* delete(phi), delete(c2); */
}

void
mittleff6 (acb_t res, const acb_t z, void * ctx)
{
		ctx_t* p = (ctx_t*) ctx;
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m] Called with parameters:\n\t    \033[1;32malpha\033[0m = %g\n\t    \033[1;32mbeta\033[0m  = %g\n\t    \033[1;32mz\033[0m = %+.14e%+.14e*I\n ",
                 __func__, arbtod(p->alpha), arbtod(p->beta), acb_real_d(z), acb_imag_d(z));
#endif 	
    acb_zero(res);
    /* num_t phi, c2; */
    /* phi = new(num), c2 = new(num); */
    /* num_set_d(phi, 2.0*M_PI*num_to_d(alpha)/3.0); */
    /* num_set_d(c2, 0.0); */
    /* mittleff5_6(res, alpha, beta, z, acc, phi, c2); */
    /* delete(phi), delete(c2); */
}


static int
compute_kmax_asymptotic (const acb_t z, void * ctx)
{
		int kmax;
		ctx_t* p = (ctx_t*) ctx;
		arb_t res, aux1, aux2;

		arb_init(res);
		arb_init(aux1);
		arb_init(aux2);
		
		/*
		 * Python code
		 *     kmax = int(ceil((1.0/alpha)*abs(z)**(1.0/alpha)) + 1.0)
		 * which can be written as
		 *     kmax = int( ceil(fac1) + 1.0 )
		 *     fac1 = aux1 * aux2**aux1
		 *     aux1 = (1.0/alpha)
		 *     aux2 = abs(z)
		 */
		arb_inv(aux1, p->alpha, p->prec);
		acb_abs(aux2, z, p->prec);
		arb_pow(res, aux2, aux1, p->prec);
		arb_ceil(res, res, p->prec);

		kmax = ((int) arbtod(res)) + 1;

		arb_clear(res);
		arb_clear(aux1);
		arb_clear(aux2);

		return kmax;
}

/* rhs of equation (2.3) */
void
asymptotic_series (acb_t res, const acb_t z, void * ctx)
{
		ctx_t* p = (ctx_t*) ctx;
		int k, kmax;
		acb_t fac1, aux1, aux2, kk;

		acb_init(fac1);
		acb_init(aux1);
		acb_init(aux2);
		acb_init(kk);

		kmax = compute_kmax_asymptotic(z, ctx);
    
#ifdef DEBUG
		log_info("\n[\033[1;33m%s\033[0m]\n\t    Summing %d terms\n",
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
		acb_zero(res);
		for (k = 1; k <= kmax; k++) {
				/* printf("k = %d\n", k); */
				acb_set_d(kk, ((double) -1.0 * k));
				acb_pow(aux1, z, kk, p->prec);
				acb_neg(aux1, aux1);
				acb_mul_arb(aux2, kk, p->alpha, p->prec);
				acb_add_arb(aux2, aux2, p->beta, p->prec);
				acb_rgamma(aux2, aux2, p->prec);
				acb_mul(fac1, aux1, aux2, p->prec);
				acb_add(res, res, fac1, p->prec);
		}
		acb_clear(kk);
		acb_clear(fac1);
		acb_clear(aux1);
		acb_clear(aux2);
}
