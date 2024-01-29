/*
 * This file is part of mittleffc (https://github.com/mittleff/mittleffc).
 *
 * mittleffc is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * mittleffc is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * mittleffc. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file mittleff.c
 * @brief Implementation of the main functions of the library.
 */
#include "mittleff.h"

#include "types.h"

#include <math.h>

#include "partition.h"
#include "algorithm.h"

#include "flintutils.h"

#include <stdbool.h>

#ifdef DEBUG
#include "log.h"
#endif

static void compute_mittleff (acb_t res, const acb_t z, void * ctx);

/* static void */
/* show_context (void * ctx) */
/* { */
/* 		ctx_t* p = (ctx_t*) ctx; */
/* 		printf("alpha = %g\n", arbtod(p->alpha)); */
/* 		/\* printf("beta = %g\n", arbtod(p->beta)); *\/ */
/* 		/\* printf("prec = %g\n", p->prec); *\/ */
/* } */

/* Main function of the library */
int
mittleff_cmplx (double* ret,
          const double a,
          const double b,
          const double x,
          const double y,
          const unsigned short int prec)
{
		ctx_t ctx;
		acb_t z, res;

		acb_init(z);
		acb_init(res);
		arb_init(ctx.alpha);
		arb_init(ctx.beta);

		arb_set_d(ctx.alpha, a);
		arb_set_d(ctx.beta, b);
		acb_set_d_d(z, x, y);
		ctx.prec = prec;
    
		compute_mittleff(res, z, &ctx);

		acbtod(ret, res);

		acb_clear(z);
		arb_clear(res);
		arb_clear(ctx.alpha);
		arb_clear(ctx.beta);
    
		return 0;
}



static void
compute_mittleff (acb_t res, const acb_t z, void * ctx)
{
		ctx_t* p = (ctx_t*) ctx;
		
    
		arb_t zero, one, two, half;

		arb_init(zero);
		arb_init(one);
		arb_init(two);
		arb_init(half);

		arb_zero(zero);
		arb_one(one);
		arb_set_d(two, 2.0);
		arb_set_d(half, 0.5);
    
		/* Test special cases */
		if (acb_is_zero(z)) { /* z = 0 */
				acb_t b;
				acb_init(b);
				acb_set_arb(b, p->beta);
				acb_hypgeom_rgamma(res, b, p->prec);

				/* printf("z = "); */
				/* acb_printn(z, 10, 10); */
				/* printf("\n"); */
				/* printf("alpha = "); */
				/* arb_printd(p->alpha, 10); */
				/* printf("\n"); */
				/* printf("beta = "); */
				/* arb_printd(p->beta, 10); */
				/* printf("prec = %d", p->prec);				 */
				/* printf("\n"); */
				
				/* printf("res = "); */
				/* acb_printd(res, 10); */
				/* printf("\n"); */
				
				acb_clear(b);
		}
		else if (arb_eq(p->alpha, one) && arb_eq(p->beta, one)) { 
				acb_exp(res, z, p->prec);
		}
		else if (arb_eq(p->alpha, two) && arb_eq(p->beta, one)) { 
				acb_sqrt(res, z, p->prec);
				acb_cosh(res, res, p->prec);
		}
		else if (arb_eq(p->alpha, half) && arb_eq(p->beta, one)) {
				acb_t exp_z2, erfc_z;
        
				acb_init(exp_z2);
				acb_init(erfc_z);
				acb_pow_arb(exp_z2, z, two, p->prec);
				acb_exp(exp_z2, exp_z2, p->prec);
        
				acb_neg(erfc_z, z);
				acb_hypgeom_erfc(erfc_z, erfc_z, p->prec);
        
				acb_mul(res, exp_z2, erfc_z, p->prec);

				acb_clear(exp_z2);
				acb_clear(erfc_z);
		}
		else if (arb_eq(p->alpha, two) && arb_eq(p->beta, two)) { 
				acb_t n, d;

				acb_init(n);
				acb_init(d);

				acb_sqrt(n, z, p->prec);
				acb_sqrt(d, z, p->prec);
				acb_sinh(n, n, p->prec);
				acb_div(res, n, d, p->prec);

				acb_clear(n);
				acb_clear(d);
		}
		else if (in_region_G0(z, ctx)) {
				mittleff0(res, z, ctx);
		}
		else if (arb_gt(p->alpha, one)) { /* apply recursive relation (2.2) */ 
				int m, h;
				arb_t alphap, one_over_2mp1;
				acb_t zp, aux;
        
				arb_init(alphap);
				arb_init(one_over_2mp1);
				acb_init(zp);
				acb_init(aux);

				m = (int) (ceil((arbtod(p->alpha) - 1.0)/2.0) + 1.0);
				arb_set_d(one_over_2mp1, 1.0/(2.0 * m + 1.0));
				acb_zero(res);
        
				for (h = -m; h <= m; h++) {
						arb_mul(alphap, p->alpha, one_over_2mp1, p->prec);

						acb_pow_arb(zp, z, one_over_2mp1, p->prec);

						acb_set_d_d(aux, 0.0, 2.0 * M_PI * h * arbtod(one_over_2mp1));
						acb_exp(aux, aux, p->prec);
						acb_pow_arb(zp, z, one_over_2mp1, p->prec);
						acb_mul(zp, zp, aux, p->prec);

						ctx_t new_ctx = {
								.alpha = alphap,
								.beta = p->beta,
								.prec = p->prec
						};
						compute_mittleff(aux, zp, &new_ctx);

						acb_add(res, res, aux, p->prec);
				}
				acb_mul_arb(res, res, one_over_2mp1, p->prec);

				arb_clear(one_over_2mp1);
				arb_clear(alphap);
				acb_clear(zp);
				acb_clear(aux);
		} else { /* alpha <= 1 */
				if (in_region_G1(z, ctx)) 
						mittleff1(res, z, ctx);
				else if (in_region_G2(z, ctx))            
						mittleff2(res, z, ctx);
				else if (in_region_G3(z, ctx))              
						mittleff3(res, z, ctx);
				else if (in_region_G4(z, ctx))             
						mittleff4(res, z, ctx);
				else if (in_region_G5(z, ctx))             
						mittleff5(res, z, ctx);
				else if (in_region_G6(z, ctx))       
						mittleff6(res, z, ctx);
				else
						fprintf(stderr, "None of the regions");
		}
}
