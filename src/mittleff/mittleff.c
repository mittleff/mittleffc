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
//#include <stdio.h>
//#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <assert.h>

#include "partition.h"
#include "algorithm.h"

#include "flintutils.h"

#include <gsl/gsl_math.h>
#include <stdbool.h>

#ifdef DEBUG
#include "log.h"
#endif

// Converts an arb_t number to double.
/* static double */
/* arbtod (const arb_t x) */
/* { */
/*     return arf_get_d(arb_midref(x), ARF_RND_NEAR); */
/* } */

static void
_mittleff (acb_t res,
           const arb_t alpha,
           const arb_t beta,
           const acb_t z,
           const arb_t tol)
{
    arb_t zero, one, two, half;
    slong prec = 64;

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
        acb_hypgeom_rgamma(res, z, prec);
    }
    else if (arb_eq(alpha, one) && arb_eq(beta, one)) { /* exp(z) */
        acb_exp(res, z, prec);
    }
    else if (arb_eq(alpha, two) && arb_eq(beta, one)) { /* cosh(sqrt(z)) */
        acb_sqrt(res, z, prec);
        acb_cosh(res, res, prec);
    }
    else if (arb_eq(alpha, half) && arb_eq(beta, one)) { /* exp(z^2)*erfc(-z) */
        acb_t exp_z2, erfc_z;
        
        acb_init(exp_z2);
        acb_init(erfc_z);
        acb_pow(exp_z2, z, two, prec);
        acb_exp(exp_z2, exp_z2);
        
        acb_neg(erfc_z, z);
        acb_hypgeom_erfc(erfc_z, erfc_z, prec);
        
        acb_mul(res, exp_z2, erfc_z, prec);

        acb_clear(exp_z2);
        acb_clear(erfc_z);
    }
    else if (arb_eq(alpha, two) && arb_eq(beta, two)) { /* sinh(sqrt(z))/sqrt(z) */
        acb_t n, d;

        acb_init(n);
        acb_init(d);

        acb_sqrt(n, z, prec);
        acb_sqrt(d, z, prec);
        acb_sinh(n, n, prec);
        acb_div(res, n, d, prec);

        acb_clear(n);
        acb_clear(d);
    }
    else if (in_region_G0(z)) {
        mittleff0(res, alpha, beta, z, tol);
    }
    else if (arb_gt(alpha, one)) { /* apply recursive relation (2.2) */ 
        int m, h;
        arb_t alphap, one_over_2mp1;
        acb_t zp, aux;
        
        arb_init(alphap);
        arb_init(one_over_2mp1);
        acb_init(zp);
        acb_init(aux);

        m = (int) (ceil((arbtod(alpha) - 1.0)/2.0) + 1.0);
        arb_set_d(one_over_2mp1, 1.0/(2.0 * m + 1.0));
        acb_zero(res);
        
        for (h = -m; h <= m; h++) {
            arb_mul(alphap, alpha, one_over_2mp1, prec);

            acb_pow_arb(zp, z, one_over_2mp1, prec);

            acb_set_d_d(aux, 0.0, 2.0 * M_PI * h * arbtod(one_over_2mp1));
            acb_exp(aux, aux, prec);
            acb_pow_arb(zp, z, one_over_2mp1, prec);
            acb_mul(zp, zp, aux, prec);

            _mittleff(aux, alphap, beta, zp, tol);

            acb_add(res, res, aux, prec);
        }
        acb_mul_arb(res, res, one_over_2mp1, prec);

        arb_clear(one_over_2mp1);
        arb_clear(alphap);
        acb_clear(zp);
        acb_clear(aux);
    } else { /* alpha <= 1 */
        if (in_region_G1(z, alpha, tol)) {            
            mittleff1(res, alpha, beta, z, tol);
        }
        else if (in_region_G2(z, alpha, tol)) {           
            mittleff2(res, alpha, beta, z, tol);
        }
        else if (in_region_G3(z, alpha, tol)) {             
            mittleff3(res, alpha, beta, z, tol);
        }
        else if (in_region_G4(z, alpha, tol)) {            
            mittleff4(res, alpha, beta, z, tol);
        }
        else if (in_region_G5(z, alpha, tol)) {            
            mittleff5(res, alpha, beta, z, tol);
        }
        else if (in_region_G6(z, alpha, tol)) {            
            mittleff6(res, alpha, beta, z, tol);
        }
        else
            fprintf(stderr, "None of the regions");
    }
}

/* Main function of the library */
int
mittleff_cmplx (double* res,
          const double alpha,
          const double beta,
          const double x, const double y,
          const double tol)
{
    assert (alpha > 0.0);

    arb_t _alpha, _beta, _tol;
    acb_t _z;

    arb_init(_alpha);
    arb_init(_beta);
    arb_init(_tol);
    acb_init(_res);

    arb_set_d(_alpha, alpha);
    arb_set_d(_beta, beta);
    arb_set_d(_tol, tol);
    acb_set_d_d(_z, x, y);

    _mittleff(_res, _alpha, _beta, _z, _tol);

    //num_to_d_d(res, _res);
    
    arb_clear(_alpha);
    arb_clear(_beta);
    arb_clear(_tol);
    acb_clear(_z)
    arb_clear(_res);

    return EXIT_SUCCESS;
}
