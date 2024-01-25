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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <assert.h>

#include "log.h"
#include "num.h"
// #include "new.h"
#include "partition.h"
#include "algorithm.h"

#include <gsl/gsl_math.h>
#include <stdbool.h>

// #define MAX2(a,b) ((a>b)?(a):(b))
// #define MAX3(a,b,c) MAX2(a,MAX2(b,c))

static void
_mittleff (num_t res,
           const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
    /* log_trace("[%s] alpha=%g, beta=%g, z=%g%+g, acc=%g", __func__, */
    /*           num_to_d(alpha), num_to_d(beta), num_real_d(z), num_imag_d(z), num_to_d(acc)); */
    /* Test special cases */
    if (num_is_zero(z)) /* z = 0 */
    {
        num_rgamma(res, beta);
        log_trace("[%s] z = 0.0, res = 1/Gamma(%g) = %g",
                  __func__,
                  num_to_d(beta),
                  num_to_d(res));
    }
    else if (num_eq_d(alpha, 1.0) && num_eq_d(beta, 1.0)) /* exp(z) */
    {
        num_exp(res, z);
        log_trace("[%s] alpha=beta=1, res = exp(%g+%g*I) = %g+%g*I",
                  __func__,
                  num_to_complex(z),
                  num_to_complex(res));
    }
    else if (num_eq_d(alpha, 2.0) && num_eq_d(beta, 1.0)) /* cosh(sqrt(z)) */
    {
        num_sqrt(res, z);
        num_cosh(res, res);
        log_trace("[%s] alpha=2, beta=1, res = cosh(sqrt(%g+%g*I)) = %g+%g*I",
                  __func__,
                  num_to_complex(z),
                  num_to_complex(res));
    }
    else if (num_eq_d(alpha, 0.5) && num_eq_d(beta, 1.0)) /* exp(z^2)*erfc(-z) */
    {
        num_t exp_z2, erfc_z;
        exp_z2 = num_init(), erfc_z = num_init();
        num_pow_d(exp_z2, z, 2.0);
        num_exp(exp_z2, exp_z2);
        log_trace("[%s] exp(z^2) = %g+%g*I",
                  __func__,
                  num_to_complex(exp_z2));
        
        num_neg(erfc_z, z);
        num_erfc(erfc_z, erfc_z);
        log_trace("[%s] erfc(-z) = %g+%g*I",
                  __func__,
                  num_to_complex(erfc_z));
        
        num_mul(res, exp_z2, erfc_z);
        num_clear(exp_z2), num_clear(erfc_z);

        log_trace("[%s] alpha=0.5, beta=1, res = exp(z^2)*erfc(-z) = %g+%g*I",
                  __func__,
                  num_to_complex(res));
    }
    else if (num_eq_d(alpha, 2.0) && num_eq_d(beta, 2.0)) /* sinh(sqrt(z))/sqrt(z) */
    {
        num_t n, d;
        n = num_init();
        d = num_init();
        num_sqrt(n, z);
        num_sqrt(d, z);
        num_sinh(n, n);
        num_div(res, n, d);
        num_clear(n), num_clear(d);
    }
    else if (in_region_G0(z))
        mittleff0(res, alpha, beta, z, acc);
    else if (num_gt_d(alpha, 1.0)) /* apply recursive relation (2.2) */
    {
        /* log_trace("[%s] alpha = %+.5e > 1, applying recursive relation", __func__, alpha); */
        /* const int m = (int) (ceil((alpha - 1)/2.0) + 1); */
        /* const double one_over_2mp1 = 1.0/(2.0 * m + 1.0); */

        /* num_t sum, th, tmp, exp_th, newz; */
        /* sum = num_init(), th = num_init(), z = num_init(), tmp = num_init(), newz = num_init(), exp_th = num_init(); */
        /* num_zero(sum); */
        /* for (int h = -m; h <= m; h++) */
        /* { */
        /*     log_trace("[%s] h = %d", __func__, h); */
        /*     num_set_d_d(th, 0.0, 2.0 * M_PI * h * one_over_2mp1); */
        /*     num_exp(exp_th, th); */
        /*     num_pow_d(newz, z, one_over_2mp1); */
        /*     num_mul(newz, newz, exp_th); */
        /*     const double complex _newz = num_to_complex(newz); */
        /*     _mittleff(_res, alpha * one_over_2mp1, beta, creal(_newz), cimag(_newz), acc); */
        /*     num_set_d_d(tmp, _res[0], _res[1]); */
        /*     num_add(sum, sum, tmp); */
        /* } */
        /* num_mul_d(sum, sum, 1.0/(2.0*m + 1.0)); */
        /* num_set(_res, sum); */
        /* num_clear(sum), num_clear(th), num_clear(tmp), num_clear(exp_th); */
    }
    else /* alpha <= 1 */
    {
        if (in_region_G1(z, alpha, acc))
            mittleff1(res, alpha, beta, z, acc);
        else if (in_region_G2(z, alpha, acc))
            mittleff2(res, alpha, beta, z, acc);
        else if (in_region_G3(z, alpha, acc))
            mittleff3(res, alpha, beta, z, acc);
        else if (in_region_G4(z, alpha, acc))
            mittleff4(res, alpha, beta, z, acc);
        else if (in_region_G5(z, alpha, acc))
            mittleff5(res, alpha, beta, z, acc);
        else if (in_region_G6(z, alpha, acc))
            mittleff6(res, alpha, beta, z, acc);
        else
            fprintf(stderr, "None of the regions");
    }
    //num_print(res, true);
    
    /* num_set(res, _res); */
    /* num_clear(_res); */
}

/* Main function of the library */
int
mittleff_cmplx (double* res,
          const double alpha,
          const double beta,
          const double x, const double y,
          const double acc)
{
    log_trace("[%s] ==============================", __func__);
    log_trace("[%s] alpha=%g, beta=%g, z=%g%+g, acc=%g", __func__, alpha, beta, x, y, acc);
    assert(alpha > 0);

    num_t _res, _alpha, _beta, _z, _acc;

    _alpha = num_init(), _beta = num_init(),_z = num_init(), _acc = num_init();
    _res = num_init();
    
    num_set_d(_alpha, alpha), num_set_d(_beta, beta);
    num_set_d_d(_z, x, y), num_set_d(_acc, acc);

    _mittleff(_res, _alpha, _beta, _z, _acc);
    //num_print(_res, true);
    
    //log_trace("[%s] Computed", __func__);
    num_to_d_d(res, _res);
    
    num_clear(_res);
    num_clear(_alpha);
    num_clear(_beta);
    num_clear(_z);
    num_clear(_acc);
    log_trace("[%s] Done", __func__);
    log_trace("[%s] ==============================", __func__);

    return EXIT_SUCCESS;
}
