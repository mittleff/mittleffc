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
#include "new.h"
#include "partition.h"
#include "algorithm.h"

#include <gsl/gsl_math.h>
#include <stdbool.h>

#define MAX2(a,b) ((a>b)?(a):(b))
#define MAX3(a,b,c) MAX2(a,MAX2(b,c))

static void
_mittleff (num_t res,
           const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc)
{
    log_trace("[%s] alpha=%g, beta=%g, z=%g%+g, acc=%g", __func__,
              num_to_d(alpha), num_to_d(beta), num_real_d(z), num_imag_d(z), num_to_d(acc));
    if (in_region_G0(z))
        mittleff0(res, alpha, beta, z, acc);
    else if (num_gt_d(alpha, 1.0)) /* apply recursive relation (2.2) */
    {
        /* log_trace("[%s] alpha = %+.5e > 1, applying recursive relation", __func__, alpha); */
        /* const int m = (int) (ceil((alpha - 1)/2.0) + 1); */
        /* const double one_over_2mp1 = 1.0/(2.0 * m + 1.0); */

        /* num_t sum, th, tmp, exp_th, newz; */
        /* sum = new(num), th = new(num), z = new(num), tmp = new(num), newz = new(num), exp_th = new(num); */
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
        /* delete(sum), delete(th), delete(tmp), delete(exp_th); */
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
    num_print(res, true);
    log_trace("[%s] Done", __func__);
    /* num_set(res, _res); */
    /* delete(_res); */
}

/* Main function of the library */
int
mittleff_cmplx (double* res,
          const double alpha,
          const double beta,
          const double x, const double y,
          const double acc)
{
    log_trace("[%s] alpha=%g, beta=%g, z=%g%+g, acc=%g", __func__, alpha, beta, x, y, acc);
    assert(alpha > 0);

    num_t _res, _alpha, _beta, _z, _acc;

    _alpha = new(num), _beta = new(num),_z = new(num), _acc = new(num);
    _res = new(num);
    
    num_set_d(_alpha, alpha), num_set_d(_beta, beta);
    num_set_d_d(_z, x, y), num_set_d(_acc, acc);

    _mittleff(_res, _alpha, _beta, _z, _acc);
    //num_print(_res, true);
    
    //log_trace("[%s] Computed", __func__);
    num_to_d_d(res, _res);
    
    delete(_res);
    delete(_alpha);
    delete(_beta);
    delete(_z);
    delete(_acc);
    log_trace("[%s] Done", __func__);

    return EXIT_SUCCESS;
}
