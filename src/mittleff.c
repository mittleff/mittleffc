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
/* #include "mittleff.h" */
#include "num.h"
#include "new.h"
#include "partition.h"
#include "algorithm.h"

#include <gsl/gsl_math.h>
#include <stdbool.h>

#define MAX2(a,b) ((a>b)?(a):(b))
#define MAX3(a,b,c) MAX2(a,MAX2(b,c))

#define GSL_EPSILON_FCMP 1e-8
static bool
gt (const double x, const double y)
{
    int ret = gsl_fcmp(x, y, GSL_EPSILON_FCMP);
    return (ret == 1) ? true : false;
}


static void
mittleff (double* res,
          const double alpha,
          const double beta,
          const double x, const double y,
          const double acc)
{
    log_trace("[%s] (alpha, beta, z, acc) = %+.5e %+.5e (%+.5e%+.5ej) %.5e", __func__, alpha, beta, x, y, acc);
    double val[2] = { 0.0, 0.0 };
    
    if (in_region_G0(x, y))
    {
        log_trace("[%s] z = (%+.5e%+.5ej) in G_0: |z| = %g", __func__, x, y, sqrt(x*x + y*y));
        mittleff0(val, alpha, beta, x, y, acc);
    }
    else if (gt(alpha, 1.0)) /* apply recursive relation (2.2) */
    {
        log_trace("[%s] alpha = %+.5e > 1, applying recursive relation", __func__, alpha);
        double _res[2] = { 0.0, 0.0 };
        const int m = (int) (ceil((alpha - 1)/2.0) + 1);
        const double one_over_2mp1 = 1.0/(2.0 * m + 1.0);

        num_t sum, th, z, tmp, exp_th, newz;
        sum = new(num), th = new(num), z = new(num), tmp = new(num), newz = new(num), exp_th = new(num);
        num_set_d_d(z, x, y);
        num_zero(sum);
        for (int h = -m; h <= m; h++)
        {
            log_trace("[%s] h = %d", __func__, h);
            num_set_d_d(th, 0.0, 2.0 * M_PI * h * one_over_2mp1);
            num_exp(exp_th, th);
            num_pow_d(newz, z, one_over_2mp1);
            num_mul(newz, newz, exp_th);            
            const double complex _newz = num_to_complex(newz);
            mittleff(_res, alpha * one_over_2mp1, beta, creal(_newz), cimag(_newz), acc);
            num_set_d_d(tmp, _res[0], _res[1]);
            num_add(sum, sum, tmp);
        }
        num_mul_d(sum, sum, 1.0/(2.0*m + 1.0));
        num_to_d_d(val, sum);
        log_trace("res = %g %g", val[0], val[1]);
        delete(sum), delete(th), delete(z), delete(tmp), delete(exp_th);
    }
    else /* alpha <= 1 */
    {
        if (in_region_G1(x, y, alpha, acc))
        {
            log_trace("[%s] z = (%+.5e%+.5ej) in G_1", __func__, x, y);
            mittleff1(val, alpha, beta, x, y, acc);
        }
        else if (in_region_G2(x, y, alpha, acc))
        {
            log_trace("region 2\n");
            /* res = mittleff2(alpha, beta, z, acc); */
        }
        else if (in_region_G3(x, y, alpha, acc))
        {            
            log_trace("region 3\n");
            /* res = mittleff3(alpha, beta, z, acc); */
        }
        else if (in_region_G4(x, y, alpha, acc))
        {
            log_trace("region 4\n");

            /* res = mittleff4(alpha, beta, z, acc); */
        }
        else if (in_region_G5(x, y, alpha, acc))
        {
            log_trace("[%s] z = (%+.5e%+.5ej) in G_5", __func__, x, y);
            mittleff5(val, alpha, beta, x, y, acc);
        }
        else if (in_region_G6(x, y, alpha, acc))
        {
            log_trace("[%s] z = (%+.5e%+.5ej) in G_6", __func__, x, y);
            mittleff6(val, alpha, beta, x, y, acc);
        }
        else
        {
            log_trace("None of the regions");
            /* res = num_from_d(-1.0); */
        }
        
    }

    /* delete(zero); delete(one); delete(two); */
    res[0] = val[0], res[1] = val[1];
    //log_trace("res = %g %g", res[0], res[1]);
}

/* Main function of the library */
int
mittleff_cmplx (double* res,
                const double alpha,
                const double beta,
                const double x, const double y,
                const double acc)
{
    assert(alpha > 0);
    
    mittleff(res, alpha, beta, x, y, acc);

    return EXIT_SUCCESS;
}
