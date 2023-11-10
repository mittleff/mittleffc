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
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <assert.h>

/* #include "mittleff.h" */
#include "num.h"
#include "new.h"
#include "partition.h"
#include "algorithm.h"

#include <gsl/gsl_math.h>
#include <stdbool.h>

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
    double val[2] = { 0.0, 0.0 };
    
    if (in_region_G0(x, y))
    {
        mittleff0(val, alpha, beta, x, y, acc);
    }
    else if (gt(alpha, 1.0))
    {
        /* /\* apply recursive relation (2.2) *\/ */
        /* const int m = (int) (ceil((alpha - 1)/2.0) + 1); */
        /* const num_t m_num = new(num, (double)m, 0.0); */
        /* num_t sum = zero; */
        /* for (int h = -m; h <= m; h++) */
        /* { */
        /*     num_t znew = num_mul(num_pow(z, num_div(one, num_add(num_mul(two, m_num), one))), */
        /*                          num_exp(num_div(new(num, 0.0, 2.0 * M_PI * h), num_add(num_mul(two, m_num), one)))); */
        /*     num_t ml = mittleff(new(num, (double)_alpha/(2*m + 1), 0.0), beta, znew, acc); */
        /*     sum = num_add(sum, ml); */
        /* } */
        /* res = num_mul(num_div(one, num_add(num_mul(two, m_num), one)), sum); */
        /* delete(m_num); */
    }
    else /* alpha <= 1 */
    {
        if (in_region_G1(x, y, alpha, acc))
        {
            /* res = mittleff1(alpha, beta, z, acc); */
        }
        else if (in_region_G2(x, y, alpha, acc))
        {          
            /* res = mittleff2(alpha, beta, z, acc); */
        }
        else if (in_region_G3(x, y, alpha, acc))
        {            

            /* res = mittleff3(alpha, beta, z, acc); */
        }
        else if (in_region_G4(x, y, alpha, acc))
        {
            

            /* res = mittleff4(alpha, beta, z, acc); */
        }
        else if (in_region_G5(x, y, alpha, acc))
        {
            

            /* res = mittleff5(alpha, beta, z, acc); */
        }
        else if (in_region_G6(x, y, alpha, acc))
        {
            

            /* res = mittleff6(alpha, beta, z, acc); */
        }
        else
        {
            

            /* res = num_from_d(-1.0); */
        }
        
    }

    /* delete(zero); delete(one); delete(two); */

    res[0] = val[0], res[1] = val[1];
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
