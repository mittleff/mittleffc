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

/* #include "log.h" */

/* #include "mittleff.h" */
#include "num.h"
#include "new.h"
/* #include "partition.h" */
/* #include "algorithm.h" */

void
mittleff0 (double* res,
           const double a,
           const double b,
           const double x,
           const double y,
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

static void
mittleff (double* res,
          const double alpha,
          const double beta,
          const double x, const double y,
          const double acc)
{
    double val[2] = { 0.0, 0.0 };
    
    mittleff0(val, alpha, beta, x, y, acc);

    res[0] = val[0], res[1] = val[1];
    /* num_t zero = num_from_d(0.0); */
    /* num_t one = num_from_d(1.0); */
    /* num_t two = num_from_d(2.0); */
    /* num_t res = num_from_d(0.0); */
    
    /* if (in_region_G0(z)) */
    /* { */
    /*     res = mittleff0(alpha, beta, z, acc); */
    /* } */
    /* else if (num_gt(alpha, one)) */
    /* { */
    /*     // apply recursive relation (2.2) */
    /*     log_info("[%s] applying recursive relation", __func__); */
    /*     const int m = (int) (ceil((_alpha - 1)/2.0) + 1); */
    /*     const num_t m_num = new(num, (double)m, 0.0); */
    /*     num_t sum = zero; */
    /*     for (int h = -m; h <= m; h++) */
    /*     { */
    /*         num_t znew = num_mul(num_pow(z, num_div(one, num_add(num_mul(two, m_num), one))), */
    /*                              num_exp(num_div(new(num, 0.0, 2.0 * M_PI * h), num_add(num_mul(two, m_num), one)))); */
    /*         num_t ml = mittleff(new(num, (double)_alpha/(2*m + 1), 0.0), beta, znew, acc); */
    /*         sum = num_add(sum, ml); */
    /*     } */
    /*     res = num_mul(num_div(one, num_add(num_mul(two, m_num), one)), sum); */
    /*     delete(m_num); */
    /* } */
    /* else /\* alpha <= 1 *\/ */
    /* { */
    /*     log_info("[%s] applying main algorithm", __func__); */

    /*     if (in_region_G1(z, alpha, acc)) */
    /*     { */
    /*         /\* log_info("[%s] z=(%+.5e, %+.5e) in region G1", *\/ */
    /*         /\*          __func__, num_real_d(z), num_imag_d(z)); *\/ */

    /*         res = mittleff1(alpha, beta, z, acc); */
    /*     } */
    /*     else if (in_region_G2(z, alpha, acc)) */
    /*     { */
    /*         /\* log_info("[%s] z=(%+.5e, %+.5e) in region G2", *\/ */
    /*         /\*          __func__, num_real_d(z), num_imag_d(z)); *\/ */

    /*         res = mittleff2(alpha, beta, z, acc); */
    /*     } */
    /*     else if (in_region_G3(z, alpha, acc)) */
    /*     { */
    /*         /\* log_info("[%s] z=(%+.5e, %+.5e) in region G3", *\/ */
    /*         /\*          __func__, num_real_d(z), num_imag_d(z)); *\/ */

    /*         res = mittleff3(alpha, beta, z, acc); */
    /*     } */
    /*     else if (in_region_G4(z, alpha, acc)) */
    /*     { */
    /*         /\* log_info("[%s] z=(%+.5e, %+.5e) in region G4", *\/ */
    /*         /\*          __func__, num_real_d(z), num_imag_d(z)); *\/ */

    /*         res = mittleff4(alpha, beta, z, acc); */
    /*     } */
    /*     else if (in_region_G5(z, alpha, acc)) */
    /*     { */
    /*         /\* log_info("[%s] z=(%+.5e, %+.5e) in region G5", *\/ */
    /*         /\*          __func__, num_real_d(z), num_imag_d(z)); *\/ */

    /*         res = mittleff5(alpha, beta, z, acc); */
    /*     } */
    /*     else if (in_region_G6(z, alpha, acc)) */
    /*     { */
    /*         /\* log_info("[%s] z=(%+.5e, %+.5e) in region G6", *\/ */
    /*         /\*          __func__, num_real_d(z), num_imag_d(z)); *\/ */

    /*         res = mittleff6(alpha, beta, z, acc); */
    /*     } */
    /*     else */
    /*     { */
    /*         /\* log_info("[%s] z=(%+.5e, %+.5e) in no region [G0, G6]", *\/ */
    /*         /\*          __func__, num_real_d(z), num_imag_d(z)); *\/ */

    /*         res = num_from_d(-1.0); */
    /*     } */
        
    /* } */

    /* delete(zero); delete(one); delete(two); */

    /* return res; */
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
