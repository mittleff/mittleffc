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

#include "log.h"

#include "mittleff.h"
#include "num.h"
#include "new.h"
#include "partition.h"
#include "algorithm.h"

static num_t
mittleff (const num_t alpha,
          const num_t beta,
          const num_t z,
          const num_t acc)
{
    const double _alpha = num_to_d(alpha);
    const double _beta = num_to_d(beta);
    const double _acc = num_to_d(acc);
         
    /* log_trace("[%s] (alpha, beta, z, acc) = (%g, %g, %g%+gj, %g)", */
    /*           __func__, */
    /*           _alpha, */
    /*           _beta, */
    /*           num_real_d(z), */
    /*           num_imag_d(z), */
    /*           _acc); */

    num_t zero = num_from_d(0.0);
    num_t one = num_from_d(1.0);
    num_t two = num_from_d(2.0);
    num_t res = num_from_d(0.0);
    
    if (in_region_G0(z))
    {
        /* log_info("[%s] z=(%+.5e, %+.5e) in region G0. Applying Taylor series", */
        /*          __func__, num_real_d(z), num_imag_d(z)); */
        res = mittleff0(alpha, beta, z, acc);
    }
    else if (num_gt(alpha, one))
    {
        // apply recursive relation (2.2)
        log_info("[%s] applying recursive relation", __func__);
        const int m = (int) (ceil((_alpha - 1)/2.0) + 1);
        const num_t m_num = new(num, (double)m, 0.0);
        num_t sum = zero;
        for (int h = -m; h <= m; h++)
        {
            num_t znew = num_mul(num_pow(z, num_div(one, num_add(num_mul(two, m_num), one))),
                                 num_exp(num_div(new(num, 0.0, 2.0 * M_PI * h), num_add(num_mul(two, m_num), one))));
            num_t ml = mittleff(new(num, (double)_alpha/(2*m + 1), 0.0), beta, znew, acc);
            sum = num_add(sum, ml);
        }
        res = num_mul(num_div(one, num_add(num_mul(two, m_num), one)), sum);
        delete(m_num);
    }
    else /* alpha <= 1 */
    {
        log_info("[%s] applying main algorithm", __func__);

        if (in_region_G1(z, alpha, acc))
        {
            /* log_info("[%s] z=(%+.5e, %+.5e) in region G1", */
            /*          __func__, num_real_d(z), num_imag_d(z)); */

            res = mittleff1(alpha, beta, z, acc);
        }
        else if (in_region_G2(z, alpha, acc))
        {
            /* log_info("[%s] z=(%+.5e, %+.5e) in region G2", */
            /*          __func__, num_real_d(z), num_imag_d(z)); */

            res = mittleff2(alpha, beta, z, acc);
        }
        else if (in_region_G3(z, alpha, acc))
        {
            /* log_info("[%s] z=(%+.5e, %+.5e) in region G3", */
            /*          __func__, num_real_d(z), num_imag_d(z)); */

            res = mittleff3(alpha, beta, z, acc);
        }
        else if (in_region_G4(z, alpha, acc))
        {
            /* log_info("[%s] z=(%+.5e, %+.5e) in region G4", */
            /*          __func__, num_real_d(z), num_imag_d(z)); */

            res = mittleff4(alpha, beta, z, acc);
        }
        else if (in_region_G5(z, alpha, acc))
        {
            /* log_info("[%s] z=(%+.5e, %+.5e) in region G5", */
            /*          __func__, num_real_d(z), num_imag_d(z)); */

            res = mittleff5(alpha, beta, z, acc);
        }
        else if (in_region_G6(z, alpha, acc))
        {
            /* log_info("[%s] z=(%+.5e, %+.5e) in region G6", */
            /*          __func__, num_real_d(z), num_imag_d(z)); */

            res = mittleff6(alpha, beta, z, acc);
        }
        else
        {
            /* log_info("[%s] z=(%+.5e, %+.5e) in no region [G0, G6]", */
            /*          __func__, num_real_d(z), num_imag_d(z)); */

            res = num_from_d(-1.0);
        }
        
    }

    delete(zero); delete(one); delete(two);

    return res;
}

/* Main function of the library */
int
mittleff_cmplx (double* res,
                const double alpha,
                const double beta,
                const double re_z, const double im_z,
                const double acc)
{
    assert(alpha > 0);
    
    num_t _alpha = num_from_d(alpha);
    num_t _beta  = num_from_d(beta);
    num_t _z     = num_from_d_d(re_z, im_z);
    num_t _acc   = num_from_d(acc);
    
    num_t ml = mittleff(_alpha, _beta, _z, _acc);
    double complex z = num_to_complex(ml);
    
    res[0] = creal(z);
    res[1] = cimag(z);
    
    delete(_alpha); delete(_beta);
    delete(_z); delete(_acc);
    delete(ml);

    return EXIT_SUCCESS;
}
