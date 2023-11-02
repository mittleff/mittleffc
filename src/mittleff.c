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
#include <complex.h>
#include <assert.h>

#include "log.h"

#include "mittleff.h"
#include "num/new.h"
#include "num/num.h"

static num_t
mittleff (const num_t alpha,
          const num_t beta,
          const num_t z,
          const num_t acc)
{
    log_trace("[%s] (alpha, beta, z, acc) = (%g, %g, %g%+gj, %g)",
              __func__,
              num_to_double(alpha),
              num_to_double(beta),
              num_to_complex(z), acc);
    num_t zero = new(num, 0.0, 0.0);
    num_t one = new(num, 1.0, 0.0);    
    num_t sum = num_rgamma(beta);
    num_t prev = num_rgamma(beta);
    num_t curr = zero;
    num_t k = one;
    while (true)
    {
        log_debug("[%s] k = %d,\t partial_sum = %+.8e", __func__,
                  (int)num_to_double(k),
                  num_to_double(sum));
        curr = num_mul(
            num_pow(z, k),
            num_rgamma(
                num_add(
                    num_mul(alpha, k),
                    beta)));
        sum = num_add(sum, curr);
        if (num_le(num_abs(num_sub(curr, prev)), acc))
            break;
        prev = curr;
        k = num_add(k, one);
    }
    log_trace("[%s] Summed %d terms: %+.8e", __func__,
                  (int)num_to_double(k),
                  num_to_double(sum));
    
    delete(zero); delete(one); delete(prev); delete(curr); delete(k);
    return sum;
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
    
    num_t _alpha = new(num, alpha, 0.0);
    num_t _beta  = new(num, beta,  0.0);
    num_t _z     = new(num, re_z, im_z);
    num_t _acc   = new(num, acc, 0.0);
    const double complex ml = num_to_complex(mittleff(_alpha, _beta, _z, _acc));
    res[0] = creal(ml); res[1] = cimag(ml);
    delete(_alpha); delete(_beta); delete(_z); delete(_acc);
    
    return EXIT_SUCCESS;
}
