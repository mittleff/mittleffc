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
#include "partition/partition.h"
#include "algorithm/algorithm.h"

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
              num_real_d(z),
              num_imag_d(z),
              acc);

    num_t one = new(num, 1.0, 0.0);
    if (in_region_G0(z))
    {
        log_info("[%s] z=(%+.5e, %+.5e) in region G0. Applying Taylor series",
                 __func__, num_real_d(z), num_imag_d(z));
        return mittleff0(alpha, beta, z, acc);
    }
    else if (num_gt(alpha, one))
    {
        return new(num, 1.0, 0.0);
    }
    else
    {
        return new(num, 0.0, 0.0);
    }
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
    num_t ml = mittleff(_alpha, _beta, _z, _acc);    
    res[0] = num_real_d(ml); res[1] = num_imag_d(ml);
    delete(_alpha); delete(_beta); delete(_z); delete(_acc);
    
    return EXIT_SUCCESS;
}
