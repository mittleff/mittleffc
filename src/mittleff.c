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
    const double _alpha = num_to_double(alpha);
    const double _beta = num_to_double(beta);
    const double _acc = num_to_double(acc);
         
    log_trace("[%s] (alpha, beta, z, acc) = (%g, %g, %g%+gj, %g)",
              __func__,
              _alpha,
              _beta,
              num_real_d(z),
              num_imag_d(z),
              _acc);

    num_t zero = new(num, 0.0, 0.0);
    num_t one = new(num, 1.0, 0.0);
    num_t two = new(num, 2.0, 0.0);
    num_t res = zero;
    
    if (in_region_G0(z))
    {
        log_info("[%s] z=(%+.5e, %+.5e) in region G0. Applying Taylor series",
                 __func__, num_real_d(z), num_imag_d(z));
        res = mittleff0(alpha, beta, z, acc);
    }
    else if (num_gt(alpha, one))
    {
        // apply recursive relation (2.2)
        log_info("[%s] applying recursive relation", __func__);
        const int m = (int) (ceil((_alpha - 1)/2.0) + 1);
        num_t m_num = new(num, (double)m, 0.0);

        num_t sum = zero;
        for (int h = -m; h <= m; h++)
        {
            num_t znew = num_mul(num_pow(z, num_div(one, num_add(num_mul(two, m_num), one))),
                                 num_exp(num_div(new(num, 0.0, 2.0 * M_PI * h), num_add(num_mul(two, m_num), one))));
            num_t ml = mittleff(new(num, (double)_alpha/(2*m + 1), 0.0), beta, znew, acc);
            sum = num_add(sum, ml);
            log_debug("[%s] %g %g", __func__, num_real_d(sum), num_imag_d(sum));
        }
        delete(m_num);
        log_info("[%s] %g %g", __func__, num_real_d(sum), num_imag_d(sum));

        return num_mul(num_div(one, num_add(num_mul(two, m_num), one)), sum);


        
        /* m_num = Numeric(m, 0.0) */
        /* s = zero */
        /* for h in range(-m, m+1):                        */
        /*     znew = num_to_complex(num_mul(num_pow(z, num_div(one, num_add(num_mul(two, m_num), one))), */
        /*             num_exp(num_div(Numeric(0.0, 2.0 * M_PI * h), num_add(num_mul(two, m_num), one))))) */
        /*     ml = _wrap_mittleff(_alpha/(2*m + 1), _beta, znew, _tol) */
        /*     s = num_add(s, ml)  */
        /* return num_mul(num_div(one, num_add(num_mul(two, m_num), one)), s) */
        /* return new(num, 1.0, 0.0); */
    }
    else /* alpha <= 1 */
    {
        log_info("[%s] applying main algorithm", __func__);

        /* compute parameters */

         /* double C0 = pow(1.3, 1 - num_to_double(p.alpha))/(M_PI * sin(M_PI * num_to_double(p.alpha))); /\* Equation (5.3) *\/ */
        /* p.r1 = new(num, pow(-2.0 * log(num_to_double(tol)/C0), num_to_double(p.alpha)), 0.0); /\* Equation (4.21) *\/ */

        if (in_region_G1(z, alpha, acc))
        {
            log_info("[%s] z=(%+.5e, %+.5e) in region G1",
                     __func__, num_real_d(z), num_imag_d(z));

            return new(num, 1.0, 0.0);
        }
        else if (in_region_G2(z, alpha, acc))
        {
            log_info("[%s] z=(%+.5e, %+.5e) in region G2",
                     __func__, num_real_d(z), num_imag_d(z));

            return new(num, 2.0, 0.0);
        }
        else if (in_region_G3(z, alpha, acc))
        {
            log_info("[%s] z=(%+.5e, %+.5e) in region G3",
                     __func__, num_real_d(z), num_imag_d(z));

            return new(num, 3.0, 0.0);
        }
        else if (in_region_G4(z, alpha, acc))
        {
            log_info("[%s] z=(%+.5e, %+.5e) in region G4",
                     __func__, num_real_d(z), num_imag_d(z));

            return new(num, 4.0, 0.0);
        }
        else if (in_region_G5(z, alpha, acc))
        {
            log_info("[%s] z=(%+.5e, %+.5e) in region G5",
                     __func__, num_real_d(z), num_imag_d(z));

            return new(num, 5.0, 0.0);
        }
        else if (in_region_G6(z, alpha, acc))
        {
            log_info("[%s] z=(%+.5e, %+.5e) in region G6",
                     __func__, num_real_d(z), num_imag_d(z));

            return new(num, 6.0, 0.0);
        }
        else
        {
            log_info("[%s] z=(%+.5e, %+.5e) in no region [G0, G6]",
                     __func__, num_real_d(z), num_imag_d(z));

            return new(num, -1.0, 0.0);
        }
        
    }

    delete(zero); delete(one);

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
    
    num_t _alpha = new(num, alpha, 0.0);
    num_t _beta  = new(num, beta,  0.0);
    num_t _z     = new(num, re_z, im_z);
    num_t _acc   = new(num, acc, 0.0);
    num_t ml = mittleff(_alpha, _beta, _z, _acc);    
    res[0] = num_real_d(ml); res[1] = num_imag_d(ml);
    delete(_alpha); delete(_beta); delete(_z); delete(_acc);
    
    return EXIT_SUCCESS;
}
