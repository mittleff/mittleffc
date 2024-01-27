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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <assert.h>

#include "num.h"
#include "new.h"
#include "partition.h"
#include "algorithm.h"

#include <gsl/gsl_math.h>
#include <stdbool.h>

#ifdef DEBUG
#include "log.h"
#endif

static void
_mittleff (num_t res,
           const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t tol)
{
#ifdef DEBUG
    log_info("\n[\033[1;33m%s\033[0m] Calling with parameters:\n\t    \033[1;32malpha\033[0m = %g\n\t    \033[1;32mbeta\033[0m  = %g\n\t    \033[1;32mz\033[0m = %+.14e%+.14e*I\n\t    \033[1;32mtol\033[0m = %g\n",
             __func__,
             num_to_d(alpha), num_to_d(beta), num_real_d(z), num_imag_d(z), num_to_d(tol));
#endif
    /* Test special cases */
    if (num_is_zero(z)) /* z = 0 */
    {
        num_rgamma(res, beta);
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m] \033[1;31mSpecial Case: z = 0.0\033[0m\n\t    \033[1;32mres\033[0m = 1/Gamma(%g) = %+.14e%+.14e*I\n",
                 __func__, num_to_d(beta), num_real_d(res), num_imag_d(res));
#endif        
    }
    else if (num_eq_d(alpha, 1.0) && num_eq_d(beta, 1.0)) /* exp(z) */
    {
        num_exp(res, z);
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m] \033[1;31mSpecial Case: alpha = beta = 1\033[0m\n\t    \033[1;32mres\033[0m = exp(z) = %+.14e%+.14e*I\n",
                 __func__, num_real_d(res), num_imag_d(res));
#endif        
    }
    else if (num_eq_d(alpha, 2.0) && num_eq_d(beta, 1.0)) /* cosh(sqrt(z)) */
    {
        num_sqrt(res, z);
        num_cosh(res, res);
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m] \033[1;31mSpecial Case: alpha = 2, beta = 1\033[0m\n\t    \033[1;32mres\033[0m = cosh(sqrt(z)) = %+.14e%+.14e*I\n",
                 __func__, num_real_d(res), num_imag_d(res));
#endif         
    }
    else if (num_eq_d(alpha, 0.5) && num_eq_d(beta, 1.0)) /* exp(z^2)*erfc(-z) */
    {
        num_t exp_z2, erfc_z;
        exp_z2 = new(num), erfc_z = new(num);
        num_pow_d(exp_z2, z, 2.0);
        num_exp(exp_z2, exp_z2);
        /* log_trace("[%s] exp(z^2) = %g+%g*I", */
        /*           __func__, */
        /*           num_to_complex(exp_z2)); */
        
        num_neg(erfc_z, z);
        num_erfc(erfc_z, erfc_z);
        /* log_trace("[%s] erfc(-z) = %g+%g*I", */
        /*           __func__, */
        /*           num_to_complex(erfc_z)); */
        
        num_mul(res, exp_z2, erfc_z);
        delete(exp_z2), delete(erfc_z);
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m] \033[1;31mSpecial Case: alpha = 0.5, beta = 1\033[0m\n\t    \033[1;32mres\033[0m = exp(z^2)*erfc(-z) = %+.14e%+.14e*I\n",
                 __func__, num_real_d(res), num_imag_d(res));
#endif         

        /* log_trace("[%s] alpha=0.5, beta=1, res = exp(z^2)*erfc(-z) = %g+%g*I", */
        /*           __func__, */
        /*           num_to_complex(res)); */
    }
    else if (num_eq_d(alpha, 2.0) && num_eq_d(beta, 2.0)) /* sinh(sqrt(z))/sqrt(z) */
    {
        num_t n, d;
        n = new(num);
        d = new(num);
        num_sqrt(n, z);
        num_sqrt(d, z);
        num_sinh(n, n);
        num_div(res, n, d);
        delete(n), delete(d);
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m] \033[1;31mSpecial Case: alpha = beta = 2\033[0m\n\t    \033[1;32mres\033[0m = sinh(sqrt(z))/sqrt(z) = %+.14e%+.14e*I\n",
                 __func__, num_real_d(res), num_imag_d(res));
#endif         
    }
    else if (in_region_G0(z))
    {
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m]\n\t z = %+.14e%+.14e*I located in region G0. Applying Taylor series.",
                 __func__, num_real_d(z), num_imag_d(z));
#endif        
        mittleff0(res, alpha, beta, z, tol);
    }
    else if (num_gt_d(alpha, 1.0)) /* apply recursive relation (2.2) */
    {
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m] \033[1;31mSpecial Case: alpha > 1\033[0m\n\t    Applying recursive relation",
                 __func__);
#endif         
        const int m = (int) (ceil((num_to_d(alpha) - 1)/2.0) + 1);

        num_t one_over_2mp1, alphap, zp, aux;

        one_over_2mp1 = new(num);
        alphap = new(num);
        zp = new(num), aux = new(num);

        num_set_d(one_over_2mp1, 1.0/(2.0 * m + 1.0));

        num_set_d(res, 0.0);
        for (int h = -m; h <= m; h++)
        {
            num_mul(alphap, alpha, one_over_2mp1);

            num_pow(zp, z, one_over_2mp1);

            num_set_d_d(aux, 0.0, 2.0 * M_PI * h * num_to_d(one_over_2mp1));
            num_exp(aux, aux);
            num_pow(zp, z, one_over_2mp1);
            num_mul(zp, zp, aux);

            _mittleff(aux, alphap, beta, zp, tol);

            num_add(res, res, aux);
        }
        num_mul(res, res, one_over_2mp1);

        delete(one_over_2mp1);
        delete(alphap);
        delete(zp), delete(aux);
        /* num_t sum, th, tmp, exp_th, newz; */
        /* sum = new(num), th = new(num), tmp = new(num), newz = new(num), exp_th = new(num); */
        /* num_zero(sum); */
        /* for (int h = -m; h <= m; h++) */
        /* { */
        /*     log_trace("[%s] h = %d", __func__, h); */
            
        /*     const double complex _newz = num_to_complex(newz); */

        /*     num_set_d_d(tmp, _res[0], _res[1]); */
        /*     num_add(sum, sum, tmp); */
        /* } */
        /* num_mul(sum, sum, one_over_2mp1); */
        /* num_set_num(res, sum); */
        /* delete(sum), delete(th), delete(tmp), delete(exp_th); */
    }
    else /* alpha <= 1 */
    {
        if (in_region_G1(z, alpha, tol))
        {
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m]\n\t z = %+.14e%+.14e*I located in region G1.",
                 __func__, num_real_d(z), num_imag_d(z));
#endif            
            mittleff1(res, alpha, beta, z, tol);
        }
        else if (in_region_G2(z, alpha, tol))
        {
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m]\n\t z = %+.14e%+.14e*I located in region G2.",
                 __func__, num_real_d(z), num_imag_d(z));
#endif            
            mittleff2(res, alpha, beta, z, tol);
        }
        else if (in_region_G3(z, alpha, tol))
        {
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m]\n\t z = %+.14e%+.14e*I located in region G3.",
                 __func__, num_real_d(z), num_imag_d(z));
#endif              
            mittleff3(res, alpha, beta, z, tol);
        }
        else if (in_region_G4(z, alpha, tol))
        {
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m]\n\t z = %+.14e%+.14e*I located in region G4.",
                 __func__, num_real_d(z), num_imag_d(z));
#endif             
            mittleff4(res, alpha, beta, z, tol);
        }
        else if (in_region_G5(z, alpha, tol))
        {
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m]\n\t z = %+.14e%+.14e*I located in region G5.",
                 __func__, num_real_d(z), num_imag_d(z));
#endif            
            mittleff5(res, alpha, beta, z, tol);
        }
        else if (in_region_G6(z, alpha, tol))
        {
#ifdef DEBUG
        log_info("\n[\033[1;33m%s\033[0m]\n\t z = %+.14e%+.14e*I located in region G6.",
                 __func__, num_real_d(z), num_imag_d(z));
#endif            
            mittleff6(res, alpha, beta, z, tol);
        }
        else
            fprintf(stderr, "None of the regions");
    }
    //num_print(res, true);
    
    /* num_set(res, _res); */
    /* delete(_res); */
}

/* Main function of the library */
int
mittleff_cmplx (double* res,
          const double alpha,
          const double beta,
          const double x, const double y,
          const double tol)
{
    #ifdef DEBUG
    log_info("\n[\033[1;33m%s\033[0m] Calling with parameters:\n\t\t \033[1;32malpha\033[0m = %g\n\t\t \033[1;32mbeta\033[0m  = %g\n\t\t \033[1;32mz\033[0m = %+.14e%+.14e*I\n\t\t \033[1;32mtol\033[0m = %g\n",
             __func__, alpha, beta, x, y, tol);
    #endif
    
    assert (alpha > 0.0);

    num_t _res, _alpha, _beta, _z, _tol;

    _alpha = new(num), _beta = new(num),_z = new(num), _tol = new(num);
    _res = new(num);
    
    num_set_d(_alpha, alpha), num_set_d(_beta, beta);
    num_set_d_d(_z, x, y), num_set_d(_tol, tol);

    _mittleff(_res, _alpha, _beta, _z, _tol);
    num_to_d_d(res, _res);
    
    delete(_res);
    delete(_alpha), delete(_beta);
    delete(_z), delete(_tol);

    return EXIT_SUCCESS;
}
