/*
 * This file is part of integration.c (https://github.com/padawanphysicist/integration.c).
 *
 * integration.c is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * integration.c is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * integration.c. If not, see <https://www.gnu.org/licenses/>.
 */

/** 
 * @file quad.c
 * @brief Implementation of the main routines.
 */

#include "quad.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "new.h"
#include "num.h"

#include "utils.h"

#include "flint/mag.h"
#include "flint/acb.h"
#include "flint/acb_calc.h"

static double
arbtod(arb_t x)
{
    return arf_get_d(arb_midref(x), ARF_RND_NEAR);
}


static void
acbtonum (num_t res, acb_t x)
{
    arb_t re, im;
    arb_init(re), arb_init(im);
    acb_get_real(re, x), acb_get_imag(im, x);
    num_set_d_d (res, arbtod(re), arbtod(im));
    arb_clear(re), arb_clear(im);
}



/* f(z) = sin(z) */
int
f_integrand (acb_ptr res, const acb_t z, void * params, slong order, slong prec)
{
    //log_trace("calling wrapper function?");
    
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    //struct my_f_params * params = (struct my_f_params *)p;
    num_function_t * F = (num_function_t *) params;

    num_t x, y;
    x = new(num), y = new(num);
    acbtonum(x, z);
    (F->function)(y, x, (F->params)); // TODO
    
    acb_set_d_d(res, num_real_d (y), num_imag_d (y));

    delete(x), delete(y);

    return 0;
}


static void
quad_gauss_legendre (num_t res,
                     num_function_t F,
                     const num_t from,
                     const num_t to)
{
    //log_trace("calling gauss legendre integration");
    acb_t ret, a, b, aux;
    mag_t tol;
    slong prec, goal;
    acb_calc_integrate_opt_t options;

    acb_calc_integrate_opt_init(options);

    prec = 63, goal = prec;

    acb_init(ret), acb_init(a), acb_init(b), acb_init(aux);
    mag_init(tol);

    mag_set_d(tol, 1.0e-15);

    /* Set integration limits */
    acb_set_d_d(a, num_real_d(from), num_imag_d(from));
    acb_set_d_d(b, num_real_d(to), num_imag_d(to));

    //log_trace("before integrate with ARB");
    
    acb_calc_integrate(ret, f_integrand, &F, a, b, goal, tol, options, prec);
    /* log_trace("after integrate with ARB"); */
    acbtonum (res, ret);
       
    acb_clear(ret), acb_clear(a), acb_clear(b), acb_clear(aux);
    mag_clear(tol);
}

void
quad (num_t res,
      num_function_t F,
      const num_t from,
      const num_t to,
      const int method)
{
    UNUSED(method);
    /* log_trace("before call routine integration"); */
    quad_gauss_legendre (res, F, from, to);

    /* switch (method) { */
    /* case 0: */
    /*     // tanh-sinh quadrature */
    /*     num_t tolerance, est_err, num_eval; */
    /*     tolerance = new(num), est_err = new(num), num_eval = new(num); */
    /*     num_set_d(tolerance, 1.0e-15); */
    /*     quad_tanhsinh_quad(res, f, ctx, from, to, tolerance, est_err, num_eval); */
    /*     delete(tolerance), delete(est_err), delete(num_eval); */
    /*     break; */
    /* default: */
    /*     num_t n; */
    /*     n = new(num); */
    /*     num_set_d(n, 100.0); */
    /*     integration_qsimp(res, f, ctx, from, to, n); */
    /*     delete(n); */
    /* } */
}
