/*
 * This file is part of mittleffc (https://github.com/mittleff/mittleffc).
 *
 * mittleffc is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * mittleffc is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * mittleffc. If not, see <https://www.gnu.org/licenses/>.
 */

/** 
 * @file integrate.c
 * @brief Implementation for the routines integrating the functions on complex plane.
 */
#include "integrate.h"

#include <flintutils.h>

#include <flint/acb_calc.h>

#include "types.h"

#include <math.h>
#include <assert.h>

/* Parameters for integrating B */
typedef struct {
    arb_t alpha;
    arb_t beta;
    acb_t z;
    arb_t c;
} integ_params_t;

void
fn_A (acb_t res,
   const acb_t z,
   const arb_t alpha,
   const arb_t beta,
   const arb_t x)
{
    /* (1.0/alpha) * z**((1.0 - beta)/alpha) * exp(z**(1.0/alpha) * cos(x/alpha)) */

    acb_t fac1, fac2, fac3, aux;

    acb_init(fac1);
    acb_init(fac2);
    acb_init(fac3);
    acb_init(aux);

    acb_set_arb(fac1, alpha);

    acb_set_d(fac2, 1.0);
    acb_sub_arb(fac2, fac2, beta, PREC);
    acb_div_arb(fac2, fac2, alpha, PREC);
    acb_pow(fac2, z, fac2, PREC);

    acb_pow(fac3, z, fac1, PREC);
    acb_set_arb(aux, x);
    acb_div_arb(aux, aux, alpha, PREC);
    acb_cos(aux, aux, PREC);
    acb_mul(fac3, fac3, aux, PREC);
    acb_exp(fac3, fac3, PREC);

    acb_mul(res, fac1, fac2, PREC);
    acb_mul(res, res, fac3, PREC);

    acb_clear(fac1);
    acb_clear(fac2);
    acb_clear(fac3);
    acb_init(aux);
}

static void
omega(acb_t res,
      const arb_t x,
      const arb_t y,
      const arb_t alpha,
      const arb_t beta)
{
    /* x**(1.0/alpha) * sin(y/alpha) + y * (1.0 + (1.0 - beta)/alpha) */
    arb_t _res, fac1, fac2, aux1, aux2;

    arb_init(fac1);
    arb_init(fac2);
    arb_init(aux1);
    arb_init(aux2);
    arb_init(_res);
    
    arb_inv(aux2, alpha, PREC);

    arb_mul(aux1, aux2, y, PREC);
    arb_sin(aux1, aux1, PREC);
    arb_set(fac1, x);
    arb_pow(fac1, fac1, aux2, PREC);
    arb_mul(fac1, fac1, aux1, PREC);

    arb_sub(aux1, aux1, beta, PREC);
    arb_div(aux1, aux1, alpha, PREC);
    arb_set_d(fac2, 1.0);
    arb_add(fac2, fac2, aux1, PREC);
    arb_mul(fac2, fac2, y, PREC);

    arb_add(_res, fac1, fac2, PREC);
    acb_set_arb(res, _res);

    arb_clear(fac1);
    arb_clear(fac2);
    arb_clear(aux1);
    arb_clear(aux2);
    arb_clear(_res);
}

int
fn_B (acb_ptr res, const acb_t r, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */
    
    integ_params_t* p = (integ_params_t*) param;

    arb_t pi, rr;
    acb_t w, aval, fac1, fac2, aux, aux2, phi;

    acb_init(w);
    arb_init(pi);
    acb_init(aval);
    acb_init(fac1);
    acb_init(fac2);
    acb_init(aux);
    acb_init(aux2);
    arb_init(rr);
    acb_init(phi);

    acb_set_arb(phi, p->c);

    fn_A (aval, r, p->alpha, p->beta, p->c /* phi */);
    acb_get_real(rr, r);
    omega(w, rr, p->c /* phi */, p->alpha, p->beta);

    acb_sub(aux, w, phi, prec);
    acb_sin(aux, aux, prec);
    acb_mul(fac1, r, aux, prec);
    acb_sin(aux, w, prec);
    acb_mul(aux, aux, p->z, prec);
    acb_sub(fac1, fac1, aux, prec);

    acb_sqr(fac2, r, prec);
    acb_sqr(aux, p->z, prec);
    acb_add(fac2, fac2, aux, prec);
    acb_cos(aux2, phi, prec);
    acb_mul(aux2, aux2, p->z, prec);
    acb_mul(aux2, aux2, r, prec);
    acb_set_d(aux, -2.0);
    acb_mul(aux, aux, aux2, prec);
    acb_add(fac2, fac2, aux, prec);

    arb_const_pi(pi, prec);

    acb_set_arb(res, pi);
    acb_inv(res, res, prec);
    acb_mul(res, res, aval, prec);
    acb_mul(res, res, fac1, prec);
    acb_div(res, res, fac2, prec);
    
    acb_clear(w);
    arb_clear(pi);
    acb_clear(aux);
    acb_clear(aux2);
    acb_clear(fac1);
    acb_clear(fac2);
    arb_clear(rr);
    acb_clear(phi);

    return 0;
}

/* void */
/* B (acb_t res, */
/*    const num_t r, */
/*    const num_t alpha, */
/*    const num_t beta, */
/*    const num_t z, */
/*    const num_t phi) */
/* { */
/*     num_t n, d, fac1, fac2; */

/*     n = new(num), d = new(num); */
/*     fac1 = new(num), fac2 = new(num); */

/*     /\* numerator *\/ */
/*     omega(fac1, r, phi, alpha, beta); */
/*     num_sub(fac1, fac1, phi); */
/*     num_sin(fac1, fac1); */
/*     num_mul(fac1, fac1, r); */
/*     omega(fac2, r, phi, alpha, beta); */
/*     //num_set_d(fac2, omega(r, phi, alpha, beta)); */
/*     num_mul(fac2, fac2, z); */
/*     num_sub(n, fac1, fac2); */

/*     /\* denominator *\/ */
/*     num_pow_d(d, r, 2.0); */
/*     num_pow_d(fac1, z, 2.0); */
/*     num_set_d(fac2, -2.0*num_to_d(r)*cos(num_to_d(phi))); */
/*     num_mul(fac2, fac2, z); */
/*     num_add(d, d, fac1); */
/*     num_add(d, d, fac2); */
    
/*     A(res, r, alpha, beta, phi); */
/*     num_mul(res, res, n); */
/*     num_div(res, res, d); */
/*     num_mul_d(res, res, 1.0/M_PI); */

/*     delete(n), delete(d), delete(fac1), delete(fac2); */
/* } */

int
fn_C (acb_ptr res, const acb_t ph, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */
    integ_params_t* p = (integ_params_t*) param;

    arb_t pi, phi;
    acb_t w, aval, fac1, fac2, aux, aux2, I, rho;
    
    acb_init(w);
    arb_init(pi);
    acb_init(aval);
    acb_init(fac1);
    acb_init(fac2);
    acb_init(aux);
    acb_init(aux2);
    acb_init(I);
    arb_init(phi);
    acb_init(rho);

    arb_const_pi(pi, prec);
    acb_onei(I);
    acb_get_real(phi, ph);
    acb_set_arb(rho, p->c);

    fn_A (aval, rho, p->alpha, p->beta, phi);
    omega(w, p->c /* rho */, phi, p->alpha, p->beta);

    acb_mul(aux, w, I, prec);
    acb_exp(fac1, aux, prec);

    acb_mul_arb(aux, I, phi, prec);
    acb_exp(fac2, aux, prec);
    acb_mul(fac2, fac2, rho, prec);
    acb_sub(fac2, fac2, p->z, prec);

    acb_set_d(res, 2.0);
    acb_mul_arb(res, res, pi, prec);
    acb_inv(res, res, prec);
    acb_mul_arb(res, res, p->c /* rho */, prec);
    acb_mul(res, res, aval, prec);
    acb_mul(res, res, fac1, prec);
    acb_div(res, res, fac2, prec);

    acb_inv(res, res, prec);
    acb_mul(res, res, aval, prec);
    acb_mul(res, res, fac1, prec);
    acb_div(res, res, fac2, prec);
    
    acb_clear(w);
    arb_clear(pi);
    acb_clear(I);
    acb_clear(aux);
    acb_clear(aux2);
    acb_clear(fac1);
    acb_clear(fac2);
    arb_clear(phi);

    return 0;
}

/* void */
/* C (num_t res, */
/*    const num_t phi, */
/*    const num_t alpha, */
/*    const num_t beta, */
/*    const num_t z, */
/*    const num_t rho) */
/* { */
/*     num_t n, d, fac1, fac2, J; */

/*     n = new(num), d = new(num), J = new(num); */
/*     fac1 = new(num), fac2 = new(num); */

/*     num_set_d_d(J, 0.0, 1.0); */

/*     /\* numerator *\/ */
/*     omega(fac1, rho, phi, alpha, beta); */
/*     //num_set_d(fac1, omega(rho, phi, alpha, beta)); */
/*     num_sin(fac1, fac1); */
/*     num_mul(fac1, fac1, J); */
/*     omega(fac2, rho, phi, alpha, beta); */
/*     num_cos(fac2, fac2); */
/*     //num_set_d(fac2, cos(omega(rho, phi, alpha, beta))); */
/*     num_add(n, fac1, fac2); */

/*     /\* denominator *\/ */
/*     num_mul(fac1, J, phi); */
/*     //num_set_d_d(fac1, 0.0, phi); */
/*     num_exp(fac1, fac1); */
/*     num_mul(fac1, fac1, rho); */
/*     num_sub(d, fac1, z); */
        
/*     A(res, rho, alpha, beta, phi); */
/*     num_mul(res, res, n); */
/*     num_div(res, res, d); */
/*     num_mul_d(res, res, num_to_d(rho)/(2.0 * M_PI)); */

/*     delete(n), delete(d), delete(J), delete(fac1), delete(fac2); */
/* } */

void
quadb (acb_t res,
       const acb_t z,
       const arb_t phi,
       const acb_t from,
       const acb_t to,
       void * ctx)
{
    ctx_t* p = (ctx_t*) ctx;
    
    slong goal;
    mag_t tol;
    acb_calc_integrate_opt_t options;
    /* integ_params_t par = { */
    /*     .alpha = p->alpha, */
    /*     .beta = p->beta, */
    /*     .z = z, */
    /*     .c = phi */
    /* }; */
    integ_params_t par;

    arb_init(par.alpha);
    arb_init(par.beta);
    acb_init(par.z);
    arb_init(par.c);

    arb_set(par.alpha, p->alpha);
    arb_set(par.beta, p->beta);
    acb_set(par.z, z);
    arb_set(par.c, phi);

    goal = p->prec;
    mag_init(tol);

    acb_calc_integrate_opt_init(options);
    mag_set_ui_2exp_si(tol, 1, -p->prec);

    acb_calc_integrate(res, fn_B, &par, from, to, goal, tol, options, p->prec);

    arb_clear(par.alpha);
    arb_clear(par.beta);
    acb_clear(par.z);
    arb_clear(par.c);
}

/* void */
/* integrate_B (num_t res, */
/*              const num_t alpha, */
/*              const num_t beta, */
/*              const num_t z, */
/*              const num_t phi, */
/*              const num_t from, */
/*              const num_t to, */
/*              const num_t acc) */
/* { */
/*     UNUSED(acc); */
/*     parameters_B p = { .alpha = alpha, .beta = beta, .z = z, .phi = phi }; */

/*     num_function_t F; */
/*     F.function = &f_wrap_B; */
/*     F.params = &p; */

/*     quad(res, F, from, to, 0); */
/* } */

void
quadc (acb_t res,
       const acb_t z,
       const arb_t rho,
       const acb_t from,
       const acb_t to,
       void * ctx)
{
    ctx_t* p = (ctx_t*) ctx;
    slong goal;
    mag_t tol;
    acb_calc_integrate_opt_t options;
    /* integ_params_t par  = { */
    /*     .alpha = p->alpha, */
    /*     .beta = p->beta, */
    /*     .z = z, */
    /*     .c = rho */
    /* }; */
    integ_params_t par;

    arb_init(par.alpha);
    arb_init(par.beta);
    acb_init(par.z);
    arb_init(par.c);

    arb_set(par.alpha, p->alpha);
    arb_set(par.beta, p->beta);
    acb_set(par.z, z);
    arb_set(par.c, rho);

    goal = p->prec;
    mag_init(tol);

    acb_calc_integrate_opt_init(options);
    mag_set_ui_2exp_si(tol, 1, -p->prec);

    acb_calc_integrate(res, fn_C, &par, from, to, goal, tol, options, p->prec);

    arb_clear(par.alpha);
    arb_clear(par.beta);
    acb_clear(par.z);
    arb_clear(par.c);
}

/* void */
/* integrate_C (num_t res, */
/*              const num_t alpha, */
/*              const num_t beta, */
/*              const num_t z, */
/*              const num_t rho, */
/*              const num_t from, */
/*              const num_t to) */
/* { */
/*     parameters_C p = { .alpha = alpha, .beta = beta, .z = z, .rho = rho }; */
    
/*     num_function_t F; */
/*     F.function = &f_wrap_C; */
/*     F.params = &p; */
    
/*     quad(res, F, from, to, 0); */
/* } */

/* void */
/* f_wrap_B(num_t res, */
/*          const num_t z, */
/*          void * ctx) */
/* { */
/*     parameters_B* p = (parameters_B *) ctx; */
/*     B(res, z, p->alpha, p->beta, p->z, p->phi); */
/* } */

/* void */
/* f_wrap_C(num_t res, */
/*          const num_t z, */
/*          void * ctx) */
/* { */
/*     parameters_C* p = (parameters_C *) ctx; */
/*     C(res, z, p->alpha, p->beta, p->z, p->rho); */
/* } */
