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

//#include "quad.h"
//#include "new.h"
//#include "num.h"
//#include "utils.h"

#include <math.h>
#include <assert.h>

/* Parameters for integrating B */
typedef struct {
    arb_t alpha;
    arb_t beta;
    acb_t z;
    arb_t c;
} integ_params_t;

/* /\* Parameters for integrating B *\/ */
/* typedef struct { */
/*     num_t alpha; */
/*     num_t beta; */
/*     num_t z; */
/*     num_t phi; */
/* } parameters_B; */

/* /\* Parameters for integrating C *\/ */
/* typedef struct { */
/*     num_t alpha; */
/*     num_t beta; */
/*     num_t z; */
/*     num_t rho; */
/* } parameters_C; */

/* void */
/* f_wrap_B(num_t res, */
/*          const num_t z, */
/*          void * ctx); */

/* void */
/* f_wrap_C(num_t res, */
/*          const num_t z, */
/*          void * ctx); */

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

    /* acb_t zp, cosxa, expz; */
    /* acb_init(zp); */
    /* acb_init(cosxa)= new(num), cosxa = new(num), expz = new(num); */

    /* /\* compute z**((1-beta)/alpha) *\/ */
    /* num_set_num(zp, z); */
    /* num_pow_d(zp, zp, (1.0 - num_to_d(beta))/num_to_d(alpha)); */
    /* /\* compute cos(x/alpha) *\/ */
    /* num_div(cosxa, x, alpha); */
    /* //num_set_d(cosxa, cosxa); */
    /* num_cos(cosxa, cosxa); */
    /* /\* compute exp(z**(1/alpha) * cos(x/alpha)) *\/ */
    /* //num_set_d_d(expz, rez, imz); */
    /* num_pow_d(expz, z, 1.0/num_to_d(alpha)); */
    /* num_mul(expz, expz, cosxa); */
    /* num_exp(expz, expz); */

    /* num_mul(res, zp, expz); */
    /* num_mul_d(res, res, 1.0/num_to_d(alpha)); */
       
    /* delete(zp), delete(cosxa), delete(expz); */
}

static void
omega(arb_t res,
      const arb_t x,
      const arb_t y,
      const arb_t alpha,
      const arb_t beta)
{
    /* x**(1.0/alpha) * sin(y/alpha) + y * (1.0 + (1.0 - beta)/alpha) */
    arb_t fac1, fac2, aux1, aux2;

    arb_init(fac1);
    arb_init(fac2);
    arb_init(aux1);
    arb_init(aux2);
    
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

    arb_add(res, fac1, fac2, PREC);

    arb_clear(fac1);
    arb_clear(fac2);
    arb_clear(aux1);
    arb_clear(aux2);
}

int
fn_B (acb_ptr res, const acb_t r, void * param, slong order, slong prec)
{
    integ_params_t* p = (integ_params_t*) param;

    arb_t w, pi;
    acb_t aval, fac1, fac2, aux, aux2;

    arb_init(w);
    arb_init(pi);
    acb_init(aval);
    acb_init(fac1);
    acb_init(fac2);
    acb_init(aux);
    acb_init(aux2);

    fn_A (aval, r, p->alpha, p->beta, p->c /* phi */);
    omega(w, r, p->c /* phi */, p->alpha, p->beta);

    acb_sub(aux, w, p->c, PREC);
    acb_sin(aux, aux, PREC);
    acb_mul(fac1, r, aux, PREC);
    acb_sin(aux, w, PREC);
    acb_mul(aux, aux, p->z, PREC);
    acb_sub(fac1, fac1, aux, PREC);

    acb_sqr(fac2, r, PREC);
    acb_sqr(aux, p->z, PREC);
    acb_add(fac2, fac2, aux, PREC);
    acb_cos(aux2, p->c, PREC);
    acb_mul(aux2, aux2, p->z, PREC);
    acb_mul(aux2, aux2, r, PREC);
    acb_set_d(aux, -2.0);
    acb_mul(aux, aux, aux2, PREC);
    acb_add(fac2, fac2, aux, PREC);

    arb_const_pi(pi, PREC);

    acb_set_arb(res, pi);
    acb_inv(res, res, PREC);
    acb_mul(res, res, aval, PREC);
    acb_mul(res, res, fac1, PREC);
    acb_div(res, res, fac2, PREC);
    
    arb_clear(w);
    arb_clear(pi);
    acb_clear(aux);
    acb_clear(aux2);
    acb_clear(fac1);
    acb_clear(fac2);
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
    integ_params_t* p = (integ_params_t*) param;

    arb_t w, pi, phi;
    acb_t aval, fac1, fac2, aux, aux2, I;
    
    arb_init(w);
    arb_init(pi);
    acb_init(aval);
    acb_init(fac1);
    acb_init(fac2);
    acb_init(aux);
    acb_init(aux2);
    acb_init(I);
    arb_init(phi);

    arb_const_pi(pi, PREC);
    acb_onei(I);
    acb_get_real(phi, ph);

    fn_A (aval, p->c /* rho */, p->alpha, p->beta, phi);
    omega(w, p->c /* rho */, phi, p->alpha, p->beta);

    acb_mul(aux, w, I, PREC);
    acb_exp(fac1, aux, PREC);

    acb_mul(aux, phi, I, PREC);
    acb_exp(fac2, aux, PREC);
    acb_mul(fac2, fac2, p->c /* rho */, PREC);
    acb_sub(fac2, fac2, p->z, PREC);

    acb_set_d(res, 2.0);
    acb_mul(res, pi, res, PREC);
    acb_inv(res, res, PREC);
    acb_mul(res, res, p->c /* rho */, PREC);
    acb_mul(res, res, aval, PREC);
    acb_mul(res, res, fac1, PREC);
    acb_div(res, res, fac2, PREC);

    acb_inv(res, res, PREC);
    acb_mul(res, res, aval, PREC);
    acb_mul(res, res, fac1, PREC);
    acb_div(res, res, fac2, PREC);
    
    arb_clear(w);
    arb_clear(pi);
    acb_clear(I);
    acb_clear(aux);
    acb_clear(aux2);
    acb_clear(fac1);
    acb_clear(fac2);
    arb_clear(phi);
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
       const arb_t alpha,
       const arb_t beta,
       const acb_t z,
       const arb_t phi,
       const acb_t from,
       const acb_t to)
{
    slong prec, goal;
    mag_t tol;
    acb_calc_integrate_opt_t options;
    integ_params_t p  = {
        .alpha = alpha,
        .beta = beta,
        .z = z,
        .c = phi
    };

    prec = PREC;
    goal = prec;
    mag_init(tol);

    acb_calc_integrate_opt_init(options);
    mag_set_ui_2exp_si(tol, 1, -prec);

    acb_calc_integrate(res, fn_B, &p, from, to, goal, tol, options, prec);
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
       const arb_t alpha,
       const arb_t beta,
       const acb_t z,
       const arb_t rho,
       const acb_t from,
       const acb_t to)
{
    slong prec, goal;
    mag_t tol;
    acb_calc_integrate_opt_t options;
    integ_params_t p  = {
        .alpha = alpha,
        .beta = beta,
        .z = z,
        .c = rho
    };

    prec = PREC;
    goal = prec;
    mag_init(tol);

    acb_calc_integrate_opt_init(options);
    mag_set_ui_2exp_si(tol, 1, -prec);

    acb_calc_integrate(res, fn_C, &p, from, to, goal, tol, options, prec);
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
