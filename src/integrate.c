#include <assert.h>

#include "integrate.h"

#include "log.h"
#include "new.h"
#include "num.h"

#include "flint/profiler.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_dirichlet.h"
#include "acb_modular.h"
#include "acb_calc.h"

/* Parameters for integrating B */
typedef struct {
    double alpha;
    double beta;
    double x;
    double y;
    double phi;
} parameters_B;

/* Parameters for integrating C */
typedef struct {
    double alpha;
    double beta;
    double x;
    double y;
    double rho;
} parameters_C;

// Converts an arb_t number to double.
static double
arbtod (const arb_t x)
{
    return arf_get_d(arb_midref(x), ARF_RND_NEAR);
}

void
num_set_acb (num_t res, const acb_t x)
{
    arb_t re, im;    
    arb_init(re); arb_init(im);
    acb_get_real(re, x);
    acb_get_imag(im, x);
    num_set_d_d(res, arbtod(re), arbtod(im));
    arb_clear(re), arb_clear(im);
}

static int
f_wrap_B(acb_ptr res, const acb_t z, void * param, slong order, slong prec);

static int
f_wrap_C(acb_ptr res, const acb_t z, void * params, slong order, slong prec);

/* num_t */
/* A (const num_t z, */
/*    const num_t alpha, */
/*    const num_t beta, */
/*    const num_t x) */
/* { */
/*     const double _alpha = num_to_d(alpha); */
/*     const double _beta = num_to_d(beta); */
    
/*     return num_mul( */
/*         num_inv(alpha), */
/*         num_mul( */
/*             num_pow(z, new(num, (1 - _beta)/_alpha, 0.0)), */
/*             num_exp(num_mul(num_pow(z, num_inv(alpha)), num_cos(num_div(x, alpha)))))); */
/* } */

void
A (num_t res,
   const double rez,
   const double imz,
   const double alpha,
   const double beta,
   const double x)
{
    num_t zp, cosxa, expz;
    zp = new(num), cosxa = new(num), expz = new(num);

    /* compute z**((1-beta)/alpha) */
    num_set_d_d(zp, rez, imz);
    num_pow_d(zp, zp, (1.0 - beta)/alpha);
    /* compute cos(x/alpha) */
    num_set_d(cosxa, x/alpha);
    num_cos(cosxa, cosxa);
    /* compute exp(z**(1/alpha) * cos(x/alpha)) */
    num_set_d_d(expz, rez, imz);
    num_pow_d(expz, expz, 1.0/alpha);
    num_mul(expz, expz, cosxa);
    num_exp(expz, expz);

    num_mul(res, zp, expz);
    num_mul_d(res, res, 1.0/alpha);
       
    delete(zp), delete(cosxa), delete(expz);
}

double
omega(const double x, const double y, const double alpha, const double beta)
{
    return pow(x, 1/alpha)*sin(y/alpha) + y*(1 + (1 - beta)/alpha);
}

/* num_t */
/* B (const num_t r, */
/*    const num_t alpha, */
/*    const num_t beta, */
/*    const num_t z, */
/*    const num_t phi) */
/* { */
/*     num_t two = new(num, 2.0, 0.0); */
/*     num_t one_over_pi = new(num, 1.0/M_PI, 0.0); */
/*     num_t num, den; */

/*     num = num_sub(num_mul(r, num_sin(num_sub(omega(r, phi, alpha, beta), phi))), num_mul(z, num_sin(omega(r, phi, alpha, beta)))); */
/*     den = num_sub( */
/*         num_add(num_pow(r, two), num_pow(z, two)), */
/*         num_mul(two, num_mul(r, num_mul(z, num_cos(phi))))); */

/*     num_t res =  num_mul( */
/*         one_over_pi, */
/*         num_mul(A(r, alpha, beta, phi), num_div(num, den))); */

/*     delete(two); delete(num); delete(den); delete(one_over_pi); */
    
/*     return res; */
/* } */

void
B (num_t res,
   const double r,
   const double alpha,
   const double beta,
   const double x,
   const double y,
   const double phi)
{
    num_t z, n, d, fac1, fac2;

    z = new(num), n = new(num), d = new(num);
    fac1 = new(num), fac2 = new(num);

    num_set_d_d(z, x, y);

    /* numerator */
    num_set_d(fac1, omega(r, phi, alpha, beta) - phi);
    num_sin(fac1, fac1);
    num_mul_d(fac1, fac1, r);
    num_set_d(fac2, omega(r, phi, alpha, beta));
    num_mul(fac2, fac2, z);
    num_sub(n, fac1, fac2);

    /* denominator */
    num_set_d(d, r*r);
    num_pow_d(fac1, z, 2.0);
    num_set_d(fac2, -2.0*r*cos(phi));
    num_mul(fac2, fac2, z);
    num_add(d, d, fac1);
    num_add(d, d, fac2);
    
    A(res, r, 0.0, alpha, beta, phi);
    num_mul(res, res, n);
    num_div(res, res, d);
    num_mul_d(res, res, 1.0/M_PI);

    delete(z), delete(n), delete(d), delete(fac1), delete(fac2);
}

void
C (num_t res,
   const double phi,
   const double alpha,
   const double beta,
   const double x,
   const double y,
   const double rho)
{
    num_t z, n, d, fac1, fac2, J;

    z = new(num), n = new(num), d = new(num), J = new(num);
    fac1 = new(num), fac2 = new(num);

    num_onei(J);
    num_set_d_d(z, x, y);

    /* numerator */
    num_set_d(fac1, omega(rho, phi, alpha, beta));
    num_sin(fac1, fac1);
    num_mul(fac1, fac1, J);
    num_set_d(fac2, cos(omega(rho, phi, alpha, beta)));
    num_add(n, fac1, fac2);

    /* denominator */
    num_set_d_d(fac1, 0.0, phi);
    num_exp(fac1, fac1);
    num_mul_d(fac1, fac1, rho);
    num_sub(d, fac1, z);
        
    A(res, rho, 0.0, alpha, beta, phi);
    num_mul(res, res, n);
    num_div(res, res, d);
    num_mul_d(res, res, rho/(2.0*M_PI));

    delete(z), delete(n), delete(d), delete(J), delete(fac1), delete(fac2);
}

void
integrate_B (num_t res,
             const double alpha,
             const double beta,
             const double x,
             const double y,
             const double phi,
             const double from,
             const double to)
{
    acb_t _res, t, _from, _to;
    mag_t tol;
    slong prec, goal;
    acb_calc_integrate_opt_t options;

    acb_calc_integrate_opt_init(options);

    prec = 53;
    goal = 20;
    
    acb_init(_from);
    acb_init(_to);
    acb_init(_res);
    acb_init(t);
    mag_init(tol);

    parameters_B p = { .alpha = alpha, .beta = beta, .x = x, .y = y, .phi = phi };

    acb_set_d(_from, from);
    acb_set_d(  _to, to);
    int status = acb_calc_integrate(_res, f_wrap_B, &p, _from, _to, goal, tol, options, prec);
    num_set_acb(res, _res);
    acb_clear(_res);
    acb_clear(t);
    acb_clear(_from);
    acb_clear(_to);
    mag_clear(tol);

    flint_cleanup_master();
}

void
integrate_C (num_t res,
             const double alpha,
             const double beta,
             const double x,
             const double y,
             const double rho,
             const double from,
             const double to)
{
    acb_t _res, t, _from, _to;
    mag_t tol;
    slong prec, goal;
    acb_calc_integrate_opt_t options;

    acb_calc_integrate_opt_init(options);

    prec = 53;
    goal = 20;
    
    acb_init(_from);
    acb_init(_to);
    acb_init(_res);
    acb_init(t);
    mag_init(tol);

    parameters_C p = { .alpha = alpha, .beta = beta, .x = x, .y = y, .rho = rho };

    acb_set_d(_from, from);
    acb_set_d(  _to, to);
    int status = acb_calc_integrate(_res, f_wrap_C, &p, _from, _to, goal, tol, options, prec);

    num_set_acb(res, _res);
    acb_clear(_res);
    acb_clear(t);
    acb_clear(_from);
    acb_clear(_to);
    mag_clear(tol);

    flint_cleanup_master();
}

static int
f_wrap_B(acb_ptr res, const acb_t z, void * params, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    parameters_B* p = (parameters_B *) params;
    
    num_t _z = new(num);
    num_set_acb(_z, z);
    assert(num_is_real(_z));
    const double r = num_to_d(_z);
    delete(_z);
    
    num_t _res;
    _res = new(num);
    B(_res, r, p->alpha, p->beta, p->x, p->y, p->phi);
    const double complex __res = num_to_complex(_res);
    delete(_res);
    
    acb_set_d_d(res, creal(__res), cimag(__res));

    return 0;
}

static int
f_wrap_C(acb_ptr res, const acb_t z, void * params, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    parameters_C* p = (parameters_C *) params;

    num_t _z = new(num);
    num_set_acb(_z, z);
    assert(num_is_real(_z));
    const double r = num_to_d(_z);
    delete(_z);

    num_t _res;
    _res = new(num);
    C(_res, r, p->alpha, p->beta, p->x, p->y, p->rho);
    const double complex __res = num_to_complex(_res);
    delete(_res);
    
    acb_set_d_d(res, creal(__res), cimag(__res));
        
    return 0;
}
