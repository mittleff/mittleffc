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
    num_t alpha;
    num_t beta;
    num_t z;
    num_t phi;
} parameters_B;

/* Parameters for integrating C */
typedef struct {
    num_t alpha;
    num_t beta;
    num_t z;
    num_t rho;
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
   const num_t z,
   const num_t alpha,
   const num_t beta,
   const num_t x)
{
    num_t zp, cosxa, expz;
    zp = new(num), cosxa = new(num), expz = new(num);

    /* compute z**((1-beta)/alpha) */
    num_set(zp, z);
    num_pow_d(zp, zp, (1.0 - num_to_d(beta))/num_to_d(alpha));
    /* compute cos(x/alpha) */
    num_div(cosxa, x, alpha);
    //num_set_d(cosxa, cosxa);
    num_cos(cosxa, cosxa);
    /* compute exp(z**(1/alpha) * cos(x/alpha)) */
    //num_set_d_d(expz, rez, imz);
    num_pow_d(expz, z, 1.0/num_to_d(alpha));
    num_mul(expz, expz, cosxa);
    num_exp(expz, expz);

    num_mul(res, zp, expz);
    num_mul_d(res, res, 1.0/num_to_d(alpha));
       
    delete(zp), delete(cosxa), delete(expz);
}

static void
omega(num_t res,
      const num_t x,
      const num_t y,
      const num_t alpha,
      const num_t beta)
{
    num_t one_over_alpha, aux, fac1, fac2;

    one_over_alpha = new(num);
    aux = new(num);
    fac1 = new(num), fac2 = new(num);
    
    num_inv(one_over_alpha, alpha);

    num_div(aux, y, alpha);
    num_sin(aux, aux);
    num_pow(fac1, x, one_over_alpha);
    num_mul(fac1, fac1, aux);

    num_set_d(aux, 1.0);
    num_sub(aux, aux, beta);
    num_add(aux, aux, alpha);
    num_div(aux, aux, alpha);
    num_mul(fac2, aux, y);
    
    //return pow(x, 1/alpha)*sin(y/alpha) + y*(1 + (1 - beta)/alpha);

    num_add(res, fac1, fac2);
    delete(one_over_alpha);
    delete(aux), delete(fac1), delete(fac2);
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
   const num_t r,
   const num_t alpha,
   const num_t beta,
   const num_t z,
   const num_t phi)
{
    num_t n, d, fac1, fac2;

    n = new(num), d = new(num);
    fac1 = new(num), fac2 = new(num);

    /* numerator */
    omega(fac1, r, phi, alpha, beta);
    num_sub(fac1, fac1, phi);
    num_sin(fac1, fac1);
    num_mul(fac1, fac1, r);
    omega(fac2, r, phi, alpha, beta);
    //num_set_d(fac2, omega(r, phi, alpha, beta));
    num_mul(fac2, fac2, z);
    num_sub(n, fac1, fac2);

    /* denominator */
    num_pow_d(d, r, 2.0);
    num_pow_d(fac1, z, 2.0);
    num_set_d(fac2, -2.0*num_to_d(r)*cos(num_to_d(phi)));
    num_mul(fac2, fac2, z);
    num_add(d, d, fac1);
    num_add(d, d, fac2);
    
    A(res, r, alpha, beta, phi);
    num_mul(res, res, n);
    num_div(res, res, d);
    num_mul_d(res, res, 1.0/M_PI);

    delete(n), delete(d), delete(fac1), delete(fac2);
}

void
C (num_t res,
   const num_t phi,
   const num_t alpha,
   const num_t beta,
   const num_t z,
   const num_t rho)
{
    num_t n, d, fac1, fac2, J;

    n = new(num), d = new(num), J = new(num);
    fac1 = new(num), fac2 = new(num);

    num_onei(J);

    /* numerator */
    omega(fac1, rho, phi, alpha, beta);
    //num_set_d(fac1, omega(rho, phi, alpha, beta));
    num_sin(fac1, fac1);
    num_mul(fac1, fac1, J);
    omega(fac2, rho, phi, alpha, beta);
    num_cos(fac2, fac2);
    //num_set_d(fac2, cos(omega(rho, phi, alpha, beta)));
    num_add(n, fac1, fac2);

    /* denominator */
    num_mul(fac1, J, phi);
    //num_set_d_d(fac1, 0.0, phi);
    num_exp(fac1, fac1);
    num_mul(fac1, fac1, rho);
    num_sub(d, fac1, z);
        
    A(res, rho, alpha, beta, phi);
    num_mul(res, res, n);
    num_div(res, res, d);
    num_mul_d(res, res, num_to_d(rho)/(2.0 * M_PI));

    delete(n), delete(d), delete(J), delete(fac1), delete(fac2);
}

void
integrate_B (num_t res,
             const num_t alpha,
             const num_t beta,
             const num_t z,
             const num_t phi,
             const num_t from,
             const num_t to,
             const num_t acc)
{
    log_trace("[%s] alpha=%g, beta=%g, z=%g+%g acc=%g",
              __func__,
              num_to_d(alpha),
              num_to_d(beta),
              num_to_complex(z),
              num_to_d(acc));
    acb_t _res, t, _from, _to;
    mag_t tol;
    slong prec, goal;
    acb_calc_integrate_opt_t options;

    acb_calc_integrate_opt_init(options);

    prec = 63;
    goal = prec;
    
    acb_init(_from);
    acb_init(_to);
    acb_init(_res);
    acb_init(t);
    mag_init(tol);

    mag_set_d(tol, num_to_d(acc));//1e-15);
    parameters_B p = { .alpha = alpha, .beta = beta, .z = z, .phi = phi };

    acb_set_d(_from, num_to_d(from));
    acb_set_d(  _to, num_to_d(to));
    int status = acb_calc_integrate(_res, &f_wrap_B, &p, _from, _to, goal, tol, options, prec);
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
             const num_t alpha,
             const num_t beta,
             const num_t z,
             const num_t rho,
             const num_t from,
             const num_t to)
{
    acb_t _res, t, _from, _to;
    mag_t tol;
    slong prec, goal;
    acb_calc_integrate_opt_t options;

    acb_calc_integrate_opt_init(options);

    prec = 53;
    goal = prec;
    
    acb_init(_from);
    acb_init(_to);
    acb_init(_res);
    acb_init(t);
    mag_init(tol);

    mag_set_d(tol, 1e-16);
    parameters_C p = { .alpha = alpha, .beta = beta, .z = z, .rho = rho };

    acb_set_d(_from, num_to_d(from));
    acb_set_d(  _to, num_to_d(to));
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
    //const double r = num_to_d(_z);
    
    
    num_t _res;
    _res = new(num);
    B(_res, _z, p->alpha, p->beta, p->z, p->phi);
    const double complex __res = num_to_complex(_res);
    printf("%g+%g\n", __res);
    delete(_res);
    delete(_z);
    
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
    //const double r = num_to_d(_z);
    

    num_t _res;
    _res = new(num);
    C(_res, _z, p->alpha, p->beta, p->z, p->rho);
    const double complex __res = num_to_complex(_res);
    delete(_res);
    delete(_z);
    
    acb_set_d_d(res, creal(__res), cimag(__res));
        
    return 0;
}
