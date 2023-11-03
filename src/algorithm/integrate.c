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
    double a;
} parameters_B;

// Converts an arb_t number to double.
static double
arbtod (const arb_t x)
{
    return arf_get_d(arb_midref(x), ARF_RND_NEAR);
}

static num_t
num_from_acb (const acb_t x)
{
    arb_t re, im;
    
    arb_init(re); arb_init(im);
    acb_get_real(re, x);
    acb_get_imag(im, x);
    return new(num, arbtod(re), arbtod(im));
}



static int
f_sin(acb_ptr res, const acb_t z, void * param, slong order, slong prec);

num_t
A (const num_t z,
   const num_t alpha,
   const num_t beta,
   const num_t x)
{
    const double _alpha = num_to_double(alpha);
    const double _beta = num_to_double(beta);
    
    return num_mul(
        num_inverse(alpha),
        num_mul(
            num_pow(z, new(num, (1 - _beta)/_alpha, 0.0)),
            num_exp(num_mul(num_pow(z, num_inverse(alpha)), num_cos(num_div(x, alpha))))));
}

num_t
omega(const num_t x, const num_t y, const num_t alpha, const num_t beta)
{
    num_t one = new(num, 1.0, 0.0);
    num_t two = new(num, 2.0, 0.0);
    return num_add(
        num_mul(num_pow(x, num_inverse(alpha)), num_sin(num_div(y, alpha))),
        num_mul(y, num_add(one, num_div(num_sub(one, beta), alpha))));
}

num_t
B (const num_t r,
   const num_t alpha,
   const num_t beta,
   const num_t z,
   const num_t phi)
{
    num_t two = new(num, 2.0, 0.0);
    num_t num, den;

    num = num_sub(num_mul(r, num_sin(num_sub(omega(r, phi, alpha, beta), phi))), num_mul(z, num_sin(omega(r, phi, alpha, beta))));
    den = num_sub(
        num_add(num_pow(r, two), num_pow(z, two)),
        num_mul(two, num_mul(r, num_mul(z, num_cos(phi)))));
    
    return num_mul(
        new(num, 1.0/M_PI, 0.0),
        num_mul(A(r, alpha, beta, phi), num_div(num, den)));
}


num_t
integrate_B (const num_t alpha,
             const num_t beta,
             const num_t z,
             const num_t phi,
             const num_t from,
             const num_t to)
{
    acb_t _res, t, _from, _to;
    mag_t tol;
    slong prec, goal;
    acb_calc_integrate_opt_t options;

    acb_calc_integrate_opt_init(options);

    prec = 53;
    goal = 10;
    
    acb_init(_from);
    acb_init(_to);
    acb_init(_res);
    acb_init(t);
    mag_init(tol);

    parameters_B p = { .a = 3.0 };

    acb_set_d(_from, num_to_double(from));
    acb_set_d(  _to, num_to_double(to));
    int status = acb_calc_integrate(_res, f_sin, &p, _from, _to, goal, tol, options, prec);

    num_t res = num_from_acb(_res);
    log_trace("===================>%d %g %g", status, num_real_d(res), num_imag_d(res));

    acb_clear(_res);
    acb_clear(t);
    acb_clear(_from);
    acb_clear(_to);
    mag_clear(tol);

    flint_cleanup_master();
                
    return res;
}

/* f(z) = sin(z) */
static int
f_sin(acb_ptr res, const acb_t z, void * params, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    parameters_B p = *(parameters_B *) params;

    acb_t a, x2, y;
    acb_init(a); acb_init(x2); acb_init(y);
    acb_set_d(a, p.a);
    acb_set_d(y, 2);
    acb_pow(x2, z, y, 53);
    
    acb_mul(res, a, x2, 53); // INtegration should return 1

    return 0;
}
