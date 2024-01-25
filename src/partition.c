#include "partition.h"
#include "log.h"
#include "num.h"
//#include "new.h"

#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_math.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

// #define R0 0.95
// #define GSL_EPSILON_FCMP 1e-15

static bool diskp (const num_t z, const num_t r);
static bool closure_diskp (const num_t z, const num_t r0);
static bool wedgep (const num_t z, const num_t phi1, const num_t phi2);
static bool closure_wedgep (const num_t z, const num_t phi1, const num_t phi2);
static bool is_between (const num_t c, const num_t a, const num_t b, const bool eq);

/* Equation (4.21) */
//static double
static void
compute_r1 (num_t r1, const num_t alpha, const num_t acc)
{
    num_t one, pi, c0, n, d, aux;

    pi = num_init(), one = num_init();
    num_set_d(pi, M_PI);
    num_set_d(one, 1.0);

    /* Compute C0 */
    /* pow(1.3, 1.0 - alpha) */
    n = num_init(), aux = num_init();
    num_set_d(n, 1.3);
    num_sub(aux, one, alpha);
    num_pow(n, n, aux);
    /* M_PI * sin(M_PI * alpha) */
    d = num_init();
    num_mul(d, pi, alpha);
    num_sin(d, d);
    num_mul(d, d, pi);
    c0 = num_init();
    num_div(c0, n, d);

    num_div(aux, acc, c0);
    num_log(aux, aux);
    num_mul_d(aux, aux, -2.0);
    num_pow(r1, aux, alpha);
    
    num_clear(pi), num_clear(one), num_clear(c0), num_clear(n), num_clear(d), num_clear(aux);    
}

bool
in_region_G0 (const num_t z)
{
    bool ret;
    num_t r0;

    r0 = num_init();
    num_set_d(r0, 0.95);
    ret = closure_diskp(z, r0);
    num_clear(r0);

    return ret;
}

bool
in_region_G1 (const num_t z, const num_t alpha, const num_t acc)
{
    bool ret;
    num_t r1, delta, phi1, phi2;

    r1 = num_init();
    compute_r1(r1, alpha, acc);
    
    /* Compute delta */
    delta = num_init();
    num_set_d(delta, M_PI * num_to_d(alpha)/8.0);
    
    /* Compute phi1, phi2 */
    phi1 = num_init(), phi2 = num_init();
    num_set_d(phi1, - M_PI * num_to_d(alpha) + num_to_d(delta) );
    num_set_d(phi2, + M_PI * num_to_d(alpha) - num_to_d(delta));

    ret = (!diskp(z, r1)) && wedgep(z, phi1, phi2);
    
    num_clear(r1), num_clear(delta), num_clear(phi1), num_clear(phi2);

    return ret;
}

bool
in_region_G2 (const num_t z, const num_t alpha, const num_t acc)
{
    bool ret;
    num_t r1, deltat, phi1, phi2;

    r1 = num_init();
    compute_r1(r1, alpha, acc);

    /* Compute deltat */
    deltat = num_init();
    num_set_d(deltat, MIN(M_PI * num_to_d(alpha)/8.0, M_PI * (num_to_d(alpha) + 1.0)/2.0));

    /* Compute phi1, phi2 */
    phi1 = num_init(), phi2 = num_init();
    num_set_d(phi1, M_PI * num_to_d(alpha) + num_to_d(deltat));
    num_set_d(phi2, -num_to_d(phi1));

    ret = (!diskp(z, r1)) && wedgep(z, phi1, phi2);

    num_clear(r1), num_clear(deltat), num_clear(phi1), num_clear(phi2);

    return ret;
}

bool
in_region_G3 (const num_t z, const num_t alpha, const num_t acc)
{
    bool ret;
    num_t r1, delta, deltat, phi1, phi2;

    r1 = num_init();
    compute_r1(r1, alpha, acc);

    /* Compute delta, deltat */
    delta = num_init(), deltat = num_init();
    num_set_d(delta, M_PI * num_to_d(alpha)/8.0);
    num_set_d(deltat, MIN(M_PI * num_to_d(alpha)/8.0, M_PI * (num_to_d(alpha) + 1.0)/2.0));

    /* Compute phi1, phi2 */
    phi1 = num_init(), phi2 = num_init();
    num_set_d(phi1, M_PI * num_to_d(alpha) - num_to_d(delta));
    num_set_d(phi2, M_PI * num_to_d(alpha) + num_to_d(deltat));

    ret = (!diskp(z, r1)) && closure_wedgep(z, phi1, phi2);

    num_clear(r1), num_clear(deltat), num_clear(phi1), num_clear(phi2);

    return ret;
}

bool
in_region_G4 (const num_t z, const num_t alpha, const num_t acc)
{
    bool ret;
    num_t r1, delta, deltat, phi1, phi2;

    r1 = num_init();
    compute_r1(r1, alpha, acc);

    /* Compute delta, deltat */
    delta = num_init(), deltat = num_init();
    num_set_d(delta, M_PI * num_to_d(alpha)/8.0);
    num_set_d(deltat, MIN(M_PI * num_to_d(alpha)/8.0, M_PI * (num_to_d(alpha) + 1.0)/2.0));

    /* Compute phi1, phi2 */
    phi1 = num_init(), phi2 = num_init();
    num_set_d(phi1, -(M_PI * num_to_d(alpha) + num_to_d(deltat)));
    num_set_d(phi2, -M_PI * num_to_d(alpha) + num_to_d(delta));

    ret = (!diskp(z, r1)) && closure_wedgep(z, phi1, phi2);

    num_clear(r1), num_clear(delta), num_clear(deltat), num_clear(phi1), num_clear(phi2);

    return ret;
}

bool
in_region_G5 (const num_t z, const num_t alpha, const num_t acc)
{
    bool ret;
    num_t r0, r1, phi1, phi2;

    r0 = num_init();
    num_set_d(r0, 0.95);

    r1 = num_init();
    compute_r1(r1, alpha, acc);
   
    /* Compute phi1, phi2 */
    phi1 = num_init(), phi2 = num_init();
    num_set_d(phi1, -5.0 * M_PI * num_to_d(alpha)/6.0);
    num_set_d(phi2, +5.0 * M_PI * num_to_d(alpha)/6.0);

    ret = diskp(z, r1) && (closure_wedgep(z, phi1, phi2) && !diskp(z, r0));

    num_clear(r0), num_clear(r1), num_clear(phi1), num_clear(phi2);

    return ret;    
}

bool
in_region_G6 (const num_t z, const num_t alpha, const num_t acc)
{
    bool ret;
    num_t r0, r1, phi1, phi2;

    r0 = num_init();
    num_set_d(r0, 0.95);

    r1 = num_init();
    compute_r1(r1, alpha, acc);
   
    /* Compute phi1, phi2 */
    phi1 = num_init(), phi2 = num_init();
    num_set_d(phi1, +5.0 * M_PI * num_to_d(alpha)/6.0);
    num_set_d(phi2, -5.0 * M_PI * num_to_d(alpha)/6.0);

    ret = diskp(z, r1) && (closure_wedgep(z, phi1, phi2) && !diskp(z, r0));

    num_clear(r0), num_clear(r1), num_clear(phi1), num_clear(phi2);

    return ret;
}


static bool
_diskp (const num_t z, const num_t r, const bool closure)
{
    num_t absz;
    bool ret;
    
    absz = num_init();
    num_abs(absz, z);
    if (closure == true)
        ret = diskp(z, r) || num_eq(absz, r);
    else
        ret = num_lt(absz, r);
    num_clear(absz);
    
    return ret;
}

static bool
diskp (const num_t z, const num_t r)
{
    return _diskp(z, r, false);
}

static bool
closure_diskp (const num_t z, const num_t r)
{
    return _diskp(z, r, true);
}

static bool
_wedgep (const num_t z, const num_t phi1, const num_t phi2, const double closure)
{
    bool res;
    num_t argz;

    argz = num_init();
    num_arg(argz, z);
    res = is_between(argz, phi1, phi2, closure);
    num_clear(argz);
    return res;
}

static bool
wedgep (const num_t z, const num_t phi1, const num_t phi2)
{
    return _wedgep(z, phi1, phi2, false);
}

static bool
closure_wedgep (const num_t z, const num_t phi1, const num_t phi2)
{
    return _wedgep(z, phi1, phi2, true);
}

/* Checks if _c is between _a and _b */
static bool
is_between (const num_t c, const num_t a, const num_t b, const bool eq)
{
    bool ret;
    num_t _a, _b, _c, n;
    
    _a = num_init(), _b = num_init(), _c = num_init(), n = num_init();
    num_set_d(n, 2.0 * M_PI);

    num_fmod(_a, a, n);
    num_fmod(_b, b, n);
    num_fmod(_c, c, n);

    if (num_lt(_a, _b))
    {
        if (eq)
            ret = (num_le(_a, _c) && num_le(_c, _b)) ? true : false;
        else
            ret = (num_lt(_a, _c) && num_lt(_c, _b)) ? true : false;
    }
    else
    { /* b < a */
        /* if in [b, a] then not in [a, b] */
        if (eq)
            ret = (num_le(_b, _c) && num_le(_c, _a)) ? false : true;
        else
            ret = (num_lt(_b, _c) && num_lt(_c, _a)) ? false : true;
    }

    num_clear(_a), num_clear(_b), num_clear(_c), num_clear(n);

    return ret;
}


