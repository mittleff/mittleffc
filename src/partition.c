#include "partition.h"
#include "log.h"
#include "num.h"
#include "new.h"

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
static double
compute_r1 (num_t r1, const num_t alpha, const num_t acc)
{
    num_t one, pi, c0, n, d, aux;

    pi = new(num), one = new(num);
    num_set_d(pi, M_PI);
    num_one(one);

    /* Compute C0 */
    /* pow(1.3, 1.0 - alpha) */
    n = new(num), aux = new(num);
    num_set_d(n, 1.3);
    num_sub(aux, one, alpha);
    num_pow(n, n, aux);
    /* M_PI * sin(M_PI * alpha) */
    d = new(num);
    num_mul(d, pi, alpha);
    num_sin(d, d);
    num_mul(d, d, pi);
    c0 = new(num);
    num_div(c0, n, d);

    num_div(aux, acc, c0);
    num_log(aux, aux);
    num_mul_d(aux, aux, -2.0);
    num_pow(r1, aux, alpha);
    
    delete(pi), delete(one), delete(c0), delete(n), delete(d), delete(aux);    
}

bool
in_region_G0 (const num_t z)
{
    bool ret;
    num_t r0;

    r0 = new(num);
    num_set_d(r0, 0.95);
    ret = closure_diskp(z, r0);
    delete(r0);

    return ret;
}

bool
in_region_G1 (const num_t z, const num_t alpha, const num_t acc)
{
    bool ret;
    num_t r1, delta, phi1, phi2;

    r1 = new(num);
    compute_r1(r1, alpha, acc);
    
    /* Compute delta */
    delta = new(num);
    num_set_d(delta, M_PI * num_to_d(alpha)/8.0);
    
    /* Compute phi1, phi2 */
    phi1 = new(num), phi2 = new(num);
    num_set_d(phi1, - M_PI * num_to_d(alpha) + num_to_d(delta) );
    num_set_d(phi2, + M_PI * num_to_d(alpha) - num_to_d(delta));

    ret = (!diskp(z, r1)) && wedgep(z, phi1, phi2);
    
    delete(r1), delete(delta), delete(phi1), delete(phi2);

    return ret;
}

bool
in_region_G2 (const num_t z, const num_t alpha, const num_t acc)
{
    bool ret;
    num_t r1, deltat, phi1, phi2;

    r1 = new(num);
    compute_r1(r1, alpha, acc);

    /* Compute deltat */
    deltat = new(num);
    num_set_d(deltat, MIN(M_PI * num_to_d(alpha)/8.0, M_PI * (num_to_d(alpha) + 1.0)/2.0));

    /* Compute phi1, phi2 */
    phi1 = new(num), phi2 = new(num);
    num_set_d(phi1, M_PI * num_to_d(alpha) + num_to_d(deltat));
    num_set_d(phi2, -num_to_d(phi1));

    ret = (!diskp(z, r1)) && wedgep(z, phi1, phi2);

    delete(r1), delete(deltat), delete(phi1), delete(phi2);

    return ret;
}

bool
in_region_G3 (const num_t z, const num_t alpha, const num_t acc)
{
    bool ret;
    num_t r1, delta, deltat, phi1, phi2;

    r1 = new(num);
    compute_r1(r1, alpha, acc);

    /* Compute delta, deltat */
    delta = new(num), deltat = new(num);
    num_set_d(delta, M_PI * num_to_d(alpha)/8.0);
    num_set_d(deltat, MIN(M_PI * num_to_d(alpha)/8.0, M_PI * (num_to_d(alpha) + 1.0)/2.0));

    /* Compute phi1, phi2 */
    phi1 = new(num), phi2 = new(num);
    num_set_d(phi1, M_PI * num_to_d(alpha) - num_to_d(delta));
    num_set_d(phi2, M_PI * num_to_d(alpha) + num_to_d(deltat));

    ret = (!diskp(z, r1)) && closure_wedgep(z, phi1, phi2);

    delete(r1), delete(deltat), delete(phi1), delete(phi2);

    return ret;

    
    /* const double delta = M_PI * alpha/8.0 ; */
    /* const double deltat = MIN(M_PI * alpha/8.0, M_PI * (alpha + 1)/2); */
    /* const double phi1 = M_PI * alpha - delta; */
    /* const double phi2 = M_PI * alpha + deltat; */
    /* const double r1 = compute_r1(alpha, acc); */
    /* const bool ret = (!diskp(z, r1)) && closure_wedgep(z, phi1, phi2); */
    /* return ret; */
}

bool
in_region_G4 (const num_t z, const num_t alpha, const num_t acc)
{
    bool ret;
    num_t r1, delta, deltat, phi1, phi2;

    r1 = new(num);
    compute_r1(r1, alpha, acc);

    /* Compute delta, deltat */
    delta = new(num), deltat = new(num);
    num_set_d(delta, M_PI * num_to_d(alpha)/8.0);
    num_set_d(deltat, MIN(M_PI * num_to_d(alpha)/8.0, M_PI * (num_to_d(alpha) + 1.0)/2.0));

    /* Compute phi1, phi2 */
    phi1 = new(num), phi2 = new(num);
    num_set_d(phi1, -(M_PI * num_to_d(alpha) + num_to_d(deltat)));
    num_set_d(phi2, -M_PI * num_to_d(alpha) + num_to_d(delta));

    ret = (!diskp(z, r1)) && closure_wedgep(z, phi1, phi2);

    delete(r1), delete(delta), delete(deltat), delete(phi1), delete(phi2);

    return ret;
    
    /* const double delta = M_PI * alpha/8.0; */
    /* const double deltat = MIN(M_PI*alpha/8.0, M_PI*(alpha + 1)/2); */
    /* const double phi1 = -(M_PI * alpha +  deltat); */
    /* const double phi2 = delta - M_PI * alpha; */
    /* const double r1 = compute_r1(alpha, acc);     */
    /* const bool ret = (!diskp(z, r1)) && closure_wedgep(z, phi1, phi2); */
    /* return ret; */
}

bool
in_region_G5 (const num_t z, const num_t alpha, const num_t acc)
{
    bool ret;
    num_t r0, r1, phi1, phi2;

    r0 = new(num);
    num_set_d(r0, 0.95);

    r1 = new(num);
    compute_r1(r1, alpha, acc);
   
    /* Compute phi1, phi2 */
    phi1 = new(num), phi2 = new(num);
    num_set_d(phi1, -5.0 * M_PI * num_to_d(alpha)/6.0);
    num_set_d(phi2, +5.0 * M_PI * num_to_d(alpha)/6.0);

    ret = diskp(z, r1) && (closure_wedgep(z, phi1, phi2) && !diskp(z, r0));

    delete(r0), delete(r1), delete(phi1), delete(phi2);

    return ret;    
    /* const double phi1 = -5.0 * M_PI * alpha/6.0; */
    /* const double phi2 = -phi1; */
    /* const double r1 = compute_r1(alpha, acc); */
    /* const bool ret = diskp(z, r1) && (closure_wedgep(z, phi1, phi2) && !diskp(z, R0)); */
    /* return ret; */
}

bool
in_region_G6 (const num_t z, const num_t alpha, const num_t acc)
{
    bool ret;
    num_t r0, r1, phi1, phi2;

    r0 = new(num);
    num_set_d(r0, 0.95);

    r1 = new(num);
    compute_r1(r1, alpha, acc);
   
    /* Compute phi1, phi2 */
    phi1 = new(num), phi2 = new(num);
    num_set_d(phi1, +5.0 * M_PI * num_to_d(alpha)/6.0);
    num_set_d(phi2, -5.0 * M_PI * num_to_d(alpha)/6.0);

    ret = diskp(z, r1) && (closure_wedgep(z, phi1, phi2) && !diskp(z, r0));

    delete(r0), delete(r1), delete(phi1), delete(phi2);

    return ret;
    
    /* const double phi1 = 5.0 * M_PI * alpha/6.0; */
    /* const double phi2 = -5.0 * M_PI * alpha/6.0; */
    /* const double r1 = compute_r1(alpha, acc); */
    /* const bool ret = diskp(z, r1) && (closure_wedgep(z, phi1, phi2) && !diskp(z, R0)); */
    /* return ret; */
}


static bool
_diskp (const num_t z, const num_t r, const bool closure)
{
    num_t absz;
    bool ret;
    
    absz = new(num);
    num_abs(absz, z);
    if (closure == true)
        ret = diskp(z, r) || num_eq(absz, r);
    else
        ret = num_lt(absz, r);
    delete(absz);
    
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
    /* num_t absz; */
    /* bool ret; */
    
    /* absz = new(num); */
    /* num_abs(absz, z); */
    /* ret = diskp(z, r) || num_eq(absz, r); */
    /* delete(absz); */

    /* return ret; */
}

static bool
_wedgep (const num_t z, const num_t phi1, const num_t phi2, const double closure)
{
    bool res;
    num_t argz;

    argz = new(num);
    num_arg(argz, z);
    res = is_between(argz, phi1, phi2, closure);
    delete(argz);
    return res;
}

static bool
wedgep (const num_t z, const num_t phi1, const num_t phi2)
{
    return _wedgep(z, phi1, phi2, false);
    /* bool res; */
    /* num_t argz; */

    /* argz = new(num); */
    /* num_set_d_d(z, x, y); */
    /* num_arg(argz, z); */
    /* res = is_between(num_to_d(argz), phi1, phi2, false); */
    /* delete(argz); */
    /* return res; */
}

static bool
closure_wedgep (const num_t z, const num_t phi1, const num_t phi2)
{
    return _wedgep(z, phi1, phi2, true);
    /* bool res; */
    /* num_t z = new(num); */
    
    /* num_set_d_d(z, x, y); */
    /* num_arg(z, z); */
    /* res = is_between(num_to_d(z), phi1, phi2, true); */
    /* delete(z); */
    /* return res; */
}

/* Checks if _c is between _a and _b */
static bool
is_between (const num_t c, const num_t a, const num_t b, const bool eq)
{
    bool ret;
    num_t _a, _b, _c, n;
    /* double n = 2.0 * M_PI; */
    
    _a = new(num), _b = new(num), _c = new(num), n = new(num);
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

    delete(_a), delete(_b), delete(_c), delete(n);

    
    /* num_set_d(a, fmod(_a, n)); */
    /* num_set_d(b, fmod(_b, n)); */
    /* num_set_d(c, fmod(_c, n)); */

    /* if (num_lt (a, b)) */
    /* { */
    /*     if (eq) */
    /*         ret = (num_le(a, c) && num_le(c, b)) ? true : false; */
    /*     else */
    /*         ret = (num_lt(a, c) && num_lt(c, b)) ? true : false; */
    /* } */
    /* else */
    /* { /\* b < a *\/ */
    /*     /\* if in [b, a] then not in [a, b] *\/ */
    /*     if (eq) */
    /*         ret = (num_le(b, c) && num_le(c, a)) ? false : true; */
    /*     else */
    /*         ret = (num_lt(b, c) && num_lt(c, a)) ? false : true; */
    /* } */

    /* delete(a), delete(b), delete(c); */

    
    
    /* const double n = 2.0 * M_PI; */
    /* const double a = fmod(_a, n); */
    /* const double b = fmod(_b, n); */
    /* const double c = fmod(_c, n); */
    
    /* if (d_lt (a, b)) */
    /* { */
    /*     if (eq) */
    /*         ret = (d_le(a, c) && d_le(c, b)) ? true : false; */
    /*     else */
    /*         ret = (d_lt(a, c) && d_lt(c, b)) ? true : false; */
    /* } */
    /* else */
    /* { /\* b < a *\/ */
    /*     /\* if in [b, a] then not in [a, b] *\/ */
    /*     if (eq) */
    /*         ret = (d_le(b, c) && d_le(c, a)) ? false : true; */
    /*     else */
    /*         ret = (d_lt(b, c) && d_lt(c, a)) ? false : true; */
    /* } */

    return ret;
}


