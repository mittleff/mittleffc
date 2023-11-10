#include "partition.h"

#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_math.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define R0 0.95
#define GSL_EPSILON_FCMP 1e-8

static bool
d_eq (const double x, const double y)
{
    int ret = gsl_fcmp(x, y, GSL_EPSILON_FCMP);
    return (ret == 0) ? true : false;
}

static bool
d_lt (const double x, const double y)
{
    int ret = gsl_fcmp(x, y, GSL_EPSILON_FCMP);
    return (ret == -1) ? true : false;
}

static bool
d_gt (const double x, const double y)
{
    int ret = gsl_fcmp(x, y, GSL_EPSILON_FCMP);
    return (ret == 1) ? true : false;
}

static bool
d_le (const double x, const double y)
{
    return d_lt(x, y) || d_eq(x, y);
}

static bool
d_ge (const double x, const double y)
{
    return d_gt(x, y) || d_eq(x, y);
}

static bool diskp (const double x, const double y, const double r);
static bool closure_diskp (const double x, const double y, const double r0);
static bool wedgep (const double x, const double y, const double phi1, const double phi2);
static bool is_between (const double c, const double a, const double b, const bool eq);
static bool closure_wedgep (const double x, const double y, const double phi1, const double phi2);

static double
compute_r1 (const double alpha, const double acc)
{
    /* Compute r1 */
    const double C0 = pow(1.3, 1 - alpha)/(M_PI * sin(M_PI * alpha)); /* Equation (5.3) */
    const double r1 = pow(-2.0 * log(acc/C0), alpha); /* Equation (4.21) */

    return r1;
}

bool
in_region_G0 (const double x, const double y)
{
    return closure_diskp(x, y, R0);
}

bool
in_region_G1 (const double x, const double y, const double alpha, const double acc)
{
    const double delta = M_PI * alpha/8.0;
    const double phi1 = delta - M_PI * alpha;
    const double phi2 = M_PI * alpha - delta;
    const double r1 = compute_r1(alpha, acc);    
    const bool res = (!diskp(x, y, r1)) && wedgep(x, y, phi1, phi2);
    return res;
}

bool
in_region_G2 (const double x, const double y, const double alpha, const double acc)
{
    const double deltat = MIN(M_PI*alpha/8.0, M_PI*(alpha + 1)/2);
    const double phi1 = deltat + M_PI * alpha;
    const double phi2 = -(M_PI * alpha + deltat);
    const double r1 = compute_r1(alpha, acc);
    const bool ret = (!diskp(x, y, r1)) && wedgep(x, y, phi1, phi2);    
    return ret;
}

bool
in_region_G3 (const double x, const double y, const double alpha, const double acc)
{
    const double delta = M_PI * alpha/8.0 ;
    const double deltat = MIN(M_PI * alpha/8.0, M_PI * (alpha + 1)/2);
    const double phi1 = M_PI * alpha - delta;
    const double phi2 = M_PI * alpha + deltat;
    const double r1 = compute_r1(alpha, acc);
    const bool ret = (!diskp(x, y, r1)) && closure_wedgep(x, y, phi1, phi2);
    return ret;
}

bool
in_region_G4 (const double x, const double y, const double alpha, const double acc)
{
    const double delta = M_PI * alpha/8.0;
    const double deltat = MIN(M_PI*alpha/8.0, M_PI*(alpha + 1)/2);
    const double phi1 = -(M_PI * alpha +  deltat);
    const double phi2 = delta - M_PI * alpha;
    const double r1 = compute_r1(alpha, acc);    
    const bool ret = (!diskp(x, y, r1)) && closure_wedgep(x, y, phi1, phi2);
    return ret;
}

bool
in_region_G5 (const double x, const double y, const double alpha, const double acc)
{
    const double phi1 = -5.0 * M_PI * alpha/6.0;
    const double phi2 = 5.0 * M_PI * alpha/6.0;
    const double r1 = compute_r1(alpha, acc);
    const bool ret = diskp(x, y, r1) && (closure_wedgep(x, y, phi1, phi2) && !diskp(x, y, R0));
    return ret;
}

bool
in_region_G6 (const double x, const double y, const double alpha, const double acc)
{
    const double phi1 = 5.0 * M_PI * alpha/6.0;
    const double phi2 = -5.0 * M_PI * alpha/6.0;
    const double r1 = compute_r1(alpha, acc);
    const bool ret = diskp(x, y, r1) && (closure_wedgep(x, y, phi1, phi2) && !diskp(x, y, R0));
    return ret;
}


static bool
diskp (const double x, const double y, const double r)
{
    const double abs_z = sqrt(x*x + y*y);
    return d_lt(abs_z, r);
}

static bool
closure_diskp (const double x, const double y, const double r)
{
    const double abs_z = sqrt(x*x + y*y);
    return diskp(x, y, r) || d_eq(abs_z, r);
}

static bool
wedgep (const double x, const double y, const double phi1, const double phi2)
{
    const double arg_z = atan(y/x); 
    return is_between(arg_z, phi1, phi2, false);
}

static bool
closure_wedgep (const double x, const double y, const double phi1, const double phi2)
{
    const double arg_z = atan(y/x);    
    return is_between(arg_z, phi1, phi2, true);
}

static bool
is_between (const double _c, const double _a, const double _b, const bool eq)
{
    bool ret = true;
    
    const double n = 2.0 * M_PI;
    const double a = fmod(_a, n);
    const double b = fmod(_b, n);
    const double c = fmod(_c, n);
    
    if (d_lt (a, b))
    {
        if (eq)
            ret = (d_le(a, c) && d_le(c, b)) ? true : false;
        else
            ret = (d_lt(a, c) && d_lt(c, b)) ? true : false;
    }
    else
    { /* b < a */
        /* if in [b, a] then not in [a, b] */
        if (eq)
            ret = (d_le(b, c) && d_le(c, a)) ? false : true;
        else
            ret = (d_lt(b, c) && d_lt(c, a)) ? false : true;
    }

    return ret;
}


