#include "partition.h"

#include "log.h"

#include "new.h"
#include "num.h"

#include <math.h>


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define R0 0.95

static bool diskp (const num_t z, const num_t r);
static bool closure_diskp (const num_t z, const num_t r);
static bool wedgep (const num_t z, const num_t phi1, const num_t phi2);
static bool is_between (const num_t _c, const num_t _a, const num_t _b, const bool eq);
static bool closure_wedgep (const num_t z, const num_t phi1, const num_t phi2);

bool
in_region_G0 (const num_t z)
{
    /* log_trace("[%s] called with parameters: z=(%+.5e, %+.5e)", */
    /*           __func__, num_real_d(z), num_imag_d(z)); */
    return closure_diskp(z, new(num, R0, 0.0));
}

bool
in_region_G1 (const num_t z, const num_t alpha, const num_t acc)
{
    const double _alpha = num_to_d(alpha);
    const num_t pi = new(num, M_PI, 0.0);
    const num_t delta = new(num, M_PI * _alpha/8.0, 0.0);
    const num_t phi1 = num_sub(delta, num_mul(pi, alpha));
    const num_t phi2 = num_sub(num_mul (pi, alpha), delta);

    /* Compute r1 */
    const double C0 = pow(1.3, 1 - num_to_d(alpha))/(M_PI * sin(M_PI * num_to_d(alpha))); /* Equation (5.3) */
    const num_t r1 = new(num, pow(-2.0 * log(num_to_d(acc)/C0), num_to_d(alpha)), 0.0); /* Equation (4.21) */
    
    const bool ret = (!diskp(z, r1)) && wedgep(z, phi1, phi2);

    delete(pi); delete(phi1); delete(phi2);
    return ret;
}

bool
in_region_G2 (const num_t z, const num_t alpha, const num_t acc)
{
    const double _alpha = num_to_d(alpha);
    const num_t pi = new(num, M_PI, 0.0);
    const num_t deltat = new(num, MIN(M_PI*_alpha/8.0, M_PI*(_alpha + 1)/2),   0.0);
    const num_t phi1 = num_add(deltat, num_mul(pi, alpha));
    const num_t phi2 = num_neg(num_add(num_mul(pi, alpha), deltat));

    /* Compute r1 */
    const double C0 = pow(1.3, 1 - num_to_d(alpha))/(M_PI * sin(M_PI * num_to_d(alpha))); /* Equation (5.3) */
    const num_t r1 = new(num, pow(-2.0 * log(num_to_d(acc)/C0), num_to_d(alpha)), 0.0); /* Equation (4.21) */
    
    const bool ret = (!diskp(z, r1)) && wedgep(z, phi1, phi2);
    
    delete(pi); delete(phi1); delete(phi2);
    
    return ret;
}

bool
in_region_G3 (const num_t z, const num_t alpha, const num_t acc)
{
    const double _alpha = num_to_d(alpha);
    const num_t pi = new(num, M_PI, 0.0);
    const num_t delta = new(num, M_PI * _alpha/8.0, 0.0);
    const num_t deltat = new(num, MIN(M_PI*_alpha/8.0, M_PI*(_alpha + 1)/2),   0.0);
    const num_t phi1 = num_sub(num_mul(pi, alpha), delta);
    const num_t phi2 = num_add(num_mul(pi, alpha), deltat);

    /* Compute r1 */
    const double C0 = pow(1.3, 1 - num_to_d(alpha))/(M_PI * sin(M_PI * num_to_d(alpha))); /* Equation (5.3) */
    const num_t r1 = new(num, pow(-2.0 * log(num_to_d(acc)/C0), num_to_d(alpha)), 0.0); /* Equation (4.21) */

    const bool ret = (!diskp(z, r1)) && closure_wedgep(z, phi1, phi2);

    delete(pi); delete(phi1); delete(phi2);

    return ret;
}

bool
in_region_G4 (const num_t z, const num_t alpha, const num_t acc)
{
    const double _alpha = num_to_d(alpha);
    const num_t pi = new(num, M_PI, 0.0);
    const num_t delta = new(num, M_PI * _alpha/8.0, 0.0);
    const num_t deltat = new(num, MIN(M_PI*_alpha/8.0, M_PI*(_alpha + 1)/2),   0.0);
    const num_t phi1 = num_neg(num_add(num_mul (pi, alpha), deltat));
    const num_t phi2 = num_sub(delta, num_mul(pi, alpha));

     /* Compute r1 */
    const double C0 = pow(1.3, 1 - num_to_d(alpha))/(M_PI * sin(M_PI * num_to_d(alpha))); /* Equation (5.3) */
    const num_t r1 = new(num, pow(-2.0 * log(num_to_d(acc)/C0), num_to_d(alpha)), 0.0); /* Equation (4.21) */
    
    const bool ret = (!diskp(z, r1)) && closure_wedgep(z, phi1, phi2);

    delete(pi); delete(phi1); delete(phi2);

    return ret;
}

bool
in_region_G5 (const num_t z, const num_t alpha, const num_t acc)
{
    const num_t phi1 = new(num, -5.0 * M_PI * num_to_d(alpha)/6.0, 0.0);
    const num_t phi2 = new(num,  5.0 * M_PI * num_to_d(alpha)/6.0, 0.0);

    /* Compute r1 */
    const double C0 = pow(1.3, 1 - num_to_d(alpha))/(M_PI * sin(M_PI * num_to_d(alpha))); /* Equation (5.3) */
    const num_t r1 = new(num, pow(-2.0 * log(num_to_d(acc)/C0), num_to_d(alpha)), 0.0); /* Equation (4.21) */
    
    const bool ret = diskp(z, r1) && (closure_wedgep(z, phi1, phi2) && !diskp(z, new(num, R0, 0.0)));
    delete(phi1); delete(phi2);
    
    return ret;
}

bool
in_region_G6 (const num_t z, const num_t alpha, const num_t acc)
{
    const num_t phi1 = new (num,  5.0 * M_PI * num_to_d(alpha)/6.0, 0.0);
    const num_t phi2 = new (num, -5.0 * M_PI * num_to_d(alpha)/6.0, 0.0);
    /* Compute r1 */
    const double C0 = pow(1.3, 1 - num_to_d(alpha))/(M_PI * sin(M_PI * num_to_d(alpha))); /* Equation (5.3) */
    const num_t r1 = new(num, pow(-2.0 * log(num_to_d(acc)/C0), num_to_d(alpha)), 0.0); /* Equation (4.21) */
    const bool ret = diskp(z, r1) && (closure_wedgep(z, phi1, phi2) && !diskp(z, new(num, R0, 0.0)));
    delete(phi1); delete(phi2);
    return ret;
}


static bool
diskp (const num_t z, const num_t r)
{
    return num_lt(num_abs(z), r);
}

static bool
closure_diskp (const num_t z, const num_t r)
{
    return diskp(z, r) || num_eq(num_abs(z), r);
}

static bool
wedgep (const num_t z, const num_t phi1, const num_t phi2)
{
    return is_between(num_arg (z), phi1, phi2, false);
}

static bool
is_between (const num_t _c, const num_t _a, const num_t _b, const bool eq)
{
    bool ret = true;
    
    const num_t n = new(num, 2.0 * M_PI, 0.0);
    const num_t a = num_fmod (_a, n);
    const num_t b = num_fmod (_b, n);
    const num_t c = num_fmod (_c, n);
    
    if (num_lt (a, b))
    {
        if (eq)
            ret = (num_le(a, c) && num_le(c, b)) ? true : false;
        else
            ret = (num_lt(a, c) && num_lt(c, b)) ? true : false;
    }
    else
    { /* b < a */
        /* if in [b, a] then not in [a, b] */
        if (eq)
            ret = (num_le(b, c) && num_le(c, a)) ? false : true;
        else
            ret = (num_lt(b, c) && num_lt(c, a)) ? false : true;
    }

    delete(n); delete(a); delete(b); delete(c);
    return ret;
}

static bool
closure_wedgep (const num_t z, const num_t phi1, const num_t phi2)
{
    return is_between(num_arg (z), phi1, phi2, true);
}
