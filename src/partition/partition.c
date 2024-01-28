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
 * @file partition.c
 * @brief Implementation for the routines regarding the partition of the complex plane.
 */

#include "partition.h"
#include "flintutils.h"

#include <math.h>
#include <stdbool.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

static bool diskp          (const acb_t z, const arb_t r);
static bool closure_diskp  (const acb_t z, const arb_t r0);
static bool wedgep         (const acb_t z, const arb_t phi1, const arb_t phi2);
static bool closure_wedgep (const acb_t z, const arb_t phi1, const arb_t phi2);
static bool is_between     (const arb_t c, const arb_t a, const arb_t b, const bool eq);

/* Equation (4.21) */
static void
compute_r1 (arb_t r1, const arb_t alpha, const arb_t acc)
{
    arb_t one, pi, c0, n, d, aux, two;

    arb_init(one);
    arb_init(pi);
    arb_init(c0);
    arb_init(n);
    arb_init(d);
    arb_init(aux);
    arb_init(two);

    arb_one(one);
    arb_const_pi(pi, PREC);
    arb_set_d(two, 2.0);

    /* Compute C0 */
    /* pow(1.3, 1.0 - alpha) */
    arb_set_d(n, 1.3);
    arb_sub(aux, one, alpha, PREC);
    arb_pow(n, n, aux, PREC);
    /* M_PI * sin(M_PI * alpha) */
    arb_mul(d, pi, alpha, PREC);
    arb_sin(d, d, PREC);
    arb_mul(d, d, pi, PREC);
    arb_div(c0, n, d, PREC);

    arb_div(aux, acc, c0, PREC);
    arb_log(aux, aux, PREC);
    arb_mul(aux, aux, two, PREC);
    arb_neg(aux, aux);
    arb_pow(r1, aux, alpha, PREC);
    
    arb_clear(pi);
    arb_clear(one);
    arb_clear(c0);
    arb_clear(n);
    arb_clear(d);
    arb_clear(aux);
}

bool
in_region_G0 (const acb_t z)
{
    bool ret;
    arb_t r0;

    arb_init(r0);
    arb_set_d(r0, 0.95);
    ret = closure_diskp(z, r0);
    arb_clear(r0);

    return ret;
}

bool
in_region_G1 (const acb_t z, const arb_t alpha, const arb_t acc)
{
    bool ret;
    arb_t r1, delta, phi1, phi2;

    arb_init(r1);
    compute_r1(r1, alpha, acc);
    
    /* Compute delta */
    arb_init(delta);
    arb_set_d(delta, M_PI * arbtod(alpha)/8.0);
    
    /* Compute phi1, phi2 */
    arb_init(phi1);
    arb_init(phi2);
    arb_set_d(phi1, - M_PI * arbtod(alpha) + arbtod(delta) );
    arb_set_d(phi2, + M_PI * arbtod(alpha) - arbtod(delta));

    ret = (!diskp(z, r1)) && wedgep(z, phi1, phi2);
    
    arb_clear(r1);
    arb_clear(delta);
    arb_clear(phi1);
    arb_clear(phi2);

    return ret;
}

bool
in_region_G2 (const acb_t z, const arb_t alpha, const arb_t acc)
{
    bool ret;
    arb_t r1, deltat, phi1, phi2;

    arb_init(r1);
    compute_r1(r1, alpha, acc);

    /* Compute deltat */
    arb_init(deltat);
    arb_set_d(deltat, MIN(M_PI * arbtod(alpha)/8.0, M_PI * (arbtod(alpha) + 1.0)/2.0));

    /* Compute phi1, phi2 */
    arb_init(phi1);
    arb_init(phi2);
    arb_set_d(phi1, M_PI * arbtod(alpha) + arbtod(deltat));
    arb_set_d(phi2, -arbtod(phi1));

    ret = (!diskp(z, r1)) && wedgep(z, phi1, phi2);

    arb_clear(r1);
    arb_clear(deltat);
    arb_clear(phi1);
    arb_clear(phi2);

    return ret;
}

bool
in_region_G3 (const acb_t z, const arb_t alpha, const arb_t acc)
{
    bool ret;
    arb_t r1, delta, deltat, phi1, phi2;

    arb_init(r1);
    compute_r1(r1, alpha, acc);

    /* Compute delta, deltat */
    arb_init(delta);
    arb_init(deltat);
    arb_set_d(delta, M_PI * arbtod(alpha)/8.0);
    arb_set_d(deltat, MIN(M_PI * arbtod(alpha)/8.0, M_PI * (arbtod(alpha) + 1.0)/2.0));

    /* Compute phi1, phi2 */
    arb_init(phi1);
    arb_init(phi2);
    arb_set_d(phi1, M_PI * arbtod(alpha) - arbtod(delta));
    arb_set_d(phi2, M_PI * arbtod(alpha) + arbtod(deltat));

    ret = (!diskp(z, r1)) && closure_wedgep(z, phi1, phi2);

    arb_clear(r1);
    arb_clear(deltat);
    arb_clear(phi1);
    arb_clear(phi2);

    return ret;
}

bool
in_region_G4 (const acb_t z, const arb_t alpha, const arb_t acc)
{
    bool ret;
    arb_t r1, delta, deltat, phi1, phi2;

    arb_init(r1);
    compute_r1(r1, alpha, acc);

    /* Compute delta, deltat */
    arb_init(delta);
    arb_init(deltat);
    arb_set_d(delta, M_PI * arbtod(alpha)/8.0);
    arb_set_d(deltat, MIN(M_PI * arbtod(alpha)/8.0, M_PI * (arbtod(alpha) + 1.0)/2.0));

    /* Compute phi1, phi2 */
    arb_init(phi1);
    arb_init(phi2);
    arb_set_d(phi1, -(M_PI * arbtod(alpha) + arbtod(deltat)));
    arb_set_d(phi2, -M_PI * arbtod(alpha) + arbtod(delta));

    ret = (!diskp(z, r1)) && closure_wedgep(z, phi1, phi2);

    arb_clear(r1);
    arb_clear(delta);
    arb_clear(deltat);
    arb_clear(phi1);
    arb_clear(phi2);

    return ret;
}

bool
in_region_G5 (const acb_t z, const arb_t alpha, const arb_t acc)
{
    bool ret;
    arb_t r0, r1, phi1, phi2;

    arb_init(r0);
    arb_set_d(r0, 0.95);

    arb_init(r1);
    compute_r1(r1, alpha, acc);
   
    /* Compute phi1, phi2 */
    arb_init(phi1);
    arb_init(phi2);
    arb_set_d(phi1, -5.0 * M_PI * arbtod(alpha)/6.0);
    arb_set_d(phi2, +5.0 * M_PI * arbtod(alpha)/6.0);

    ret = diskp(z, r1) && (closure_wedgep(z, phi1, phi2) && !diskp(z, r0));

    arb_clear(r0);
    arb_clear(r1);
    arb_clear(phi1);
    arb_clear(phi2);

    return ret;    
}

bool
in_region_G6 (const acb_t z, const arb_t alpha, const arb_t acc)
{
    bool ret;
    arb_t r0, r1, phi1, phi2;

    arb_init(r0);
    arb_set_d(r0, 0.95);

    arb_init(r1);
    compute_r1(r1, alpha, acc);
   
    /* Compute phi1, phi2 */
    arb_init(phi1);
    arb_init(phi2);
    arb_set_d(phi1, +5.0 * M_PI * arbtod(alpha)/6.0);
    arb_set_d(phi2, -5.0 * M_PI * arbtod(alpha)/6.0);

    ret = diskp(z, r1) && (closure_wedgep(z, phi1, phi2) && !diskp(z, r0));

    arb_clear(r0);
    arb_clear(r1);
    arb_clear(phi1);
    arb_clear(phi2);

    return ret;
}


static bool
_diskp (const acb_t z, const arb_t r, const bool closure)
{
    arb_t absz;
    bool ret;
    
    arb_init(absz);
    acb_abs(absz, z, PREC);
    if (closure == true)
        ret = diskp(z, r) || arb_eq(absz, r);
    else
        ret = arb_lt(absz, r);

    arb_clear(absz);
    
    return ret;
}

static bool
diskp (const acb_t z, const arb_t r)
{
    return _diskp(z, r, false);
}

static bool
closure_diskp (const acb_t z, const arb_t r)
{
    return _diskp(z, r, true);
}

static bool
_wedgep (const acb_t z, const arb_t phi1, const arb_t phi2, const double closure)
{
    bool res;
    arb_t argz;

    arb_init(argz);

    acb_arg(argz, z, PREC);
    res = is_between(argz, phi1, phi2, closure);

    arb_clear(argz);
    
    return res;
}

static bool
wedgep (const acb_t z, const arb_t phi1, const arb_t phi2)
{
    return _wedgep(z, phi1, phi2, false);
}

static bool
closure_wedgep (const acb_t z, const arb_t phi1, const arb_t phi2)
{
    return _wedgep(z, phi1, phi2, true);
}

/* Checks if _c is between _a and _b */
static bool
is_between (const arb_t c, const arb_t a, const arb_t b, const bool eq)
{
    bool ret;
    arb_t _a, _b, _c, n;
    
    arb_init(_a);
    arb_init(_b);
    arb_init(_c);
    arb_init(n);
    arb_set_d(n, 2.0 * M_PI);

    arb_fmod(_a, a, n);
    arb_fmod(_b, b, n);
    arb_fmod(_c, c, n);

    if (arb_lt(_a, _b))
    {
        if (eq)
            ret = (arb_le(_a, _c) && arb_le(_c, _b)) ? true : false;
        else
            ret = (arb_lt(_a, _c) && arb_lt(_c, _b)) ? true : false;
    }
    else
    { /* b < a */
        /* if in [b, a] then not in [a, b] */
        if (eq)
            ret = (arb_le(_b, _c) && arb_le(_c, _a)) ? false : true;
        else
            ret = (arb_lt(_b, _c) && arb_lt(_c, _a)) ? false : true;
    }

    arb_clear(_a);
    arb_clear(_b);
    arb_clear(_c);
    arb_clear(n);

    return ret;
}

#undef MIN
#undef MAX
