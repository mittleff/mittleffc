#include <complex.h>
#include <math.h>

#include <flint/arb.h>
#include <flint/acb.h>

#include "include/partition.h"
#include "include/utils.h"


/************************************************************************/
/* Region checking functions for the Mittag-Leffler algorithm           */
/* These functions determine which computational branch to use based on */
/* the complex argument z and parameter alpha                           */
/************************************************************************/

int
region_flint (const acb_t z, const arb_t alpha, const arb_t taylor_radius, slong prec);

static int
in_region_G1_flint (const acb_t z, const arb_t alpha, const arb_t taylor_radius, slong prec);

static int
in_region_G2_flint (const acb_t z, const arb_t alpha, const arb_t taylor_radius, slong prec);

static int
in_region_G3_flint (const acb_t z, const arb_t alpha, const arb_t taylor_radius, slong prec);

static int
in_region_G4_flint (const acb_t z, const arb_t alpha, const arb_t taylor_radius, slong prec);

static int
in_region_G5_flint (const acb_t z, const arb_t alpha, const arb_t taylor_radius, slong prec);

static int
in_region_G6_flint(const acb_t z, const arb_t alpha, const arb_t taylor_radius, slong prec);

unsigned short int
region (double complex z, double alpha, double taylor_radius, slong prec)
{
    
    int result;
    arb_t alpha_arb, taylor_radius_arb;
    acb_t z_acb;
    
    arb_init(alpha_arb);
    arb_init(taylor_radius_arb);
    acb_init(z_acb);
    
    arb_set_d(alpha_arb, alpha);
    arb_set_d(taylor_radius_arb, taylor_radius);
    acb_set_d_d(z_acb, creal(z), cimag(z));
    
    result = region_flint(z_acb, alpha_arb, taylor_radius_arb, prec);

    arb_clear(alpha_arb);
    arb_clear(taylor_radius_arb);
    acb_clear(z_acb);
    
    return result;
}

int
region_flint (const acb_t z, const arb_t alpha, const arb_t taylor_radius, slong prec)
{
    arb_t mag_z, r1;
    int result = -1; // Default error value
    
    // Initialize temporary variables
    arb_init(mag_z);
    arb_init(r1);
    
    // Compute magnitude of z
    acb_abs(mag_z, z, prec);
    
    // Check if |z| <= taylor_radius
    if (arb_le(mag_z, taylor_radius)) {
        result = 0;
        goto cleanup;
    }
    
    // Compute r1
    compute_r1_flint(r1, alpha, prec);
    
    // Check if |z| >= r1
    if (arb_ge(mag_z, r1)) {
        // Check regions 1-4
        if (in_region_G1_flint(z, alpha, taylor_radius, prec)) {
            result = 1;
        } else if (in_region_G2_flint(z, alpha, taylor_radius, prec)) {
            result = 2;
        } else if (in_region_G3_flint(z, alpha, taylor_radius, prec)) {
            result = 3;
        } else if (in_region_G4_flint(z, alpha, taylor_radius, prec)) {
            result = 4;
        }
    } else {
        // Check regions 5-6
        if (in_region_G5_flint(z, alpha, taylor_radius, prec)) {
            result = 5;
        } else if (in_region_G6_flint(z, alpha, taylor_radius, prec)) {
            result = 6;
        }
    }
    
cleanup:
    // Clear temporary variables
    arb_clear(mag_z);
    arb_clear(r1);
    
    return result;
}


// Compute cross product for angular comparison: sin(a - b)
static void
cross_product (arb_t result, const arb_t x, const arb_t y, slong prec)
{
    arb_sub(result, x, y, prec);
    arb_sin(result, result, prec);
}

// Check if angle a is clockwise from angle b
static int
clockwise_p (const arb_t a, const arb_t b, int closed, slong prec)
{
    int result;
    arb_t cross;
    
    arb_init(cross);
    cross_product(cross, a, b, prec);
    
    if (closed) {
        result = arb_is_nonnegative(cross);
    } else {
        result = arb_is_positive(cross);
    }
    
    arb_clear(cross);
    
    return result;
}

// Check if angle c is between angles phi1 and phi2
static int
between_angles (const arb_t c, const arb_t phi1, const arb_t phi2, int closed, slong prec)
{
    if (clockwise_p(phi1, phi2, 0, prec)) {
        // If phi1 is clockwise from phi2, then NOT between phi2 and phi1
        return !between_angles(c, phi2, phi1, closed, prec);
    } else {
        // Check: NOT clockwise(phi1, c) AND clockwise(phi2, c)
        return !clockwise_p(phi1, c, closed, prec) && clockwise_p(phi2, c, closed, prec);
    }
}

// Check if complex number z is in an open wedge
static int
open_wedge_p (const acb_t z, const arb_t phi1, const arb_t phi2, slong prec)
{
    int result;    
    arb_t arg_z;
    
    arb_init(arg_z);
    acb_arg(arg_z, z, prec);
    
    result = between_angles(arg_z, phi1, phi2, 0, prec); // closed = 0 for open wedge
    
    arb_clear(arg_z);    
    return result;
}

// Check if complex number z is in a closed wedge
static int
closed_wedge_p (const acb_t z, const arb_t phi1, const arb_t phi2, slong prec)
{
    int result;
    arb_t arg_z;
    
    arb_init(arg_z);
    acb_arg(arg_z, z, prec);
    
    result = between_angles(arg_z, phi1, phi2, 1, prec); // closed = 1 for closed wedge
    
    arb_clear(arg_z);
    return result;
}

static int
in_region_G0_flint(const acb_t z, const arb_t alpha, const arb_t taylor_radius, slong prec) {
    int result;
    arb_t abs_z;
    
    arb_init(abs_z);
    acb_abs(abs_z, z, prec);
    
    result = arb_le(abs_z, taylor_radius);

    arb_clear(abs_z);
    
    return result;
}

// Region checking functions implementing exact Scheme algorithm
static int
in_region_G1_flint (const acb_t z, const arb_t alpha, const arb_t taylor_radius, slong prec)
{
    arb_t delta, phi1, phi2, pi_val, temp;
    
    arb_init(delta);
    arb_init(phi1);
    arb_init(phi2);
    arb_init(pi_val);
    arb_init(temp);
    
    arb_const_pi(pi_val, prec);
    arb_mul(temp, pi_val, alpha, prec);
    arb_set_d(delta, 0.125);
    arb_mul(delta, delta, temp, prec); // delta = pi * alpha * 0.125
    
    arb_neg(phi1, temp);
    arb_add(phi1, phi1, delta, prec); // phi1 = -pi * alpha + delta
    
    arb_sub(phi2, temp, delta, prec); // phi2 = pi * alpha - delta
    
    int result = open_wedge_p(z, phi1, phi2, prec);
    
    arb_clear(delta);
    arb_clear(phi1);
    arb_clear(phi2);
    arb_clear(pi_val);
    arb_clear(temp);
    
    return result;
}

static int
in_region_G2_flint (const acb_t z, const arb_t alpha, const arb_t taylor_radius, slong prec)
{
    arb_t tau, phi1, phi2, pi_val, temp1, temp2, temp3, one_plus_alpha;
    
    arb_init(tau);
    arb_init(phi1);
    arb_init(phi2);
    arb_init(pi_val);
    arb_init(temp1);
    arb_init(temp2);
    arb_init(temp3);
    arb_init(one_plus_alpha);
    
    arb_const_pi(pi_val, prec);
    
    // tau = min(0.125 * pi * alpha, 0.5 * pi * (1 + alpha))
    arb_mul(temp1, pi_val, alpha, prec); // pi * alpha
    arb_set_d(temp2, 0.125);
    arb_mul(temp2, temp2, temp1, prec); // 0.125 * pi * alpha
    
    arb_t one;
    arb_init(one);
    arb_one(one);
    arb_add(one_plus_alpha, one, alpha, prec); // 1 + alpha
    arb_set_d(temp3, 0.5);
    arb_mul(temp3, temp3, pi_val, prec);
    arb_mul(temp3, temp3, one_plus_alpha, prec); // 0.5 * pi * (1 + alpha)
    arb_clear(one);
    
    arb_min(tau, temp2, temp3, prec);
    
    arb_add(phi1, temp1, tau, prec); // phi1 = pi * alpha + tau
    arb_neg(phi2, phi1); // phi2 = -phi1
    
    int result = open_wedge_p(z, phi1, phi2, prec);
    
    arb_clear(tau);
    arb_clear(phi1);
    arb_clear(phi2);
    arb_clear(pi_val);
    arb_clear(temp1);
    arb_clear(temp2);
    arb_clear(temp3);
    arb_clear(one_plus_alpha);
    
    return result;
}

static int
in_region_G3_flint (const acb_t z, const arb_t alpha, const arb_t taylor_radius, slong prec)
{
    arb_t delta, tau, phi1, phi2, pi_val, temp1, temp2, temp3, one_plus_alpha;
    
    arb_init(delta);
    arb_init(tau);
    arb_init(phi1);
    arb_init(phi2);
    arb_init(pi_val);
    arb_init(temp1);
    arb_init(temp2);
    arb_init(temp3);
    arb_init(one_plus_alpha);
    
    arb_const_pi(pi_val, prec);
    arb_mul(temp1, pi_val, alpha, prec); // pi * alpha
    arb_set_d(delta, 0.125);
    arb_mul(delta, delta, temp1, prec); // delta = pi * alpha * 0.125
    
    // tau = min(0.125 * pi * alpha, 0.5 * pi * (1 + alpha))
    arb_t one;
    arb_init(one);
    arb_one(one);
    arb_add(one_plus_alpha, one, alpha, prec); // 1 + alpha
    arb_set_d(temp3, 0.5);
    arb_mul(temp3, temp3, pi_val, prec);
    arb_mul(temp3, temp3, one_plus_alpha, prec); // 0.5 * pi * (1 + alpha)
    arb_clear(one);
    arb_min(tau, delta, temp3, prec);
    
    arb_sub(phi1, temp1, delta, prec); // phi1 = pi * alpha - delta
    arb_add(phi2, temp1, tau, prec); // phi2 = pi * alpha + tau
    
    int result = closed_wedge_p(z, phi1, phi2, prec);
    
    arb_clear(delta);
    arb_clear(tau);
    arb_clear(phi1);
    arb_clear(phi2);
    arb_clear(pi_val);
    arb_clear(temp1);
    arb_clear(temp2);
    arb_clear(temp3);
    arb_clear(one_plus_alpha);
    
    return result;
}

static int
in_region_G4_flint (const acb_t z, const arb_t alpha, const arb_t taylor_radius, slong prec)
{
    arb_t delta, tau, phi1, phi2, pi_val, temp1, temp2, temp3, one_plus_alpha;
    
    arb_init(delta);
    arb_init(tau);
    arb_init(phi1);
    arb_init(phi2);
    arb_init(pi_val);
    arb_init(temp1);
    arb_init(temp2);
    arb_init(temp3);
    arb_init(one_plus_alpha);
    
    arb_const_pi(pi_val, prec);
    arb_mul(temp1, pi_val, alpha, prec); // pi * alpha
    arb_set_d(delta, 0.125);
    arb_mul(delta, delta, temp1, prec); // delta = pi * alpha * 0.125
    
    // tau = min(0.125 * pi * alpha, 0.5 * pi * (1 + alpha))
    arb_t one;
    arb_init(one);
    arb_one(one);
    arb_add(one_plus_alpha, one, alpha, prec); // 1 + alpha
    arb_set_d(temp3, 0.5);
    arb_mul(temp3, temp3, pi_val, prec);
    arb_mul(temp3, temp3, one_plus_alpha, prec); // 0.5 * pi * (1 + alpha)
    arb_clear(one);
    arb_min(tau, delta, temp3, prec);
    
    arb_neg(temp2, temp1); // -pi * alpha
    arb_sub(phi1, temp2, tau, prec); // phi1 = -pi * alpha - tau
    arb_add(phi2, temp2, delta, prec); // phi2 = -pi * alpha + delta
    
    int result = closed_wedge_p(z, phi1, phi2, prec);
    
    arb_clear(delta);
    arb_clear(tau);
    arb_clear(phi1);
    arb_clear(phi2);
    arb_clear(pi_val);
    arb_clear(temp1);
    arb_clear(temp2);
    arb_clear(temp3);
    arb_clear(one_plus_alpha);
    
    return result;
}

static int
in_region_G5_flint (const acb_t z, const arb_t alpha, const arb_t taylor_radius, slong prec)
{
    arb_t phi1, phi2, pi_val, abs_z, temp1, temp2;
    
    arb_init(phi1);
    arb_init(phi2);
    arb_init(pi_val);
    arb_init(abs_z);
    arb_init(temp1);
    arb_init(temp2);
    
    arb_const_pi(pi_val, prec);
    arb_mul(temp1, pi_val, alpha, prec); // pi * alpha
    
    // phi1 = -(5/6) * pi * alpha, phi2 = (5/6) * pi * alpha
    arb_set_d(temp2, -5.0/6.0);
    arb_mul(phi1, temp2, temp1, prec); // phi1 = -(5/6) * pi * alpha
    
    arb_set_d(temp2, 5.0/6.0);
    arb_mul(phi2, temp2, temp1, prec); // phi2 = (5/6) * pi * alpha
    
    acb_abs(abs_z, z, prec);
    //arb_set_d(taylor_radius, TAYLOR_RADIUS_DBL);
    
    int in_wedge = closed_wedge_p(z, phi1, phi2, prec);
    int magnitude_check = !in_region_G0_flint(z, alpha, taylor_radius, prec);//arb_gt(abs_z, taylor_radius);
    
    int result = in_wedge && magnitude_check;
    
    arb_clear(phi1);
    arb_clear(phi2);
    arb_clear(pi_val);
    arb_clear(abs_z);
    arb_clear(temp1);
    arb_clear(temp2);
    
    return result;
}

static int
in_region_G6_flint(const acb_t z, const arb_t alpha, const arb_t taylor_radius, slong prec) {
    arb_t phi1, phi2, pi_val, abs_z, temp1, temp2;
    
    arb_init(phi1);
    arb_init(phi2);
    arb_init(pi_val);
    arb_init(abs_z);
    arb_init(temp1);
    arb_init(temp2);
    
    arb_const_pi(pi_val, prec);
    arb_mul(temp1, pi_val, alpha, prec); // pi * alpha
    
    // phi1 = (5/6) * pi * alpha, phi2 = -(5/6) * pi * alpha
    arb_set_d(temp2, 5.0/6.0);
    arb_mul(phi1, temp2, temp1, prec); // phi1 = (5/6) * pi * alpha
    
    arb_set_d(temp2, -5.0/6.0);
    arb_mul(phi2, temp2, temp1, prec); // phi2 = -(5/6) * pi * alpha
    
    acb_abs(abs_z, z, prec);
    //arb_set_d(taylor_radius, TAYLOR_RADIUS_DBL);
    
    int in_wedge = open_wedge_p(z, phi1, phi2, prec);
    int magnitude_check = !in_region_G0_flint(z, alpha, taylor_radius, prec); //arb_gt(abs_z, taylor_radius);
    
    int result = in_wedge && magnitude_check;
    
    arb_clear(phi1);
    arb_clear(phi2);
    arb_clear(pi_val);
    arb_clear(abs_z);
    arb_clear(temp1);
    arb_clear(temp2);
    
    return result;
}

void
compute_r1_flint(arb_t result, const arb_t a, slong prec)
{
    arb_t acc, c0, pi, t1, t2, t3;
    
    /* Initialize temporary variables */
    arb_init(acc);
    arb_init(c0);
    arb_init(pi);
    arb_init(t1);
    arb_init(t2);
    arb_init(t3);
    
    /* Compute accuracy: acc = 2^-prec */
    // arb_set_d(acc, pow(2, -prec));
    accuracy_from_prec_flint(acc, prec);
    
    /* Compute pi */
    arb_const_pi(pi, prec);
    
    /* Compute denominator: pi * sin(pi * a) */
    arb_mul(t1, pi, a, prec);       /* t1 = pi * a */
    arb_sin(t2, t1, prec);          /* t2 = sin(pi * a) */
    arb_mul(t3, pi, t2, prec);      /* t3 = pi * sin(pi * a) */
    
    /* Compute numerator: 1.3^(1 - a) */
    arb_set_si(t1, 13);
    arb_div_si(t1, t1, 10, prec);   /* t1 = 1.3 */
    arb_sub_si(t2, a, 1, prec);
    arb_neg(t2, t2);                /* t2 = 1 - a */
    arb_pow(c0, t1, t2, prec);      /* c0 = 1.3^(1 - a) */
    
    /* Complete c0 calculation */
    arb_div(c0, c0, t3, prec);      /* c0 = numerator/denominator */
    
    /* Compute log(acc / c0) */
    arb_div(t1, acc, c0, prec);
    arb_log(t2, t1, prec);
    
    /* Compute final result */
    arb_mul_si(t3, t2, -2, prec);
    arb_pow(result, t3, a, prec);
    
    /* Clear temporary variables */
    arb_clear(acc);
    arb_clear(c0);
    arb_clear(pi);
    arb_clear(t1);
    arb_clear(t2);
    arb_clear(t3);
}


/* double */
/* compute_r1(double alpha, slong prec) */
/* { */

/*     arb_t alpha_arb, result_arb; */
/*     double result; */
    
/*     arb_init(alpha_arb); */
/*     arb_init(result_arb); */
/*     arb_set_d(alpha_arb, alpha); */
    
/*     compute_r1_flint(result_arb, alpha_arb, prec); */
/*     result = arf_get_d(arb_midref(result_arb), ARF_RND_NEAR); */
    
/*     arb_clear(alpha_arb); */
/*     arb_clear(result_arb); */
    
/*     return result; */

/* } */


/* /\* C wrapper functions for use with standard C types                    *\/ */
/* /\* These functions provide a simple interface for other languages       *\/ */

/* /\** */
/*  * C wrapper for region G0 checking */
/*  * @param z_real Real part of complex number z */
/*  * @param z_imag Imaginary part of complex number z */
/*  * @param alpha Parameter alpha */
/*  * @return 1 if z is in region G0, 0 otherwise */
/*  *\/ */
/* unsigned short int */
/* in_region_G0 (double complex z, double alpha, double r0, long int prec) */
/* { */
/*     arb_t alpha_arb, r0_arb; */
/*     acb_t z_acb; */
/*     int result; */
    
/*     arb_init(alpha_arb); */
/*     arb_init(r0_arb); */
/*     acb_init(z_acb); */
    
/*     arb_set_d(alpha_arb, alpha); */
/*     arb_set_d(r0_arb, r0); */
/*     acb_set_d_d(z_acb, creal(z), cimag(z)); */
    
/*     result = in_region_G0_flint(z_acb, alpha_arb, r0_arb, prec); */
    
/*     arb_clear(alpha_arb); */
/*     arb_clear(r0_arb); */
/*     acb_clear(z_acb); */
    
/*     return result; */

/* } */

/* /\** */
/*  * C wrapper for region G1 checking */
/*  * @param z_real Real part of complex number z */
/*  * @param z_imag Imaginary part of complex number z */
/*  * @param alpha Parameter alpha */
/*  * @return 1 if z is in region G1, 0 otherwise */
/*  *\/ */
/* unsigned short int */
/* in_region_G1 (double complex z, double alpha, double r0, slong prec) */
/* { */
    
/*     arb_t alpha_arb; */
/*     acb_t z_acb; */
/*     int result; */
    
/*     arb_init(alpha_arb); */
/*     acb_init(z_acb); */
    
/*     arb_set_d(alpha_arb, alpha); */
/*     acb_set_d_d(z_acb, creal(z), cimag(z)); */
    
/*     result = in_region_G1_flint(z_acb, alpha_arb, prec); */
    
/*     arb_clear(alpha_arb); */
/*     acb_clear(z_acb); */
    
/*     return result; */
       
/* } */

/* /\** */
/*  * C wrapper for region G2 checking */
/*  * @param z_real Real part of complex number z */
/*  * @param z_imag Imaginary part of complex number z */
/*  * @param alpha Parameter alpha */
/*  * @return 1 if z is in region G2, 0 otherwise */
/*  *\/ */
/* unsigned short int */
/* in_region_G2 (double complex z, double alpha, double r0, slong prec) */
/* { */
    
/*     arb_t alpha_arb; */
/*     acb_t z_acb; */
/*     int result; */
    
/*     arb_init(alpha_arb); */
/*     acb_init(z_acb); */
    
/*     arb_set_d(alpha_arb, alpha); */
/*     acb_set_d_d(z_acb, creal(z), cimag(z)); */
    
/*     result = in_region_G2_flint(z_acb, alpha_arb, prec); */
    
/*     arb_clear(alpha_arb); */
/*     acb_clear(z_acb); */
    
/*     return result; */
    
/* } */

/* /\** */
/*  * C wrapper for region G3 checking */
/*  * @param z_real Real part of complex number z */
/*  * @param z_imag Imaginary part of complex number z */
/*  * @param alpha Parameter alpha */
/*  * @return 1 if z is in region G3, 0 otherwise */
/*  *\/ */
/* unsigned short int */
/* in_region_G3(double complex z, double alpha, double r0, slong prec) */
/* { */
    
/*     arb_t alpha_arb; */
/*     acb_t z_acb; */
/*     int result; */
    
/*     arb_init(alpha_arb); */
/*     acb_init(z_acb); */
    
/*     arb_set_d(alpha_arb, alpha); */
/*     acb_set_d_d(z_acb, creal(z), cimag(z)); */
    
/*     result = in_region_G3_flint(z_acb, alpha_arb, prec); */
    
/*     arb_clear(alpha_arb); */
/*     acb_clear(z_acb); */
    
/*     return result; */
    
/* } */

/* /\** */
/*  * C wrapper for region G4 checking */
/*  * @param z_real Real part of complex number z */
/*  * @param z_imag Imaginary part of complex number z */
/*  * @param alpha Parameter alpha */
/*  * @return 1 if z is in region G4, 0 otherwise */
/*  *\/ */
/* unsigned short int */
/* in_region_G4(double complex z, double alpha, double r0, slong prec) */
/* { */
    
/*     arb_t alpha_arb; */
/*     acb_t z_acb; */
/*     int result; */
    
/*     arb_init(alpha_arb); */
/*     acb_init(z_acb); */
    
/*     arb_set_d(alpha_arb, alpha); */
/*     acb_set_d_d(z_acb, creal(z), cimag(z)); */
    
/*     result = in_region_G4_flint(z_acb, alpha_arb, prec); */
    
/*     arb_clear(alpha_arb); */
/*     acb_clear(z_acb); */
    
/*     return result; */
    
/* } */

/* /\** */
/*  * C wrapper for region G5 checking */
/*  * @param z_real Real part of complex number z */
/*  * @param z_imag Imaginary part of complex number z */
/*  * @param alpha Parameter alpha */
/*  * @return 1 if z is in region G5, 0 otherwise */
/*  *\/ */
/* unsigned short int */
/* in_region_G5(double complex z, double alpha, double r0, slong prec) */
/* { */
    
/*     arb_t alpha_arb, r0_arb; */
/*     acb_t z_acb; */
/*     int result; */
    
/*     arb_init(alpha_arb); */
/*     arb_init(r0_arb); */
/*     acb_init(z_acb); */
    
/*     arb_set_d(alpha_arb, alpha); */
/*     arb_set_d(r0_arb, r0); */
/*     acb_set_d_d(z_acb, creal(z), cimag(z)); */
    
/*     result = in_region_G5_flint(z_acb, alpha_arb, r0_arb, prec); */
    
/*     arb_clear(alpha_arb); */
/*     arb_clear(r0_arb); */
/*     acb_clear(z_acb); */
    
/*     return result; */

/* } */

/* /\** */
/*  * C wrapper for region G6 checking */
/*  * @param z_real Real part of complex number z */
/*  * @param z_imag Imaginary part of complex number z */
/*  * @param alpha Parameter alpha */
/*  * @return 1 if z is in region G6, 0 otherwise */
/*  *\/ */
/* unsigned short int */
/* in_region_G6(double complex z, double alpha, double r0, slong prec) */
/* { */
    
/*     int result; */
/*     arb_t alpha_arb, r0_arb; */
/*     acb_t z_acb; */
    
/*     arb_init(alpha_arb); */
/*     arb_init(r0_arb); */
/*     acb_init(z_acb); */
    
/*     arb_set_d(alpha_arb, alpha); */
/*     arb_set_d(r0_arb, r0); */
/*     acb_set_d_d(z_acb, creal(z), cimag(z)); */
    
/*     result = in_region_G6_flint(z_acb, alpha_arb, r0_arb, prec); */
    
/*     arb_clear(alpha_arb); */
/*     arb_clear(r0_arb); */
/*     acb_clear(z_acb); */
    
/*     return result; */

/* } */
