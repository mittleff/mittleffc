#include <partition.h>


/**
 * Convert bit precision to decimal accuracy
 *
 * Mathematical formula: accuracy = 10^(-ceil(prec * log10(2)))
 *
 * This function computes the decimal accuracy threshold corresponding to
 * a given bit precision. It's used throughout the algorithm to determine
 * convergence criteria and threshold values.
 *
 * Implementation uses FLINT arbitrary precision arithmetic to ensure
 * accuracy even for high precision values.
 *
 * @param accuracy Output parameter (must be initialized by caller)
 * @param prec     Precision in bits
 *
 * Example:
 *   arb_t acc;
 *   arb_init(acc);
 *   ml_accuracy_from_prec(acc, 53);  // For double precision
 *   // acc now contains approximately 1.11e-16
 *   arb_clear(acc);
 */
void ml_accuracy_from_prec(arb_t accuracy, slong prec) {
    arb_t log10_2, prec_arb, exponent, base;

    /* Initialize temporary variables */
    arb_init(log10_2);
    arb_init(prec_arb);
    arb_init(exponent);
    arb_init(base);

    /* Compute log10(2) using FLINT's logarithm function */
    /* log10(2) = log(2) / log(10) */
    arb_const_log2(log10_2, prec);      /* log(2) */
    arb_t log10;
    arb_init(log10);
    arb_set_ui(log10, 10);
    arb_log(log10, log10, prec);        /* log(10) */
    arb_div(log10_2, log10_2, log10, prec);  /* log(2)/log(10) = log10(2) */
    arb_clear(log10);

    /* Convert precision to arb_t */
    arb_set_si(prec_arb, prec);

    /* Compute exponent = prec * log10(2) */
    arb_mul(exponent, prec_arb, log10_2, prec);

    /* Apply ceiling function */
    arb_ceil(exponent, exponent, prec);

    /* Negate to get -ceil(prec * log10(2)) */
    arb_neg(exponent, exponent);

    /* Compute accuracy = 10^exponent */
    arb_set_ui(base, 10);
    arb_pow(accuracy, base, exponent, prec);

    /* Cleanup */
    arb_clear(base);
    arb_clear(exponent);
    arb_clear(prec_arb);
    arb_clear(log10_2);
}

/**
 * Check if angle a is clockwise from angle b
 *
 * Uses the cross product (sine of angle difference) to determine orientation.
 * Matches Scheme implementation: sin(a - b) > 0 means clockwise
 *
 * @param a      First angle
 * @param b      Second angle
 * @param closed If true, use >= 0 instead of > 0
 * @return 1 if clockwise, 0 otherwise
 */
static int ml_clockwise(const arb_t a, const arb_t b, int closed) {
    arb_t diff, cross;
    int result;

    arb_init(diff);
    arb_init(cross);

    /* cross = sin(a - b) */
    arb_sub(diff, a, b, 53);
    arb_sin(cross, diff, 53);

    /* Check sign */
    if (closed) {
        result = !arb_is_negative(cross);  /* >= 0 */
    } else {
        result = arb_is_positive(cross);    /* > 0 */
    }

    arb_clear(diff);
    arb_clear(cross);

    return result;
}

/**
 * Check if angle c is between phi1 and phi2 in counterclockwise direction
 *
 * This function handles angle wrapping correctly using the clockwise predicate.
 *
 * Algorithm:
 * - If phi1 is clockwise from phi2, the region wraps around, so we check
 *   that c is NOT in the opposite region
 * - Otherwise, c must be counterclockwise from phi1 and phi2 must be
 *   counterclockwise from c
 *
 * @param c      Angle to test
 * @param phi1   Lower bound angle
 * @param phi2   Upper bound angle
 * @param closed If true, include boundaries
 * @return 1 if c is between phi1 and phi2, 0 otherwise
 */
static int ml_between(const arb_t c, const arb_t phi1, const arb_t phi2, int closed) {
    int result;

    /* If phi1 is clockwise from phi2, the region wraps around */
    if (ml_clockwise(phi1, phi2, 0)) {
        /* Check that c is NOT in the opposite region */
        result = !ml_between(c, phi2, phi1, closed);
    } else {
        /* c must be counterclockwise from phi1, and phi2 must be ccw from c */
        int cond1 = !ml_clockwise(phi1, c, closed);
        int cond2 = ml_clockwise(phi2, c, closed);
        result = cond1 && cond2;
    }

    return result;
}

/**
 * Check if angle c is in open wedge (phi1, phi2)
 * Open wedge: boundaries excluded
 */
static int ml_open_wedge(const arb_t c, const arb_t phi1, const arb_t phi2) {
    return ml_between(c, phi1, phi2, 0);  /* closed=0 → open */
}

/**
 * Check if angle c is in closed wedge [phi1, phi2]
 * Closed wedge: boundaries included
 */
static int ml_closed_wedge(const arb_t c, const arb_t phi1, const arb_t phi2) {
    return ml_between(c, phi1, phi2, 1);  /* closed=1 → closed */
}

/* ========================================================================== */
/* Threshold Radius Computation                                              */
/* ========================================================================== */

/**
 * Compute threshold radius r1 between asymptotic and integral methods
 *
 * Mathematical formula:
 *   r1 = (-2 * log(acc / c0))^alpha
 *
 * where:
 *   acc = 10^(-ceil(prec * log10(2)))      [from ml_accuracy_from_prec]
 *   c0 = 1.3^(1-alpha) / (π * sin(π*alpha))
 *
 * This threshold separates regions G1-G4 (asymptotic, |z| >= r1) from
 * regions G5-G6 (integral representation, |z| < r1).
 *
 * @param r1    Output parameter for threshold radius
 * @param alpha Parameter α (uses real part)
 * @param prec  Precision in bits
 */
void ml_compute_r1(arb_t r1, const acb_t alpha, slong prec) {
    arb_t acc, c0, alpha_real, term1, term2, log_term;
    arb_t pi, sin_val, power_val;

    /* Initialize temporary variables */
    arb_init(acc);
    arb_init(c0);
    arb_init(alpha_real);
    arb_init(term1);
    arb_init(term2);
    arb_init(log_term);
    arb_init(pi);
    arb_init(sin_val);
    arb_init(power_val);

    /* Compute accuracy from precision */
    ml_accuracy_from_prec(acc, prec);

    /* Extract real part of alpha */
    acb_get_real(alpha_real, alpha);

    /* Get π constant */
    arb_const_pi(pi, prec);

    /* Compute c0 = 1.3^(1-alpha) / (π * sin(π * alpha)) */

    /* Step 1: Compute 1.3^(1-alpha) */
    arb_set_d(term1, 1.3);
    arb_sub_ui(term2, alpha_real, 1, prec);  /* alpha - 1 */
    arb_neg(term2, term2);                    /* 1 - alpha */
    arb_pow(power_val, term1, term2, prec);   /* 1.3^(1-alpha) */

    /* Step 2: Compute π * sin(π * alpha) */
    arb_mul(term2, pi, alpha_real, prec);     /* π * alpha */
    arb_sin(sin_val, term2, prec);            /* sin(π * alpha) */
    arb_mul(term2, pi, sin_val, prec);        /* π * sin(π * alpha) */

    /* Step 3: Divide to get c0 */
    arb_div(c0, power_val, term2, prec);      /* c0 = numerator / denominator */

    /* Compute log_term = -2 * log(acc / c0) */
    arb_div(log_term, acc, c0, prec);         /* acc / c0 */
    arb_log(log_term, log_term, prec);        /* log(acc / c0) */
    arb_mul_2exp_si(log_term, log_term, 1);  /* Multiply by 2 (shift left) */
    arb_neg(log_term, log_term);              /* Negate to get -2 * log(...) */

    /* Compute r1 = log_term^alpha */
    arb_pow(r1, log_term, alpha_real, prec);

    /* Cleanup */
    arb_clear(power_val);
    arb_clear(sin_val);
    arb_clear(pi);
    arb_clear(log_term);
    arb_clear(term2);
    arb_clear(term1);
    arb_clear(alpha_real);
    arb_clear(c0);
    arb_clear(acc);
}

/* ========================================================================== */
/* Region Predicates                                                          */
/* ========================================================================== */

/**
 * Region G1: Wedge around negative real axis
 *
 * Parameters:
 *   delta = π * alpha * 0.125
 *   phi1 = -π * alpha + delta
 *   phi2 = +π * alpha - delta
 *
 * Condition: z in open-wedge(phi1, phi2)
 *
 * This region uses asymptotic expansion with specific error function corrections.
 */
int ml_in_region_G1(const acb_t z, const acb_t alpha) {
    arb_t alpha_real, delta, phi1, phi2, z_arg, pi;
    int result;

    arb_init(alpha_real);
    arb_init(delta);
    arb_init(phi1);
    arb_init(phi2);
    arb_init(z_arg);
    arb_init(pi);

    /* Extract real part of alpha */
    acb_get_real(alpha_real, alpha);

    /* Get π constant */
    arb_const_pi(pi, 53);

    /* Compute delta = π * alpha * 0.125 (= π * alpha / 8) */
    arb_mul(delta, pi, alpha_real, 53);
    arb_mul_2exp_si(delta, delta, -3);  /* Multiply by 1/8 */

    /* Compute phi1 = -π * alpha + delta */
    arb_mul(phi1, pi, alpha_real, 53);
    arb_neg(phi1, phi1);                /* -π * alpha */
    arb_add(phi1, phi1, delta, 53);     /* -π * alpha + delta */

    /* Compute phi2 = +π * alpha - delta */
    arb_mul(phi2, pi, alpha_real, 53);
    arb_sub(phi2, phi2, delta, 53);     /* +π * alpha - delta */

    /* Get angle of z */
    acb_arg(z_arg, z, 53);

    /* Check if z is in open wedge (phi1, phi2) */
    result = ml_open_wedge(z_arg, phi1, phi2);

    /* Cleanup */
    arb_clear(pi);
    arb_clear(z_arg);
    arb_clear(phi2);
    arb_clear(phi1);
    arb_clear(delta);
    arb_clear(alpha_real);

    return result;
}

/**
 * Region G2: Wedge below positive real axis
 *
 * Parameters:
 *   tau = min(π * alpha * 0.125, π * (1 + alpha) * 0.5)
 *   phi1 = +π * alpha + tau
 *   phi2 = -π * alpha - tau
 *
 * Condition: z in open-wedge(phi1, phi2)
 */
int ml_in_region_G2(const acb_t z, const acb_t alpha) {
    arb_t alpha_real, tau, tau1, tau2, phi1, phi2, z_arg, pi;
    arb_t temp;
    int result;

    arb_init(alpha_real);
    arb_init(tau);
    arb_init(tau1);
    arb_init(tau2);
    arb_init(phi1);
    arb_init(phi2);
    arb_init(z_arg);
    arb_init(pi);
    arb_init(temp);

    /* Extract real part of alpha */
    acb_get_real(alpha_real, alpha);

    /* Get π constant */
    arb_const_pi(pi, 53);

    /* Compute tau1 = π * alpha * 0.125 */
    arb_mul(tau1, pi, alpha_real, 53);
    arb_mul_2exp_si(tau1, tau1, -3);  /* Multiply by 1/8 */

    /* Compute tau2 = π * (1 + alpha) * 0.5 */
    arb_add_ui(temp, alpha_real, 1, 53);     /* 1 + alpha */
    arb_mul(tau2, pi, temp, 53);             /* π * (1 + alpha) */
    arb_mul_2exp_si(tau2, tau2, -1);         /* Multiply by 1/2 */

    /* tau = min(tau1, tau2) */
    if (arb_lt(tau1, tau2)) {
        arb_set(tau, tau1);
    } else {
        arb_set(tau, tau2);
    }

    /* Compute phi1 = +π * alpha + tau */
    arb_mul(phi1, pi, alpha_real, 53);
    arb_add(phi1, phi1, tau, 53);

    /* Compute phi2 = -π * alpha - tau */
    arb_mul(phi2, pi, alpha_real, 53);
    arb_neg(phi2, phi2);                /* -π * alpha */
    arb_sub(phi2, phi2, tau, 53);       /* -π * alpha - tau */

    /* Get angle of z */
    acb_arg(z_arg, z, 53);

    /* Check if z is in open wedge (phi1, phi2) */
    result = ml_open_wedge(z_arg, phi1, phi2);

    /* Cleanup */
    arb_clear(temp);
    arb_clear(pi);
    arb_clear(z_arg);
    arb_clear(phi2);
    arb_clear(phi1);
    arb_clear(tau2);
    arb_clear(tau1);
    arb_clear(tau);
    arb_clear(alpha_real);

    return result;
}

/**
 * Region G3: Wedge around positive imaginary axis
 *
 * Parameters:
 *   delta = π * alpha * 0.125
 *   tau = min(π * alpha * 0.125, π * (1 + alpha) * 0.5)
 *   phi1 = +π * alpha - delta
 *   phi2 = +π * alpha + tau
 *
 * Condition: z in closed-wedge[phi1, phi2]
 */
int ml_in_region_G3(const acb_t z, const acb_t alpha) {
    arb_t alpha_real, delta, tau, tau1, tau2, phi1, phi2, z_arg, pi;
    arb_t temp;
    int result;

    arb_init(alpha_real);
    arb_init(delta);
    arb_init(tau);
    arb_init(tau1);
    arb_init(tau2);
    arb_init(phi1);
    arb_init(phi2);
    arb_init(z_arg);
    arb_init(pi);
    arb_init(temp);

    /* Extract real part of alpha */
    acb_get_real(alpha_real, alpha);

    /* Get π constant */
    arb_const_pi(pi, 53);

    /* Compute delta = π * alpha * 0.125 */
    arb_mul(delta, pi, alpha_real, 53);
    arb_mul_2exp_si(delta, delta, -3);  /* Multiply by 1/8 */

    /* Compute tau1 = π * alpha * 0.125 */
    arb_set(tau1, delta);  /* Same as delta */

    /* Compute tau2 = π * (1 + alpha) * 0.5 */
    arb_add_ui(temp, alpha_real, 1, 53);     /* 1 + alpha */
    arb_mul(tau2, pi, temp, 53);             /* π * (1 + alpha) */
    arb_mul_2exp_si(tau2, tau2, -1);         /* Multiply by 1/2 */

    /* tau = min(tau1, tau2) */
    if (arb_lt(tau1, tau2)) {
        arb_set(tau, tau1);
    } else {
        arb_set(tau, tau2);
    }

    /* Compute phi1 = +π * alpha - delta */
    arb_mul(phi1, pi, alpha_real, 53);
    arb_sub(phi1, phi1, delta, 53);

    /* Compute phi2 = +π * alpha + tau */
    arb_mul(phi2, pi, alpha_real, 53);
    arb_add(phi2, phi2, tau, 53);

    /* Get angle of z */
    acb_arg(z_arg, z, 53);

    /* Check if z is in closed wedge [phi1, phi2] */
    result = ml_closed_wedge(z_arg, phi1, phi2);

    /* Cleanup */
    arb_clear(temp);
    arb_clear(pi);
    arb_clear(z_arg);
    arb_clear(phi2);
    arb_clear(phi1);
    arb_clear(tau2);
    arb_clear(tau1);
    arb_clear(tau);
    arb_clear(delta);
    arb_clear(alpha_real);

    return result;
}

/**
 * Region G4: Wedge above positive real axis
 *
 * Parameters:
 *   delta = π * alpha * 0.125
 *   tau = min(π * alpha * 0.125, π * (1 + alpha) * 0.5)
 *   phi1 = -π * alpha - tau
 *   phi2 = -π * alpha + delta
 *
 * Condition: z in closed-wedge[phi1, phi2]
 */
int ml_in_region_G4(const acb_t z, const acb_t alpha) {
    arb_t alpha_real, delta, tau, tau1, tau2, phi1, phi2, z_arg, pi;
    arb_t temp;
    int result;

    arb_init(alpha_real);
    arb_init(delta);
    arb_init(tau);
    arb_init(tau1);
    arb_init(tau2);
    arb_init(phi1);
    arb_init(phi2);
    arb_init(z_arg);
    arb_init(pi);
    arb_init(temp);

    /* Extract real part of alpha */
    acb_get_real(alpha_real, alpha);

    /* Get π constant */
    arb_const_pi(pi, 53);

    /* Compute delta = π * alpha * 0.125 */
    arb_mul(delta, pi, alpha_real, 53);
    arb_mul_2exp_si(delta, delta, -3);  /* Multiply by 1/8 */

    /* Compute tau1 = π * alpha * 0.125 */
    arb_set(tau1, delta);  /* Same as delta */

    /* Compute tau2 = π * (1 + alpha) * 0.5 */
    arb_add_ui(temp, alpha_real, 1, 53);     /* 1 + alpha */
    arb_mul(tau2, pi, temp, 53);             /* π * (1 + alpha) */
    arb_mul_2exp_si(tau2, tau2, -1);         /* Multiply by 1/2 */

    /* tau = min(tau1, tau2) */
    if (arb_lt(tau1, tau2)) {
        arb_set(tau, tau1);
    } else {
        arb_set(tau, tau2);
    }

    /* Compute phi1 = -π * alpha - tau */
    arb_mul(phi1, pi, alpha_real, 53);
    arb_neg(phi1, phi1);                /* -π * alpha */
    arb_sub(phi1, phi1, tau, 53);       /* -π * alpha - tau */

    /* Compute phi2 = -π * alpha + delta */
    arb_mul(phi2, pi, alpha_real, 53);
    arb_neg(phi2, phi2);                /* -π * alpha */
    arb_add(phi2, phi2, delta, 53);     /* -π * alpha + delta */

    /* Get angle of z */
    acb_arg(z_arg, z, 53);

    /* Check if z is in closed wedge [phi1, phi2] */
    result = ml_closed_wedge(z_arg, phi1, phi2);

    /* Cleanup */
    arb_clear(temp);
    arb_clear(pi);
    arb_clear(z_arg);
    arb_clear(phi2);
    arb_clear(phi1);
    arb_clear(tau2);
    arb_clear(tau1);
    arb_clear(tau);
    arb_clear(delta);
    arb_clear(alpha_real);

    return result;
}

/**
 * Region G5: Intermediate lower region
 *
 * Parameters:
 *   phi1 = -5π * alpha / 6
 *   phi2 = +5π * alpha / 6
 *
 * Condition: z in closed-wedge[phi1, phi2] AND |z| >= taylor_radius
 *
 * This region uses integral representation.
 */
int ml_in_region_G5(const acb_t z, const acb_t alpha, arb_t taylor_radius) {
    arb_t alpha_real, phi1, phi2, z_arg, z_abs, pi;
    arb_t coeff;
    int result;

    arb_init(alpha_real);
    arb_init(phi1);
    arb_init(phi2);
    arb_init(z_arg);
    arb_init(z_abs);
    arb_init(pi);
    arb_init(coeff);

    /* Extract real part of alpha */
    acb_get_real(alpha_real, alpha);

    /* Get π constant */
    arb_const_pi(pi, 53);

    /* Compute coefficient: 5π/6 */
    arb_mul_ui(coeff, pi, 5, 53);       /* 5π */
    arb_div_ui(coeff, coeff, 6, 53);    /* 5π/6 */

    /* Compute phi1 = -5π * alpha / 6 */
    arb_mul(phi1, coeff, alpha_real, 53);
    arb_neg(phi1, phi1);

    /* Compute phi2 = +5π * alpha / 6 */
    arb_mul(phi2, coeff, alpha_real, 53);

    /* Get angle and magnitude of z */
    acb_arg(z_arg, z, 53);
    acb_abs(z_abs, z, 53);

    /* Check conditions: in wedge AND |z| >= taylor_radius */
    int wedge_result = ml_closed_wedge(z_arg, phi1, phi2);
    double z_mag_val = arf_get_d(arb_midref(z_abs), ARF_RND_NEAR);
    int mag_check = (arb_ge(z_abs, taylor_radius)); // z_mag_val >= taylor_radius);

    if (wedge_result && mag_check) {
        result = 1;
    } else {
        result = 0;
    }

    /* Cleanup */
    arb_clear(coeff);
    arb_clear(pi);
    arb_clear(z_abs);
    arb_clear(z_arg);
    arb_clear(phi2);
    arb_clear(phi1);
    arb_clear(alpha_real);

    return result;
}

/**
 * Region G6: Intermediate upper region
 *
 * Parameters:
 *   phi1 = +5π * alpha / 6
 *   phi2 = -5π * alpha / 6
 *
 * Condition: z in open-wedge(phi1, phi2) AND |z| >= taylor_radius
 *
 * This region uses integral representation.
 */
      int ml_in_region_G6(const acb_t z, const acb_t alpha, arb_t taylor_radius) {
    arb_t alpha_real, phi1, phi2, z_arg, z_abs, pi;
    arb_t coeff;
    int result;

    arb_init(alpha_real);
    arb_init(phi1);
    arb_init(phi2);
    arb_init(z_arg);
    arb_init(z_abs);
    arb_init(pi);
    arb_init(coeff);

    /* Extract real part of alpha */
    acb_get_real(alpha_real, alpha);

    /* Get π constant */
    arb_const_pi(pi, 53);

    /* Compute coefficient: 5π/6 */
    arb_mul_ui(coeff, pi, 5, 53);       /* 5π */
    arb_div_ui(coeff, coeff, 6, 53);    /* 5π/6 */

    /* Compute phi1 = +5π * alpha / 6 */
    arb_mul(phi1, coeff, alpha_real, 53);

    /* Compute phi2 = -5π * alpha / 6 */
    arb_mul(phi2, coeff, alpha_real, 53);
    arb_neg(phi2, phi2);

    /* Get angle and magnitude of z */
    acb_arg(z_arg, z, 53);
    acb_abs(z_abs, z, 53);

    /* Check conditions: in wedge AND |z| >= taylor_radius */
    if (ml_open_wedge(z_arg, phi1, phi2) &&
        // arf_get_d(arb_midref(z_abs), ARF_RND_NEAR) >= taylor_radius
        arb_ge(z_abs, taylor_radius)
        ) {
        result = 1;
    } else {
        result = 0;
    }

    /* Cleanup */
    arb_clear(coeff);
    arb_clear(pi);
    arb_clear(z_abs);
    arb_clear(z_arg);
    arb_clear(phi2);
    arb_clear(phi1);
    arb_clear(alpha_real);

    return result;
}

/**
 * Determine which region z belongs to
 *
 * Algorithm:
 * 1. If |z| <= taylor_radius -> Region G0 (for Taylor series)
 * 2. Compute threshold r1
 * 3. If |z| >= r1 -> Check regions G1-G4 (for asymptotic expansions)
 * 4. Else -> Check regions G5-G6 (for integral representations)
 *
 * @param z             Complex argument
 * @param alpha         Parameter α
 * @param taylor_radius Threshold for Region G0
 * @param prec          Precision in bits
 * @return Region number (0-6), or -1 if no region matched (error)
 */
int ml_compute_region(const acb_t z, const acb_t alpha, arb_t taylor_radius,
                      slong prec) {

  arb_t z_abs, r1;
  double z_mag;
    
  arb_init(r1);
  arb_init(z_abs);

  /* Compute |z| */
  acb_abs(z_abs, z, prec);

  /* Region G0: |z| <= taylor_radius */
  if (arb_le(z_abs, taylor_radius)) {
    arb_clear(r1);
    arb_clear(z_abs);
    return 0;
  } else {
    /* Compute threshold r1 */
    ml_compute_r1(r1, alpha, prec);

    /* Check if |z| >= r1 (asymptotic regions) or |z| < r1 (integral regions) */
    if (arb_ge(z_abs, r1)) {
      /* Asymptotic regions (G1-G4) */
      if (ml_in_region_G1(z, alpha)) {
        arb_clear(r1);
        arb_clear(z_abs);
        return 1;
      }
      if (ml_in_region_G2(z, alpha)) {
        arb_clear(r1);
        arb_clear(z_abs);
        return 2;
      }
      if (ml_in_region_G3(z, alpha)) {
        arb_clear(r1);
        arb_clear(z_abs);
        return 3;
      }
      if (ml_in_region_G4(z, alpha)) {
        arb_clear(r1);
        arb_clear(z_abs);
        return 4;
      }
    } else { /* Integral regions (G5-G6) */
      if (ml_in_region_G5(z, alpha, taylor_radius)) {
        arb_clear(r1);
        arb_clear(z_abs);
        return 5;
      }
      if (ml_in_region_G6(z, alpha, taylor_radius)) {
        arb_clear(r1);
        arb_clear(z_abs);
        return 6;
      }
    }
  }

  /* Error: no region matched */
  arb_clear(r1);
  arb_clear(z_abs);
  return -1;
}
