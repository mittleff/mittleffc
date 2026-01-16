/**
 * @file mittleff.c
 * @brief Main Mittag-Leffler function implementation
 *
 * This module implements the main entry point for computing the Mittag-Leffler
 * function E_{α,β}(z), following the algorithm described in:
 *
 * Garrappa, R. (2015). "Numerical evaluation of two and three parameter
 * Mittag-Leffler functions." SIAM Journal on Numerical Analysis, 53(3), 1350-1369.
 *
 * The implementation follows the Scheme reference implementation structure:
 * 1. Check if |z| <= taylor_radius → use Taylor series
 * 2. Check if α > 1 → apply recursion relation to reduce to α ≤ 1
 * 3. Determine region and apply appropriate method:
 *    - Asymptotic expansion for large |z|
 *    - Integral representation for intermediate |z|
 */

#include <config.h>
#include <flint/acb_hypgeom.h>
#include <mittleff/mittleff.h>
#include "../algorithm/algorithm.h"
#include "../partition/partition.h"
#include <flint/acb.h>
#include <flint/arb.h>
#include <math.h>

/* Mathematical constant π */
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286
#endif

/**
 * Auxiliary Mittag-Leffler computation
 *
 * This function assumes α ≤ 1 and dispatches to the appropriate method based on
 * the magnitude of z and its location in the complex plane.
 *
 * Based on mittleff.scm:35-41
 */
void ml_mittleff_aux(acb_t result, const acb_t alpha, const acb_t beta, const acb_t z,
                     slong prec, const arb_t taylor_radius)
{
    arb_t r1, z_abs;
    int region;

    arb_init(r1);
    arb_init(z_abs);

    /* Compute r1 threshold using partition module */
    extern void ml_compute_r1(arb_t r1, const acb_t alpha, slong prec);
    ml_compute_r1(r1, alpha, prec);

    /* Get |z| */
    acb_abs(z_abs, z, prec);

    /* Check if |z| >= r1 */
    if (arb_ge(z_abs, r1)) {
        /* Use asymptotic expansion (regions 1-4) */
        /* Create a non-const copy for ml_compute_region */
        arb_t radius_copy;
        arb_init(radius_copy);
        arb_set(radius_copy, taylor_radius);
        region = ml_compute_region(z, alpha, radius_copy, prec);
        arb_clear(radius_copy);

        /* Ensure we got an asymptotic region */
        if (region >= 1 && region <= 4) {
            ml_asymptotics(result, alpha, beta, z, region, prec);
        } else {
            /* Fallback to region 1 if detection failed */
            ml_asymptotics(result, alpha, beta, z, 1, prec);
        }
    } else {
        /* Use integral representation (regions 5-6) */
        /* Create a non-const copy for ml_compute_region */
        arb_t radius_copy;
        arb_init(radius_copy);
        arb_set(radius_copy, taylor_radius);
        region = ml_compute_region(z, alpha, radius_copy, prec);
        arb_clear(radius_copy);

        /* Ensure we got an integral region */
        if (region >= 5 && region <= 6) {
            ml_integral_rep(result, alpha, beta, z, region, prec);
        } else {
            /* Fallback to region 5 if detection failed */
            ml_integral_rep(result, alpha, beta, z, 5, prec);
        }
    }

    arb_clear(r1);
    arb_clear(z_abs);
}

/**
 * Main Mittag-Leffler function
 *
 * Computes E_{α,β}(z) using the optimal method based on parameters and argument.
 *
 * Algorithm:
 * 1. If |z| <= taylor_radius → Taylor series
 * 2. If α > 1 → Apply recursion relation:
 *    E_{α,β}(z) = (1/(2m+1)) * Σ_{l=-m}^{m} E_{α',β}(z')
 *    where:
 *      m = ceil((α - 1) / 2)
 *      α' = α / (2m + 1)
 *      z' = z^(1/(2m+1)) * exp(2πi*l/(2m+1))
 * 3. Otherwise (α ≤ 1) → Call mittleff_aux
 *
 * Based on mittleff.scm:11-33
 */
void ml_mittleff(acb_t result, const acb_t alpha, const acb_t beta, const acb_t z,
                 slong prec, const arb_t taylor_radius)
{
    arb_t z_abs, alpha_real, one, temp_arb;
    acb_t temp_acb;
    int use_taylor, use_recursion;

    arb_init(z_abs);
    arb_init(alpha_real);
    arb_init(one);
    arb_init(temp_arb);
    acb_init(temp_acb);

    arb_one(one);

    /* Get |z| and real part of α */
    acb_abs(z_abs, z, prec);
    acb_get_real(alpha_real, alpha);

    /* Check if we should use Taylor series: |z| <= taylor_radius */
    use_taylor = arb_le(z_abs, taylor_radius);

    if (use_taylor) {
      if (arb_is_zero(z_abs)) {
        /* fprintf(stderr, "|z| = 0\n"); */
        /* fflush(stderr); */
        acb_hypgeom_rgamma(result, beta, prec);
        
      }
      else {
        /* ====================================
         * Evaluate using Taylor Series (G0)
         * ==================================== */
        ml_taylor(result, alpha, beta, z, prec);
      }
    } else {
        /* Check if we should apply recursion: α > 1 */
        use_recursion = arb_gt(alpha_real, one);

        if (use_recursion) {
            /* =========================================
             * Apply Recursion Relation for α > 1
             * ========================================= */
            arb_t m_arb, two_m_plus_1, one_over_2mp1;
            acb_t a_prime, sum_result, term, z_power, exp_term, z_prime;
            slong m, l;

            arb_init(m_arb);
            arb_init(two_m_plus_1);
            arb_init(one_over_2mp1);
            acb_init(a_prime);
            acb_init(sum_result);
            acb_init(term);
            acb_init(z_power);
            acb_init(exp_term);
            acb_init(z_prime);

            /* Compute m = ceil((α - 1) / 2) */
            arb_sub_ui(m_arb, alpha_real, 1, prec);
            arb_div_ui(m_arb, m_arb, 2, prec);
            arb_ceil(m_arb, m_arb, prec);
            m = (slong)arf_get_d(arb_midref(m_arb), ARF_RND_NEAR);
            if (m < 1) m = 1;  /* Ensure m >= 1 */

            /* Compute 2m + 1 */
            arb_set_si(two_m_plus_1, 2 * m + 1);

            /* Compute 1 / (2m + 1) */
            arb_inv(one_over_2mp1, two_m_plus_1, prec);

            /* Compute α' = α / (2m + 1) */
            acb_set_arb(a_prime, one_over_2mp1);
            acb_mul(a_prime, alpha, a_prime, prec);

            /* Initialize sum */
            acb_zero(sum_result);

            /* Compute z^(1/(2m+1)) */
            acb_pow_arb(z_power, z, one_over_2mp1, prec);

            /* Sum from l = -m to m */
            for (l = -m; l <= m; l++) {
                arb_t two_pi_l_over_2mp1, angle;
                acb_t i_angle;

                arb_init(two_pi_l_over_2mp1);
                arb_init(angle);
                acb_init(i_angle);

                /* Compute angle = 2πl / (2m + 1) */
                arb_const_pi(angle, prec);
                arb_mul_2exp_si(angle, angle, 1);  /* 2π */
                arb_mul_si(angle, angle, l, prec);  /* 2πl */
                arb_div(angle, angle, two_m_plus_1, prec);  /* 2πl / (2m+1) */

                /* Compute exp(i * angle) */
                acb_set_arb(i_angle, angle);
                acb_mul_onei(i_angle, i_angle);  /* i * angle */
                acb_exp(exp_term, i_angle, prec);  /* exp(i * angle) */

                /* Compute z' = z^(1/(2m+1)) * exp(2πil/(2m+1)) */
                acb_mul(z_prime, z_power, exp_term, prec);

                /* Recursively compute E_{α',β}(z') */
                ml_mittleff_aux(term, a_prime, beta, z_prime, prec, taylor_radius);

                /* Add to sum */
                acb_add(sum_result, sum_result, term, prec);

                arb_clear(two_pi_l_over_2mp1);
                arb_clear(angle);
                acb_clear(i_angle);
            }

            /* Multiply sum by 1/(2m+1) */
            acb_mul_arb(result, sum_result, one_over_2mp1, prec);

            arb_clear(m_arb);
            arb_clear(two_m_plus_1);
            arb_clear(one_over_2mp1);
            acb_clear(a_prime);
            acb_clear(sum_result);
            acb_clear(term);
            acb_clear(z_power);
            acb_clear(exp_term);
            acb_clear(z_prime);
        } else {
            /* =============================================
             * Apply main algorithm for α ≤ 1
             * ============================================= */
            ml_mittleff_aux(result, alpha, beta, z, prec, taylor_radius);
        }
    }

    arb_clear(z_abs);
    arb_clear(alpha_real);
    arb_clear(one);
    arb_clear(temp_arb);
    acb_clear(temp_acb);
}
