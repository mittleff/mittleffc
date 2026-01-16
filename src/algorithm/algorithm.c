#include "algorithm.h"
#include "../quad/quad.h"
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>
#include <flint/acb_calc.h>
#include <flint/mag.h>
#include <math.h>

/* Constants matching the Scheme implementation */
#define ML_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286
/* 10000 Match Scheme implementation */
#define MITTLEFF_QUAD_INTEGRATION_MAX_DEPTH 10000
/* Integration tolerance */
#define MITTLEFF_QUAD_INTEGRATION_TOLERANCE 1e-10

/**
 * Helper: Compute accuracy from precision
 * accuracy = 10^(-ceil(|log10(2^(-prec))|))
 *          = 10^(-ceil(prec * log10(2)))
 */
static void accuracy_from_prec(arb_t acc, slong prec)
{
    arb_t tmp, log10_2;
    arb_init(tmp);
    arb_init(log10_2);

    /* log10(2) ≈ 0.30102999566 */
    arb_const_log10(log10_2, prec);
    arb_set_ui(tmp, 2);
    arb_log(tmp, tmp, prec);
    arb_div(log10_2, tmp, log10_2, prec);

    /* prec * log10(2) */
    arb_mul_si(tmp, log10_2, prec, prec);
    arb_ceil(tmp, tmp, prec);
    arb_neg(tmp, tmp);

    /* 10^(-ceil(prec * log10(2))) */
    arb_set_ui(acc, 10);
    arb_pow(acc, acc, tmp, prec);

    arb_clear(tmp);
    arb_clear(log10_2);
}

/**
 * Core series computation
 *
 * Computes: Σ_{k=kmin}^{kmax} z^(s*k) / Γ(a*s*k + b)
 *
 * Based on algorithm.scm:16-37
 */
void ml_series(acb_t result, const acb_t a, const acb_t b, const acb_t z,
               int asymptotic, slong prec)
{
    arb_t acc, abs_z, tmp_arb, one_minus_b_over_a;
    acb_t term, z_power, gamma_arg, z_to_1_over_a;
    slong kmin, kmax, k, s;

    arb_init(acc);
    arb_init(abs_z);
    arb_init(tmp_arb);
    arb_init(one_minus_b_over_a);
    acb_init(term);
    acb_init(z_power);
    acb_init(gamma_arg);
    acb_init(z_to_1_over_a);

    accuracy_from_prec(acc, prec);
    acb_abs(abs_z, z, prec);

    /* Determine kmin, kmax, and s based on asymptotic flag */
    if (asymptotic) {
        s = -1;
        kmin = 1;

        /* kmax = 1 + ceil((1/a) * |z|^(1/a)) */
        acb_inv(z_to_1_over_a, a, prec);
        acb_pow(z_to_1_over_a, z, z_to_1_over_a, prec);
        acb_abs(tmp_arb, z_to_1_over_a, prec);
        arb_ceil(tmp_arb, tmp_arb, prec);
        arb_add_ui(tmp_arb, tmp_arb, 1, prec);
        kmax = (slong)arf_get_d(arb_midref(tmp_arb), ARF_RND_NEAR);
    } else {
        s = 1;
        kmin = 0;

        /* Compute kmax as max of k1 and k2 */

        /* k1 = 1 + ceil(|(1-b)/a|) */
        acb_sub_ui(term, b, 1, prec);
        acb_neg(term, term);
        acb_div(term, term, a, prec);
        acb_abs(tmp_arb, term, prec);
        arb_ceil(tmp_arb, tmp_arb, prec);
        arb_add_ui(tmp_arb, tmp_arb, 1, prec);
        slong k1 = (slong)arf_get_d(arb_midref(tmp_arb), ARF_RND_NEAR);

        /* k2 = 1 + ceil(|log(acc * (|z| - 1)) / log(|z|)|) if |z| != 0 */
        slong k2 = 1;
        if (arb_is_nonzero(abs_z)) { // !arb_is_zero(abs_z)) {
        
            arb_t log_z, log_num, z_minus_1;
            arb_init(log_z);
            arb_init(log_num);
            arb_init(z_minus_1);

            /* fprintf(stdout, "|z| = "); */
            /* arb_fprintd(stdout, abs_z, 6); */
            /* fprintf(stdout, "\n"); */
            /* fflush(stdout); */

            

            arb_sub_ui(z_minus_1, abs_z, 1, prec);
            /* fprintf(stdout, "|z| - 1 = "); */
            /* arb_fprintd(stdout, z_minus_1, 6); */
            /* fprintf(stdout, "\n"); */
            /* fflush(stdout); */

            /* fprintf(stdout, "acc = "); */
            /* arb_fprintd(stdout, acc, 6); */
            /* fprintf(stdout, "\n"); */
            /* fflush(stdout); */

            arb_mul(log_num, acc, z_minus_1, prec);
            /* fprintf(stdout, "acc*(|z| - 1) = "); */
            /* arb_fprintd(stdout, log_num, 6); */
            /* fprintf(stdout, "\n"); */
            /* fflush(stdout); */

            acb_t log_num_acb;
            acb_init(log_num_acb);
            
            acb_set_arb(log_num_acb, log_num);
            acb_log(log_num_acb, log_num_acb, prec);
            // arb_log(log_num, log_num, prec);

            /* fprintf(stdout, "log(acc*(|z| - 1)) = "); */
            /* acb_fprintd(stdout, log_num_acb, 6); */
            /* fprintf(stdout, "\n"); */
            /* fflush(stdout); */

            arb_log(log_z, abs_z, prec);
            /* arb_fprintd(stdout, log_num, 6); */
            /* fprintf(stdout, "\n"); */
            /* arb_fprintd(stdout, log_z, 6); */
            /* fprintf(stdout, "\n"); */
            /* fflush(stdout); */

            if (arb_is_nonzero(log_z)) { // (!arb_is_zero(log_z)) {
              /* fprintf(stdout, "hellow\n"); */
              /* fflush(stdout); */

              acb_div_arb(log_num_acb, log_num_acb, log_z, prec);
              
              // arb_div(tmp_arb, log_num, log_z, prec);
              //acb_set_arb(term, tmp_arb);
                acb_abs(tmp_arb, log_num_acb, prec);
                arb_ceil(tmp_arb, tmp_arb, prec);
                arb_add_ui(tmp_arb, tmp_arb, 1, prec);
                k2 = (slong)arf_get_d(arb_midref(tmp_arb), ARF_RND_NEAR);
            }

            /* fprintf(stdout, "kk2= %d\n", k2); */
            /* fflush(stdout); */

            arb_clear(log_z);
            arb_clear(log_num);
            arb_clear(z_minus_1);
            acb_clear(log_num_acb);
        }

        kmax = (k1 > k2) ? k1 : k2;
        /* fprintf(stdout, "k1, k2, kmax = %d, %d, %d", k1, k2, kmax); */
        /* fflush(stdout); */
    }
    

    /* Compute the series sum */
    acb_zero(result);

    /* fprintf(stdout, "kmin, kmax = %d, %d\n", kmin, kmax); */
    /* fflush(stdout); */
    for (k = kmin; k <= kmax; k++) {
      /* fprintf(stdout, "k = %d\n", k); */
      /* fflush(stdout); */
        /* Compute z^(s*k) */
        if (s * k == 0) {
            acb_one(z_power);
        } else {
            acb_pow_si(z_power, z, s * k, prec);
        }

        /* Check if z^(s*k) is zero */
        if (acb_is_zero(z_power)) {
            continue;
        }

        /* Compute Γ(a*s*k + b) using 1/Γ(x) */
        acb_mul_si(gamma_arg, a, s * k, prec);
        acb_add(gamma_arg, gamma_arg, b, prec);

        /* Compute 1/Γ(gamma_arg) */
        acb_rgamma(term, gamma_arg, prec);

        /* Check if 1/Γ is zero */
        if (acb_is_zero(term)) {
            continue;
        }

        /* term = z^(s*k) * (1/Γ(a*s*k + b)) */
        acb_mul(term, z_power, term, prec);

        /* Add to result */
        acb_add(result, result, term, prec);
    }

    arb_clear(acc);
    arb_clear(abs_z);
    arb_clear(tmp_arb);
    arb_clear(one_minus_b_over_a);
    acb_clear(term);
    acb_clear(z_power);
    acb_clear(gamma_arg);
    acb_clear(z_to_1_over_a);
}

/**
 * Taylor series expansion
 *
 * Based on algorithm.scm:39-40
 */
void ml_taylor(acb_t result, const acb_t a, const acb_t b, const acb_t z, slong prec)
{
    ml_series(result, a, b, z, 0, prec);
}

/**
 * Helper function: c(theta) for asymptotic expansion
 *
 * c(θ) = sqrt(2 * (2*sin²(θ/2) + i*(θ - sin(θ))))
 *
 * Based on algorithm.scm:51-52
 */
static void compute_c_theta(acb_t result, const arb_t theta, slong prec)
{
    arb_t theta_half, sin_theta_half, sin_theta;
    arb_t real_part, imag_part, two;
    acb_t inside;

    arb_init(theta_half);
    arb_init(sin_theta_half);
    arb_init(sin_theta);
    arb_init(real_part);
    arb_init(imag_part);
    arb_init(two);
    acb_init(inside);

    arb_set_ui(two, 2);

    /* theta/2 */
    arb_div_ui(theta_half, theta, 2, prec);

    /* sin(theta/2) */
    arb_sin(sin_theta_half, theta_half, prec);

    /* sin(theta) */
    arb_sin(sin_theta, theta, prec);

    /* Real part: 2 * sin²(theta/2) */
    arb_sqr(real_part, sin_theta_half, prec);
    arb_mul(real_part, real_part, two, prec);

    /* Imaginary part: theta - sin(theta) */
    arb_sub(imag_part, theta, sin_theta, prec);

    /* inside = 2 * (real_part + i*imag_part) */
    acb_set_arb_arb(inside, real_part, imag_part);
    acb_mul_arb(inside, inside, two, prec);

    /* result = sqrt(inside) */
    acb_sqrt(result, inside, prec);

    arb_clear(theta_half);
    arb_clear(sin_theta_half);
    arb_clear(sin_theta);
    arb_clear(real_part);
    arb_clear(imag_part);
    arb_clear(two);
    acb_clear(inside);
}

/**
 * Helper function: ρ(i) for asymptotic expansion
 *
 * Based on algorithm.scm:53-59
 */
static void compute_rho(acb_t result, int region, const acb_t z, const acb_t a, slong prec)
{
    arb_t theta, abs_z_to_1_over_a, sqrt_half, tmp;
    acb_t c_theta, arg_erfc, z_to_1_over_a;

    if (region == 1) {
        acb_one(result);
        return;
    } else if (region == 2) {
        acb_zero(result);
        return;
    }

    /* For regions 3 and 4, we need erfc */
    arb_init(theta);
    arb_init(abs_z_to_1_over_a);
    arb_init(sqrt_half);
    arb_init(tmp);
    acb_init(c_theta);
    acb_init(arg_erfc);
    acb_init(z_to_1_over_a);

    /* theta = arg(z) */
    acb_arg(theta, z, prec);

    /* c(theta) */
    compute_c_theta(c_theta, theta, prec);

    /* |z|^(1/a) */
    acb_inv(z_to_1_over_a, a, prec);
    acb_pow(z_to_1_over_a, z, z_to_1_over_a, prec);
    acb_abs(abs_z_to_1_over_a, z_to_1_over_a, prec);

    /* sqrt(0.5) */
    arb_set_d(sqrt_half, 0.5);
    arb_sqrt(sqrt_half, sqrt_half, prec);

    /* sqrt(0.5 * |z|^(1/a)) */
    arb_mul(tmp, sqrt_half, abs_z_to_1_over_a, prec);
    arb_sqrt(tmp, tmp, prec);
    acb_set_arb(arg_erfc, tmp);

    /* c(theta) * sqrt(0.5 * |z|^(1/a)) */
    acb_mul(arg_erfc, c_theta, arg_erfc, prec);

    if (region == 3) {
        /* rho = 0.5 * erfc(+1 * arg_erfc) */
        acb_hypgeom_erfc(result, arg_erfc, prec);
    } else { /* region == 4 */
        /* rho = 0.5 * erfc(-1 * arg_erfc) */
        acb_neg(arg_erfc, arg_erfc);
        acb_hypgeom_erfc(result, arg_erfc, prec);
    }

    acb_mul_2exp_si(result, result, -1); /* multiply by 0.5 */

    arb_clear(theta);
    arb_clear(abs_z_to_1_over_a);
    arb_clear(sqrt_half);
    arb_clear(tmp);
    acb_clear(c_theta);
    acb_clear(arg_erfc);
    acb_clear(z_to_1_over_a);
}

/**
 * Asymptotic expansion
 *
 * Based on algorithm.scm:42-63
 */
void ml_asymptotics(acb_t result, const acb_t a, const acb_t b, const acb_t z,
                    int region, slong prec)
{
    acb_t rho_val, fac1, fac2, fac3, series_term, z_to_1_over_a, one_minus_b_over_a;

    acb_init(rho_val);
    acb_init(fac1);
    acb_init(fac2);
    acb_init(fac3);
    acb_init(series_term);
    acb_init(z_to_1_over_a);
    acb_init(one_minus_b_over_a);

    /* Compute ρ(region) */
    compute_rho(rho_val, region, z, a, prec);

    /* fac1 = ρ(region) / a */
    acb_div(fac1, rho_val, a, prec);

    /* fac2 = z^((1-b)/a) */
    acb_sub_ui(one_minus_b_over_a, b, 1, prec);
    acb_neg(one_minus_b_over_a, one_minus_b_over_a);
    acb_div(one_minus_b_over_a, one_minus_b_over_a, a, prec);
    acb_pow(fac2, z, one_minus_b_over_a, prec);

    /* fac3 = exp(z^(1/a)) */
    acb_inv(z_to_1_over_a, a, prec);
    acb_pow(z_to_1_over_a, z, z_to_1_over_a, prec);
    acb_exp(fac3, z_to_1_over_a, prec);

    /* result = fac1 * fac2 * fac3 */
    acb_mul(result, fac1, fac2, prec);
    acb_mul(result, result, fac3, prec);

    /* Compute asymptotic series correction */
    ml_series(series_term, a, b, z, 1, prec);

    /* result = fac1 * fac2 * fac3 - series_term */
    acb_sub(result, result, series_term, prec);

    acb_clear(rho_val);
    acb_clear(fac1);
    acb_clear(fac2);
    acb_clear(fac3);
    acb_clear(series_term);
    acb_clear(z_to_1_over_a);
    acb_clear(one_minus_b_over_a);
}

/**
 * Helper function: fn_a(r, a, b, x) for integral representation
 *
 * fn_a = (1/a) * r^((1-b)/a) * exp(r^(1/a) * cos(x/a))
 *
 * Based on algorithm.scm:65-66
 */
static void compute_fn_a(acb_t result, const acb_t r, const acb_t a, const acb_t b,
                         const arb_t x, slong prec)
{
    acb_t r_to_power, exp_term, one_over_a, one_minus_b_over_a, x_over_a, cos_term;

    acb_init(r_to_power);
    acb_init(exp_term);
    acb_init(one_over_a);
    acb_init(one_minus_b_over_a);
    acb_init(x_over_a);
    acb_init(cos_term);

    /* 1/a */
    acb_inv(one_over_a, a, prec);

    /* (1-b)/a */
    acb_sub_ui(one_minus_b_over_a, b, 1, prec);
    acb_neg(one_minus_b_over_a, one_minus_b_over_a);
    acb_div(one_minus_b_over_a, one_minus_b_over_a, a, prec);

    /* r^((1-b)/a) */
    /* TODO: When r=0 and b=1, this computes 0^0 which FLINT marks as a special value,
     * causing integration to fail with "ulp error not defined for special value".
     * This occurs when integrating from lower bound 0 in ml_integral_rep when beta <= 1. */
    acb_pow(r_to_power, r, one_minus_b_over_a, prec);

    /* x/a */
    acb_set_arb(x_over_a, x);
    acb_div(x_over_a, x_over_a, a, prec);

    /* cos(x/a) */
    acb_cos(cos_term, x_over_a, prec);

    /* r^(1/a) */
    acb_pow(exp_term, r, one_over_a, prec);

    /* r^(1/a) * cos(x/a) */
    acb_mul(exp_term, exp_term, cos_term, prec);

    /* exp(r^(1/a) * cos(x/a)) */
    acb_exp(exp_term, exp_term, prec);

    /* result = (1/a) * r^((1-b)/a) * exp_term */
    acb_mul(result, one_over_a, r_to_power, prec);
    acb_mul(result, result, exp_term, prec);

    acb_clear(r_to_power);
    acb_clear(exp_term);
    acb_clear(one_over_a);
    acb_clear(one_minus_b_over_a);
    acb_clear(x_over_a);
    acb_clear(cos_term);
}

/**
 * Helper function: ω(x, y, a, b) for integral representation
 *
 * ω = x^(1/a) * sin(y/a) + y * (1 + (1-b)/a)
 *
 * Based on algorithm.scm:68-69
 */
static void compute_omega(arb_t result, const arb_t x, const arb_t y,
                          const acb_t a, const acb_t b, slong prec)
{
    arb_t x_to_1_over_a, sin_term, y_factor, term1, term2, one_over_a, one_minus_b_over_a;
    acb_t a_inv, power_val;

    arb_init(x_to_1_over_a);
    arb_init(sin_term);
    arb_init(y_factor);
    arb_init(term1);
    arb_init(term2);
    arb_init(one_over_a);
    arb_init(one_minus_b_over_a);
    acb_init(a_inv);
    acb_init(power_val);

    /* 1/a */
    acb_inv(a_inv, a, prec);
    acb_get_real(one_over_a, a_inv);

    /* x^(1/a) */
    /* TODO: When x=0, this could produce special values that cause integration issues */
    arb_pow(x_to_1_over_a, x, one_over_a, prec);

    /* y/a */
    arb_div(sin_term, y, acb_realref(a), prec);

    /* sin(y/a) */
    arb_sin(sin_term, sin_term, prec);

    /* term1 = x^(1/a) * sin(y/a) */
    arb_mul(term1, x_to_1_over_a, sin_term, prec);

    /* (1-b)/a */
    acb_sub_ui(power_val, b, 1, prec);
    acb_neg(power_val, power_val);
    acb_div(power_val, power_val, a, prec);
    acb_get_real(one_minus_b_over_a, power_val);

    /* 1 + (1-b)/a */
    arb_add_ui(y_factor, one_minus_b_over_a, 1, prec);

    /* term2 = y * (1 + (1-b)/a) */
    arb_mul(term2, y, y_factor, prec);

    /* result = term1 + term2 */
    arb_add(result, term1, term2, prec);

    arb_clear(x_to_1_over_a);
    arb_clear(sin_term);
    arb_clear(y_factor);
    arb_clear(term1);
    arb_clear(term2);
    arb_clear(one_over_a);
    arb_clear(one_minus_b_over_a);
    acb_clear(a_inv);
    acb_clear(power_val);
}

/**
 * Helper function: fn_b(r, a, b, z, phi) for integral representation
 *
 * Based on algorithm.scm:71-74
 */
static void compute_fn_b(acb_t result, const arb_t r_val, const acb_t a, const acb_t b,
                         const acb_t z, const arb_t phi, slong prec)
{
    arb_t omega_val, w_minus_phi, sin_w_minus_phi, sin_w;
    acb_t fn_a_val, r, numerator, denominator, pi_val, r_squared, z_squared, temp;

    arb_init(omega_val);
    arb_init(w_minus_phi);
    arb_init(sin_w_minus_phi);
    arb_init(sin_w);
    acb_init(fn_a_val);
    acb_init(r);
    acb_init(numerator);
    acb_init(denominator);
    acb_init(pi_val);
    acb_init(r_squared);
    acb_init(z_squared);
    acb_init(temp);

    acb_set_arb(r, r_val);

    /* ω(r, φ, a, b) */
    compute_omega(omega_val, r_val, phi, a, b, prec);

    /* fn_a(r, a, b, φ) */
    compute_fn_a(fn_a_val, r, a, b, phi, prec);

    /* w - φ */
    arb_sub(w_minus_phi, omega_val, phi, prec);

    /* sin(w - φ) */
    arb_sin(sin_w_minus_phi, w_minus_phi, prec);

    /* sin(w) */
    arb_sin(sin_w, omega_val, prec);

    /* numerator = r * sin(w - φ) - z * sin(w) */
    acb_mul_arb(numerator, r, sin_w_minus_phi, prec);
    acb_mul_arb(temp, z, sin_w, prec);
    acb_sub(numerator, numerator, temp, prec);

    /* denominator = r² - 2rz*cos(φ) + z² */
    acb_sqr(r_squared, r, prec);
    acb_sqr(z_squared, z, prec);
    arb_cos(w_minus_phi, phi, prec); /* reuse w_minus_phi for cos(φ) */
    acb_mul_arb(temp, z, w_minus_phi, prec);
    acb_mul(temp, r, temp, prec);
    acb_mul_2exp_si(temp, temp, 1); /* multiply by 2 */
    acb_add(denominator, r_squared, z_squared, prec);
    acb_sub(denominator, denominator, temp, prec);

    /* result = (1/π) * fn_a * (numerator / denominator) */
    acb_const_pi(pi_val, prec);
    acb_inv(pi_val, pi_val, prec);
    acb_div(result, numerator, denominator, prec);
    acb_mul(result, result, fn_a_val, prec);
    acb_mul(result, result, pi_val, prec);

    arb_clear(omega_val);
    arb_clear(w_minus_phi);
    arb_clear(sin_w_minus_phi);
    arb_clear(sin_w);
    acb_clear(fn_a_val);
    acb_clear(r);
    acb_clear(numerator);
    acb_clear(denominator);
    acb_clear(pi_val);
    acb_clear(r_squared);
    acb_clear(z_squared);
    acb_clear(temp);
}

/**
 * Helper function: fn_c(ph, a, b, z, rho) for integral representation
 *
 * Based on algorithm.scm:76-78
 */
static void compute_fn_c(acb_t result, const arb_t ph, const acb_t a, const acb_t b,
                         const acb_t z, const arb_t rho_val, slong prec)
{
    arb_t omega_val;
    acb_t fn_a_val, rho, two_pi, i_omega, exp_i_omega, i_ph, exp_i_ph, denominator;

    arb_init(omega_val);
    acb_init(fn_a_val);
    acb_init(rho);
    acb_init(two_pi);
    acb_init(i_omega);
    acb_init(exp_i_omega);
    acb_init(i_ph);
    acb_init(exp_i_ph);
    acb_init(denominator);

    acb_set_arb(rho, rho_val);

    /* ω(ρ, ph, a, b) */
    compute_omega(omega_val, rho_val, ph, a, b, prec);

    /* fn_a(ρ, a, b, ph) */
    compute_fn_a(fn_a_val, rho, a, b, ph, prec);

    /* 2π */
    acb_const_pi(two_pi, prec);
    acb_mul_2exp_si(two_pi, two_pi, 1);

    /* i*ω */
    acb_set_arb(i_omega, omega_val);
    acb_mul_onei(i_omega, i_omega);

    /* exp(i*ω) */
    acb_exp(exp_i_omega, i_omega, prec);

    /* i*ph */
    acb_set_arb(i_ph, ph);
    acb_mul_onei(i_ph, i_ph);

    /* ρ * exp(i*ph) */
    acb_exp(exp_i_ph, i_ph, prec);
    acb_mul(exp_i_ph, rho, exp_i_ph, prec);

    /* denominator = ρ * exp(i*ph) - z */
    acb_sub(denominator, exp_i_ph, z, prec);

    /* result = (ρ / (2π)) * fn_a * (exp(i*ω) / denominator) */
    acb_div(result, exp_i_omega, denominator, prec);
    acb_mul(result, result, fn_a_val, prec);
    acb_mul(result, result, rho, prec);
    acb_div(result, result, two_pi, prec);

    arb_clear(omega_val);
    acb_clear(fn_a_val);
    acb_clear(rho);
    acb_clear(two_pi);
    acb_clear(i_omega);
    acb_clear(exp_i_omega);
    acb_clear(i_ph);
    acb_clear(exp_i_ph);
    acb_clear(denominator);
}

/**
 * Integration callback structure for fn_b
 */
typedef struct {
    acb_t a;
    acb_t b;
    acb_t z;
    arb_t phi;
    slong prec;
} fn_b_params_t;

/**
 * Integration callback for fn_b
 */
static int fn_b_integrand(acb_ptr out, const acb_t inp, void* param, slong order, slong prec)
{
    fn_b_params_t* p = (fn_b_params_t*)param;
    arb_t r_val;

    if (order > 1) {
        return 1; /* derivatives not implemented */
    }

    arb_init(r_val);
    acb_get_real(r_val, inp);

    compute_fn_b(out, r_val, p->a, p->b, p->z, p->phi, p->prec);

    arb_clear(r_val);
    return 0;
}

/**
 * Integration callback structure for fn_c
 */
typedef struct {
    acb_t a;
    acb_t b;
    acb_t z;
    arb_t rho;
    slong prec;
} fn_c_params_t;

/**
 * Integration callback for fn_c
 */
static int fn_c_integrand(acb_ptr out, const acb_t inp, void* param, slong order, slong prec)
{
    fn_c_params_t* p = (fn_c_params_t*)param;
    arb_t ph;

    if (order > 1) {
        return 1; /* derivatives not implemented */
    }

    arb_init(ph);
    acb_get_real(ph, inp);

    compute_fn_c(out, ph, p->a, p->b, p->z, p->rho, p->prec);

    arb_clear(ph);
    return 0;
}

/**
 * Quad integration callback for fn_b
 */
static int fn_b_quad_integrand(acb_t out, const acb_t inp, void* param, slong prec)
{
    fn_b_params_t* p = (fn_b_params_t*)param;
    arb_t r_val;

    arb_init(r_val);
    acb_get_real(r_val, inp);

    compute_fn_b(out, r_val, p->a, p->b, p->z, p->phi, p->prec);

    arb_clear(r_val);
    return 0;
}

/**
 * Quad integration callback for fn_c
 */
static int fn_c_quad_integrand(acb_t out, const acb_t inp, void* param, slong prec)
{
    fn_c_params_t* p = (fn_c_params_t*)param;
    arb_t ph;

    arb_init(ph);
    acb_get_real(ph, inp);

    compute_fn_c(out, ph, p->a, p->b, p->z, p->rho, p->prec);

    arb_clear(ph);
    return 0;
}

/**
 * Integral representation
 *
 * Based on algorithm.scm:80-104
 */
void ml_integral_rep(acb_t result, const acb_t a, const acb_t b, const acb_t z,
                     int region, slong prec)
{
    arb_t acc, phi_val, rho_b_val, rib_val, half;
    arb_t log_term, abs_z, abs_b, temp_arb, tol;
    acb_t A_i, B_i, C_i, pi_val;
    slong max_depth;

    arb_init(acc);
    arb_init(phi_val);
    arb_init(rho_b_val);
    arb_init(rib_val);
    arb_init(half);
    arb_init(log_term);
    arb_init(abs_z);
    arb_init(abs_b);
    arb_init(temp_arb);
    acb_init(A_i);
    acb_init(B_i);
    acb_init(C_i);
    acb_init(pi_val);
    arb_init(tol);

    

    /* Set integration parameters */
    max_depth = MITTLEFF_QUAD_INTEGRATION_MAX_DEPTH; /* Match Scheme implementation */
    arb_set_d(tol, MITTLEFF_QUAD_INTEGRATION_TOLERANCE); /* Integration tolerance */

    accuracy_from_prec(acc, prec);
    acb_const_pi(pi_val, prec);
    arb_set_d(half, 0.5);
    acb_abs(abs_z, z, prec);
    acb_abs(abs_b, b, prec);

    /* Φ(i) */
    if (region == 5) {
        /* Φ(5) = π * a */
        acb_mul(pi_val, pi_val, a, prec);
        acb_get_real(phi_val, pi_val);
    } else { /* region == 6 */
        /* Φ(6) = π * a * 2/3 */
        acb_mul(pi_val, pi_val, a, prec);
        acb_mul_2exp_si(pi_val, pi_val, 1); /* multiply by 2 */
        arb_div_ui(acb_realref(pi_val), acb_realref(pi_val), 3, prec);
        acb_get_real(phi_val, pi_val);
    }

    /* Ρ(β) */
    /* TODO: When beta <= 1, the integration lower bound is 0, which causes 0^0 issues
     * in compute_fn_a when beta=1 exactly. This triggers "ulp error not defined for
     * special value" in acb_calc_integrate. Consider using a small epsilon or handling
     * the special case in the integrand functions. */
    arb_t one;
    arb_init(one);
    arb_one(one);

    if (arb_gt(acb_realref(b), one)) {
        arb_set(rho_b_val, half);
    } else {
        arb_zero(rho_b_val);
    }
    // arb_set(rho_b_val, half); // FIX: avoid rho = 0
    arb_clear(one);

    /* Rib(i, β) - Complex computation, simplified version */
    /* This is a simplified calculation - the full Scheme version has many cases */
    arb_t one_val;
    arb_init(one_val);
    arb_one(one_val);
    arb_set_ui(rib_val, 1);
    arb_mul_2exp_si(rib_val, abs_z, 1); /* 2*|z| */
    arb_max(rib_val, rib_val, one_val, prec);
    arb_clear(one_val);

    /* Compute A_i */
    if (region == 5) {
        /* A_i = fn_a(z, a, b, 0) */
        arb_zero(temp_arb);
        compute_fn_a(A_i, z, a, b, temp_arb, prec);
    } else { /* region == 6 */
        acb_zero(A_i);
    }

    /* Compute B_i using integration */
    fn_b_params_t b_params;
    acb_init(b_params.a);
    acb_init(b_params.b);
    acb_init(b_params.z);
    arb_init(b_params.phi);

    acb_set(b_params.a, a);
    acb_set(b_params.b, b);
    acb_set(b_params.z, z);
    arb_set(b_params.phi, phi_val);
    b_params.prec = prec;

    acb_t lower_b, upper_b;
    acb_init(lower_b);
    acb_init(upper_b);
    acb_set_arb(lower_b, rho_b_val);
    acb_set_arb(upper_b, rib_val);

    /* Use quad adaptive Simpson integration */
    quad_complex(B_i, fn_b_quad_integrand, &b_params, lower_b, upper_b,
                 tol, max_depth, prec);

    acb_clear(lower_b);
    acb_clear(upper_b);

    /* Compute C_i using integration if β > 1 */
    arb_init(one);
    arb_one(one);
    if (arb_gt(acb_realref(b), one)) {
        fn_c_params_t c_params;
        acb_init(c_params.a);
        acb_init(c_params.b);
        acb_init(c_params.z);
        arb_init(c_params.rho);

        acb_set(c_params.a, a);
        acb_set(c_params.b, b);
        acb_set(c_params.z, z);
        arb_set(c_params.rho, half);
        c_params.prec = prec;

        acb_t lower_c, upper_c;
        acb_init(lower_c);
        acb_init(upper_c);

        arb_neg(temp_arb, phi_val);
        acb_set_arb(lower_c, temp_arb);
        acb_set_arb(upper_c, phi_val);

        /* Use quad adaptive Simpson integration */
        quad_complex(C_i, fn_c_quad_integrand, &c_params, lower_c, upper_c,
                     tol, max_depth, prec);

        acb_clear(lower_c);
        acb_clear(upper_c);
        acb_clear(c_params.a);
        acb_clear(c_params.b);
        acb_clear(c_params.z);
        arb_clear(c_params.rho);
    } else {
        acb_zero(C_i);
    }

    /* result = A_i + B_i + C_i */
    acb_add(result, A_i, B_i, prec);
    acb_add(result, result, C_i, prec);

    arb_clear(one);
    acb_clear(b_params.a);
    acb_clear(b_params.b);
    acb_clear(b_params.z);
    arb_clear(b_params.phi);
    arb_clear(acc);
    arb_clear(phi_val);
    arb_clear(rho_b_val);
    arb_clear(rib_val);
    arb_clear(half);
    arb_clear(log_term);
    arb_clear(abs_z);
    arb_clear(abs_b);
    arb_clear(temp_arb);
    acb_clear(A_i);
    acb_clear(B_i);
    acb_clear(C_i);
    acb_clear(pi_val);
    arb_clear(tol);
}
