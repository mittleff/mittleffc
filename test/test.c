#include <assert.h>

#include "mittleff.h"
#include "new.h"
#include "num.h"

#include "unity.h"
#include "log.h"

#include <complex.h>

const double acc = 1.0e-10;

#define TEST_VALUE(expected,computed) TEST_ASSERT_DOUBLE_WITHIN(acc, creal(expected), creal(computed)); TEST_ASSERT_DOUBLE_WITHIN(1e-12, cimag(expected), cimag(computed));

/* Wrapper around the library function */
double complex
mittleff (const double alpha, const double beta, const double complex z)
{
    double res[2];
    const int status = mittleff_cmplx(
        res,
        alpha, beta,
        creal(z), cimag(z),
        acc);
    assert(status == 0);
    return res[0] + res[1]*I;
}

double complex
derf (const double complex x)
{
    /* compute erf(x) */
    double erf_x[2];
    num_t res = new(num);
    num_set_d_d(res, creal(x), cimag(x));
    num_erf(res, res);
    num_to_d_d(erf_x, res);
    delete(res);
    return erf_x[0] + erf_x[1] * I;
}

void
setUp (void)
{
    // set stuff up here
}

void
tearDown (void)
{
    // clean stuff up here
}


/*******************************/
/* Tests for the main function */
/*******************************/

/* sage: z = 0.9; exp(z**2) * (1 + erf(z)) */
/* 4.03928432202983 */
void
test_G0_real_parameter (void) /* Uses G0 */
{
    const double complex x = 0.9;
    const double complex expected = cexp(x * x) * (1.0 + derf(x));
    const double complex computed = mittleff(0.5, 1.0, x);
    TEST_VALUE(expected, computed);
}

/* sage: z = 0.46 + 0.76*I; exp(z**2) * (1 + erf(z)) */
/* 0.608499571536833 + 1.21446566539644*I */
void
test_G0_complex_parameter (void) 
{
    const double complex x = 0.46 + 0.76 * I;
    const double complex expected = cexp(x * x) * (1.0 + derf(x));
    const double complex computed = mittleff(0.5, 1.0, x);
    TEST_VALUE(expected, computed);
}

/* sage: z = 5.0 - 6.0*I; exp(z**2) * (1 + erf(z)) */
/* -0.0467872933620735 - 0.0551794173038091*I */
void
test_G1_complex_parameter (void)
{
    const double complex x = 5.0 - 6.0 * I;
    const double complex expected = cexp(x * x) * (1.0 + derf(x));
    printf("%g %g\n", expected);
    const double complex computed = mittleff(0.5, 1.0, x);
    TEST_VALUE(expected, computed);
}

/* void */
/* test_region_0 (void) */
/* { */
/*     const double complex z = +5.31167097e-01 + 1.11251357e-01 * I; */
/*     const double complex expected = cexp(z); */
/*     const double complex computed = mittleff(1, 1, z); */
/*     TEST_VALUE(expected, computed); */
/* } */

/* void */
/* test_region_1 (void) */
/* { */
/*     const double complex z = +9.48638023e+01 + 3.65839418e+00 * I; */
/*     const double complex expected = cexp(z); */
/*     const double complex computed = mittleff(1, 1, z); */
/*     TEST_VALUE(expected, computed); */
/* } */

/* void */
/* test_exp_small_z (void) */
/* { */
/*     const double complex z = 0.9; */
/*     const double complex expected = cexp(z); */
/*     const double complex computed = mittleff(1, 1, z); */
/*     TEST_VALUE(expected, computed); */
/* } */

/* void */
/* test_sin_small_z (void) */
/* { */
/*     const double complex z = 0.9; */
/*     const double complex expected = csin(z); */
/*     const double complex computed = z*mittleff(2, 2, -z*z); */
/*     TEST_VALUE(expected, computed); */
/* } */

/* void */
/* test_cos_small_z (void) */
/* { */
/*     const double complex z = 0.9; */
/*     const double complex expected = ccos(z); */
/*     const double complex computed = mittleff(2, 1, -z*z); */
/*     TEST_VALUE(expected, computed); */
/* } */

/* void */
/* test_sinh_small_z (void) */
/* { */
/*     const double complex z = 0.9; */
/*     const double complex expected = csinh(z); */
/*     const double complex computed = z*mittleff(2, 2, z*z); */
/*     TEST_VALUE(expected, computed); */
/* } */

/* void */
/* test_cosh_small_z (void) */
/* { */
/*     const double complex z = 0.9; */
/*     const double complex expected = ccosh(z); */
/*     const double complex computed = mittleff(2, 1, z*z); */
/*     TEST_VALUE(expected, computed); */
/* } */

/* void */
/* test_sin_recurrence (void) /\* Uses G5 and G6 *\/ */
/* { */
/*     TEST_IGNORE(); */
/*     const double complex z = 1.001; */
/*     const double complex expected = csin(z); */
/*     const double complex computed = z*mittleff(2, 2, -z*z); */
/*     TEST_VALUE(expected, computed); */
/* } */

/* void */
/* test_cos_recurrence (void) /\* Uses G5 and G6 *\/ */
/* { */
/*     TEST_IGNORE(); */
/*     const double complex z = 1.001; */
/*     const double complex expected = ccos(z); */
/*     const double complex computed = mittleff(2, 1, -z*z); */
/*     TEST_VALUE(expected, computed); */
/* } */


/* void */
/* test_erf_small_z (void) /\* Uses G5 and G6 *\/ */
/* { */
/*     const double complex x = 0.9; */

/*     /\* compute erf(x) *\/ */
/*     double erf_x; */
/*     num_t res = new(num); */
/*     num_set_d(res, x); */
/*     num_erf(res, res); */
/*     erf_x = num_to_d(res); */
/*     delete(res); */
    
/*     const double complex expected = exp(x * x) * (1.0 + erf_x) + 0.0 * I; */
/*     const double complex computed = mittleff(0.5, 1.0, x + 0.0 * I); */
/*     TEST_VALUE(expected, computed); */
/* } */

/* void */
/* test_erf_not_small_z (void) /\* Uses G5 *\/ */
/* { */
/*     const double complex x = 1.01; */

/*     /\* compute erf(x) *\/ */
/*     double erf_x; */
/*     num_t res = new(num); */
/*     num_set_d(res, x); */
/*     num_erf(res, res); */
/*     erf_x = num_to_d(res); */
/*     delete(res); */
    
/*     const double complex expected = exp(x * x) * (1.0 + erf_x) + 0.0 * I; */
/*     const double complex computed = mittleff(0.5, 1.0, x + 0.0 * I); */
/*     TEST_VALUE(expected, computed); */
/* } */


int
main (void)
{
    UNITY_BEGIN();

    RUN_TEST(test_G0_real_parameter);
    RUN_TEST(test_G0_complex_parameter);

    RUN_TEST(test_G1_complex_parameter);
    //RUN_TEST(test_region_1);

    /* RUN_TEST(test_exp_small_z); */
    /* RUN_TEST(test_sin_small_z); */
    /* RUN_TEST(test_cos_small_z); */
    /* RUN_TEST(test_sinh_small_z); */
    /* RUN_TEST(test_cosh_small_z); */
    /* RUN_TEST(test_erf_small_z); */
    
    /* RUN_TEST(test_sin_recurrence); */
    /* RUN_TEST(test_cos_recurrence); */
    
    /* RUN_TEST(test_erf_not_small_z); */

    return UNITY_END();
}

