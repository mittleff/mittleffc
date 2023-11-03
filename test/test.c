#include <assert.h>

#include "unity.h"
#include "log.h"
#include "mittleff.h"

#include <complex.h>

const double acc = 1.0e-10;

#define TEST_VALUE(expected,computed) TEST_ASSERT_DOUBLE_WITHIN(acc, creal(expected), creal(computed)); TEST_ASSERT_DOUBLE_WITHIN(acc, cimag(expected), cimag(computed));

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

/********************************/
/* Tests for the bignum library */
/********************************/

/*******************************/
/* Tests for the main function */
/*******************************/

void
test_exp_small_z (void)
{
    const double complex z = 0.9;

    log_info("[%s] Testing whether E(1, 1, z) == exp(z)", __func__);

    const double complex expected = cexp(z);   
    const double complex computed = mittleff(1, 1, z);
    TEST_VALUE(expected, computed);
}

void
test_sin_small_z (void)
{
    const double complex z = 0.9;

    log_info("[%s] Testing whether z*E(2, 2, -z**2) == sin(z)", __func__);

    const double complex expected = csin(z);   
    const double complex computed = z*mittleff(2, 2, -z*z);
    TEST_VALUE(expected, computed);
}

void
test_cos_small_z (void)
{
    const double complex z = 0.9;

    log_info("[%s] Testing whether E(2, 1, -z**2) == cos(z)", __func__);

    const double complex expected = ccos(z);   
    const double complex computed = mittleff(2, 1, -z*z);
    TEST_VALUE(expected, computed);
}

void
test_sinh_small_z (void)
{
    const double complex z = 0.9;

    log_info("[%s] Testing whether z*E(2, 2, z**2) == sinh(z)", __func__);

    const double complex expected = csinh(z);   
    const double complex computed = z*mittleff(2, 2, z*z);
    TEST_VALUE(expected, computed);
}

void
test_cosh_small_z (void)
{
    const double complex z = 0.9;

    log_info("[%s] Testing whether E(2, 1, z**2) == cosh(z)", __func__);

    const double complex expected = ccosh(z);   
    const double complex computed = mittleff(2, 1, z*z);
    TEST_VALUE(expected, computed);
}

void
test_sin_recurrence (void) /* Uses G5 */
{
    TEST_IGNORE();
    const double complex z = 1.01;

    log_info("[%s] Testing whether z*E(2, 2, -z**2) == sin(z)", __func__);

    const double complex expected = csin(z);   
    const double complex computed = z*mittleff(2, 2, -z*z);
    TEST_VALUE(expected, computed);
}

void
test_cos_recurrence (void) /* Uses G5 and G6 */
{
    const double complex z = 1.001;

    log_info("[%s] Testing whether E(2, 1, -z**2) == cos(z)", __func__);

    const double complex expected = ccos(z);   
    const double complex computed = mittleff(2, 1, -z*z);
    TEST_VALUE(expected, computed);
}


int
main (void)
{
    UNITY_BEGIN();
    
    RUN_TEST(test_exp_small_z);
    RUN_TEST(test_sin_small_z);
    RUN_TEST(test_cos_small_z);
    RUN_TEST(test_sinh_small_z);
    RUN_TEST(test_cosh_small_z);
    RUN_TEST(test_sin_recurrence);
    RUN_TEST(test_cos_recurrence);

    return UNITY_END();
}

