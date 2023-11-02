#include <assert.h>

#include "unity.h"
#include "mittleff.h"

#include <complex.h>

/* Wrapper around the library function */
double complex
mittleff (const double alpha, const double beta, const double complex z)
{
    const double acc = 1.0e-15;
    
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

void
test_exp_small_z (void)
{
    const double complex z = 0.9;

    const double complex expected = cexp(z);   
    const double complex computed = mittleff(1, 1, z);
    TEST_ASSERT_EQUAL_DOUBLE(creal(expected), creal(computed));
    TEST_ASSERT_EQUAL_DOUBLE(cimag(expected), cimag(computed));
}

int
main (void)
{
    UNITY_BEGIN();
    
    RUN_TEST(test_exp_small_z);

    return UNITY_END();
}

