#include "unity.h"

#include "mittleff.h"
#include "new.h"
#include "num.h"
#include "partition.h"

#include <assert.h>
#include <complex.h>

num_t acc;

typedef double complex complex_t;

#define new_max(x,y) (((x) >= (y)) ? (x) : (y))
int isclose (const complex_t a, const complex_t b)
{
    /* abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol ) */
    const double rtol = 1.0e-5;
    const double atol = 1.0e-8;
    int status = (cabs(a-b) <= new_max(rtol * new_max(cabs(a), cabs(b)), fabs(atol))) ? 1 : 0;
    return status;
}

#define TEST_VALUE(e,c) test_value(e, c)

/* Wrapper around the library function */
complex_t mittleff (const double alpha, const double beta, const complex_t z)
{
    double res[2];
    const int status = mittleff_cmplx(
        res,
        alpha, beta,
        creal(z), cimag(z),
        num_to_d(acc));
    assert(status == 0);
    return res[0] + res[1]*I;
}

void test_value(const complex_t expected, const complex_t computed) {
    bool res;
    
    res = (isclose(expected, computed)) ? true : false;

    TEST_ASSERT_TRUE(res);
}

void setUp (void)
{
    // set stuff up here
    acc = new(num);
    num_set_d(acc, 1.0e-15);
}

void tearDown (void)
{
    // clean stuff up here
    delete(acc);
}

#include "test_partition.h"
#include "test_z_zero.h"
#include "test_exp.h"
#include "test_cos.h"
#include "test_cosh.h"
#include "test_sin.h"
#include "test_sinh.h"
#include "test_siam.h"
#include "test_erfc.h"
#include "test_misc.h"

int
main (void)
{
    UNITY_BEGIN();

    RUN_TEST(test_in_region_G0);
    RUN_TEST(test_in_region_G1);
    RUN_TEST(test_in_region_G2);
    RUN_TEST(test_in_region_G3);
    RUN_TEST(test_in_region_G4);
    RUN_TEST(test_in_region_G5);
    RUN_TEST(test_in_region_G6);

    RUN_TEST(test_z_zero_0);
    RUN_TEST(test_z_zero_1);
    RUN_TEST(test_z_zero_2);

    RUN_TEST(test_exp_0);
    RUN_TEST(test_exp_1);
    RUN_TEST(test_exp_2);

    RUN_TEST(test_cos_0);
    RUN_TEST(test_cos_1);
    RUN_TEST(test_cos_2);
    RUN_TEST(test_cos_3);
    RUN_TEST(test_cos_4);

    RUN_TEST(test_cosh_0);
    RUN_TEST(test_cosh_1);
    RUN_TEST(test_cosh_2);
    RUN_TEST(test_cosh_3);
    RUN_TEST(test_cosh_4);

    RUN_TEST(test_sin_0);
    RUN_TEST(test_sin_1);
    RUN_TEST(test_sin_2);
    RUN_TEST(test_sin_3);
    RUN_TEST(test_sin_4);

    RUN_TEST(test_sinh_0);
    RUN_TEST(test_sinh_1);
    RUN_TEST(test_sinh_2);
    RUN_TEST(test_sinh_3);
    RUN_TEST(test_sinh_4);

    RUN_TEST(test_erfc_0);
    RUN_TEST(test_erfc_1);
    RUN_TEST(test_erfc_2);
    RUN_TEST(test_erfc_3);
    RUN_TEST(test_erfc_4);
    RUN_TEST(test_erfc_5);
    RUN_TEST(test_erfc_6);

    RUN_TEST(test_siam_1);
    RUN_TEST(test_siam_2);
    RUN_TEST(test_siam_3);
    RUN_TEST(test_siam_4);
    RUN_TEST(test_siam_5);
    RUN_TEST(test_siam_6);
    RUN_TEST(test_siam_7);
    RUN_TEST(test_siam_8);
    RUN_TEST(test_siam_9);
    RUN_TEST(test_siam_10);
    RUN_TEST(test_siam_11);
    RUN_TEST(test_siam_12);
    RUN_TEST(test_siam_13);
    RUN_TEST(test_siam_14);
    RUN_TEST(test_siam_15);
    RUN_TEST(test_siam_16);
    RUN_TEST(test_siam_17);
    RUN_TEST(test_siam_18);

    RUN_TEST(test_misc_1);
    RUN_TEST(test_misc_2);
    RUN_TEST(test_misc_3);
    RUN_TEST(test_misc_4);
    RUN_TEST(test_misc_5);
    RUN_TEST(test_misc_6);
    RUN_TEST(test_misc_7);
    RUN_TEST(test_misc_8);
    RUN_TEST(test_misc_9);
    RUN_TEST(test_misc_10);

    return UNITY_END();
}

