#include "unity.h"

#include "num.h"

#include <assert.h>
#include <complex.h>
#include <stdbool.h>

num_t acc;

//#define TEST_VALUE(expected,computed) TEST_ASSERT_DOUBLE_WITHIN(1e-8, creal(expected), creal(computed)); TEST_ASSERT_DOUBLE_WITHIN(1e-8, cimag(expected), cimag(computed));
// #define TEST_VALUE(expected,computed) TEST_ASSERT_EQUAL_DOUBLE(creal(expected), creal(computed)); TEST_ASSERT_EQUAL_DOUBLE(cimag(expected), cimag(computed));
//#define TEST_VALUE(expected,computed) TEST_ASSERT_TRUE(test_value(expected, computed))

/* Wrapper around the library function */
/* double complex */
/* mittleff (const double alpha, const double beta, const double complex z) */
/* { */
/*     double res[2]; */
/*     const int status = mittleff_cmplx( */
/*         res, */
/*         alpha, beta, */
/*         creal(z), cimag(z), */
/*         num_to_d(acc)); */
/*     assert(status == 0); */
/*     return res[0] + res[1]*I; */
/* } */

/* bool */
/* test_value(const double complex expected, const double complex computed) */
/* { */
/*     bool res; */
/*     double tol, abs_err; */

/*     tol = 1e-13; */
/*     abs_err = fabs(expected - computed)/fabs(expected); */
    
/*     printf("expected = %+.14e%+.14ej, computed = %+.14e%+.14ej, abs_err = %.14e\n", */
/*            creal(expected), cimag(expected), */
/*            creal(computed), cimag(computed), */
/*            abs_err); */

/*     res = (abs_err < tol) ? true : false; */
    
/*     return res; */
/*     //TEST_ASSERT_DOUBLE_WITHIN(1e-10, creal(expected), creal(computed)); */
/*     //TEST_ASSERT_DOUBLE_WITHIN(1e-10, cimag(expected), cimag(computed)); */
/* } */

void
setUp (void)
{
    // set stuff up here
    /* acc = new(num); */
    /* num_set_d(acc, 1.0e-15); */
}

void
tearDown (void)
{
    // clean stuff up here
    /* delete(acc); */
}

void
test_Memory_Management (void)
{
    num_t x;

    x = num_init ();
    TEST_ASSERT_NOT_NULL(x);

    num_clear (x);
}

void
test_Initialization_and_Accessors (void)
{
    num_t x, y;

    x = num_init ();
    TEST_ASSERT_NOT_NULL(x);

    num_set_d (x, +2.0);
    TEST_ASSERT_EQUAL_DOUBLE( +2.0, num_real_d (x) );
    TEST_ASSERT_EQUAL_DOUBLE(  0.0, num_imag_d (x) );

    num_set_d_d (x, +3.0, -4.0);
    TEST_ASSERT_EQUAL_DOUBLE( +3.0, num_real_d (x) );
    TEST_ASSERT_EQUAL_DOUBLE( -4.0, num_imag_d (x) );

    y = num_init ();
    num_set_num (y, x);
    TEST_ASSERT_EQUAL_DOUBLE( +3.0, num_real_d (y) );
    TEST_ASSERT_EQUAL_DOUBLE( -4.0, num_imag_d (y) );

    num_clear (x), num_clear (y);
}

void
test_Precision_and_Comparisons (void)
{
    num_t x, y;

    x = num_init (), y = num_init ();
    TEST_ASSERT_NOT_NULL(x);
    TEST_ASSERT_NOT_NULL(y);
    TEST_ASSERT(num_is_zero(x));
    TEST_ASSERT(num_is_zero(y));

    /* x = 3.0, y = 3.0, z = x */
    num_set_d (x, 3.0), num_set_d (y, 3.0);
    TEST_ASSERT(num_eq (x, y));

    num_set_d_d (x, 2.0, -10.0);
    TEST_ASSERT_FALSE(num_eq (x, y));
    TEST_ASSERT(num_ne (x, y));

    num_set_d (x, 0.0), num_set_d (y, 0.0);
    TEST_ASSERT(num_eq (x, y));

    num_set_d (x, 2.0);
    TEST_ASSERT( num_is_real (x) );

    num_set_d_d (x, 2.0, -10.0);
    TEST_ASSERT_FALSE( num_is_real (x) );

    num_set_d (x, 2.0), num_set_d (y, -3.0);
    TEST_ASSERT_TRUE( num_gt (x, y) );
    TEST_ASSERT_TRUE( num_lt (y, x) );
    TEST_ASSERT_TRUE( num_ge (x, y) );
    TEST_ASSERT_TRUE( num_le (y, x) );

    num_set_d (x, 2.0), num_set_d (y, 2.0);
    TEST_ASSERT_FALSE( num_gt (x, y) );
    TEST_ASSERT_FALSE( num_lt (y, x) );
    TEST_ASSERT_TRUE( num_ge (x, y) );
    TEST_ASSERT_TRUE( num_le (y, x) );

    /* Issue */
    num_set_d_d (x, 4.94623e-05, +4.62743e-21);
    TEST_ASSERT_TRUE( num_is_real (x) );

    num_clear (x), num_clear(y);
}

void
test_Arithmetic (void)
{
    num_t x, y, z;

    x = num_init (), y = num_init ();
    TEST_ASSERT_NOT_NULL(x);
    TEST_ASSERT_NOT_NULL(y);

    num_set_d_d (x, 3.0, -4.0);
    TEST_ASSERT_EQUAL_DOUBLE( +3.0, num_real_d (x) );
    TEST_ASSERT_EQUAL_DOUBLE( -4.0, num_imag_d (x) );
    num_neg (x, x);
    TEST_ASSERT_EQUAL_DOUBLE( -3.0, num_real_d (x) );
    TEST_ASSERT_EQUAL_DOUBLE( +4.0, num_imag_d (x) );
    num_conj (x, x);
    TEST_ASSERT_EQUAL_DOUBLE( -3.0, num_real_d (x) );
    TEST_ASSERT_EQUAL_DOUBLE( -4.0, num_imag_d (x) );

    
    z = num_init ();
    num_set_d_d (x, +3.0, +4.0), num_set_d_d (y, +3.0, -4.0);

    /* add */
    num_add (z, x, y);
    TEST_ASSERT_EQUAL_DOUBLE( +6.0, num_real_d (z) );
    TEST_ASSERT_EQUAL_DOUBLE(  0.0, num_imag_d (z) );

    /* sub */
    num_sub (z, x, y);
    TEST_ASSERT_EQUAL_DOUBLE(  0.0, num_real_d (z) );
    TEST_ASSERT_EQUAL_DOUBLE( +8.0, num_imag_d (z) );

    /* mul */
    num_mul (z, x, y);
    TEST_ASSERT_EQUAL_DOUBLE( +25.0, num_real_d (z) );
    TEST_ASSERT_EQUAL_DOUBLE( 0.0, num_imag_d (z) );

    // Issue when multiplying large and small values

    /* div */
    num_div (z, x, y);
    TEST_ASSERT_EQUAL_DOUBLE( -0.28, num_real_d (z) );
    TEST_ASSERT_EQUAL_DOUBLE( +0.96, num_imag_d (z) );

    num_clear (x), num_clear (y), num_clear (z);
}

void
test_Special_Functions (void)
{
    num_t x, y;

    x = num_init (), y = num_init ();
    TEST_ASSERT_NOT_NULL(x);
    TEST_ASSERT_NOT_NULL(y);

    num_set_d (x, +1.0), num_rgamma (y, x);
    TEST_ASSERT_EQUAL_DOUBLE( +1.0, num_real_d (y) );
    TEST_ASSERT_EQUAL_DOUBLE(  0.0, num_imag_d (y) );

    num_set_d (x, +4.0), num_rgamma (y, x);
    TEST_ASSERT_EQUAL_DOUBLE( 1.0/6.0, num_real_d (y) );
    TEST_ASSERT_EQUAL_DOUBLE( 0.0, num_imag_d (y) );

    num_clear (x), num_clear (y);
}

/* void */
/* test_num_comparison (void) */
/* { */
/*     num_t x, y; */

/*     x = num_init (), y = num_init (); */
/*     TEST_ASSERT_NOT_NULL(x); */
/*     TEST_ASSERT_NOT_NULL(y); */

/*     /\* x = 3.0, y = 3.0, z = x *\/ */
/*     num_set_d (x, 3.0), num_set_d (y, 3.0); */

/*     num_clear (&x), num_clear (&y), num_clear (&z); */
/*     /\* num_print(x); *\/ */
/*     /\* num_print(y); *\/ */
/*     /\* num_print(z); *\/ */
/*     TEST_ASSERT_NULL(x); */
/*     TEST_ASSERT_NULL(y); */
/*     TEST_ASSERT_NULL(z); */
/* } */

/* #include "test_partition.h" */
/* #include "test_z_zero.h" */
/* #include "test_exp.h" */
/* #include "test_cos.h" */
/* #include "test_cosh.h" */
/* #include "test_sin.h" */
/* #include "test_sinh.h" */
/* #include "test_siam.h" */
/* #include "test_erfc.h" */
/* #include "test_misc.h" */

int
main (void)
{
    UNITY_BEGIN();

    RUN_TEST(test_Memory_Management);
    RUN_TEST(test_Initialization_and_Accessors);
    RUN_TEST(test_Precision_and_Comparisons);
    RUN_TEST(test_Arithmetic);
    RUN_TEST(test_Special_Functions);
    
    /* RUN_TEST(test_siam_4); */
    /* RUN_TEST(test_siam_5); */
    /* RUN_TEST(test_siam_6); */
    /* RUN_TEST(test_siam_8); */
    /* RUN_TEST(test_siam_9); */
    /* RUN_TEST(test_siam_10); */
    /* RUN_TEST(test_siam_11); */
    /* RUN_TEST(test_siam_12); */
    /* RUN_TEST(test_siam_15); */
    /* RUN_TEST(test_siam_16); */
    /* RUN_TEST(test_siam_17); */
    /* RUN_TEST(test_siam_18); */
    /* RUN_TEST(test_misc_3); */
    /* RUN_TEST(test_misc_7); */
    /* RUN_TEST(test_misc_8); */
    /* RUN_TEST(test_misc_10); */

    /* RUN_TEST(test_in_region_G0); */
    /* RUN_TEST(test_in_region_G1); */
    /* RUN_TEST(test_in_region_G2); */
    /* RUN_TEST(test_in_region_G3); */
    /* RUN_TEST(test_in_region_G4); */
    /* RUN_TEST(test_in_region_G5); */
    /* RUN_TEST(test_in_region_G6); */

    /* RUN_TEST(test_z_zero_0); */
    /* RUN_TEST(test_z_zero_1); */
    /* RUN_TEST(test_z_zero_2); */

    /* RUN_TEST(test_exp_0); */
    /* RUN_TEST(test_exp_1); */
    /* RUN_TEST(test_exp_2); */

    /* RUN_TEST(test_cos_0); */
    /* RUN_TEST(test_cos_1); */
    /* RUN_TEST(test_cos_2); */
    /* RUN_TEST(test_cos_3); */
    /* RUN_TEST(test_cos_4); */

    /* RUN_TEST(test_cosh_0); */
    /* RUN_TEST(test_cosh_1); */
    /* RUN_TEST(test_cosh_2); */
    /* RUN_TEST(test_cosh_3); */
    /* RUN_TEST(test_cosh_4); */

    /* RUN_TEST(test_sin_0); */
    /* RUN_TEST(test_sin_1); */
    /* RUN_TEST(test_sin_2); */
    /* RUN_TEST(test_sin_3); */
    /* RUN_TEST(test_sin_4); */

    /* RUN_TEST(test_sinh_0); */
    /* RUN_TEST(test_sinh_1); */
    /* RUN_TEST(test_sinh_2); */
    /* RUN_TEST(test_sinh_3); */
    /* RUN_TEST(test_sinh_4); */

    /* RUN_TEST(test_siam_1); */
    /* RUN_TEST(test_siam_2); */
    /* RUN_TEST(test_siam_3); */
    /* RUN_TEST(test_siam_4); */
    /* RUN_TEST(test_siam_5); */
    /* RUN_TEST(test_siam_6); */
    /* RUN_TEST(test_siam_7); */
    /* RUN_TEST(test_siam_8); */
    /* RUN_TEST(test_siam_9); */
    /* RUN_TEST(test_siam_10); */
    /* RUN_TEST(test_siam_11); */
    /* RUN_TEST(test_siam_12); */
    /* RUN_TEST(test_siam_13); */
    /* RUN_TEST(test_siam_14); */
    /* RUN_TEST(test_siam_15); */
    /* RUN_TEST(test_siam_16); */
    /* RUN_TEST(test_siam_17); */
    /* RUN_TEST(test_siam_18); */
    
    /* RUN_TEST(test_erfc_0); */
    /* RUN_TEST(test_erfc_1); */
    /* RUN_TEST(test_erfc_2); */
    /* RUN_TEST(test_erfc_3); */
    /* RUN_TEST(test_erfc_4); */
    /* RUN_TEST(test_erfc_5); */
    /* RUN_TEST(test_erfc_6); */

    /* RUN_TEST(test_misc_1); */
    /* RUN_TEST(test_misc_2); */
    /* RUN_TEST(test_misc_3); */
    /* RUN_TEST(test_misc_4); */
    /* RUN_TEST(test_misc_5); */
    /* RUN_TEST(test_misc_6); */
    /* RUN_TEST(test_misc_7); */
    /* RUN_TEST(test_misc_8); */
    /* RUN_TEST(test_misc_9); */
    /* RUN_TEST(test_misc_10); */


    return UNITY_END();
}

