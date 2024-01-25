#include "unity.h"

#include "mittleff.h"

#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <complex.h>

typedef double complex complex_t;

#define new_max(x,y) (((x) >= (y)) ? (x) : (y))

int isclose (const complex_t a, const complex_t b)
{
    /* abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol ) */
    const double rtol = 1.0e-5;
    const double atol = 1.0e-8;
    return (cabs(a-b) <= new_max(rtol * new_max(cabs(a), cabs(b)), fabs(atol))) ? 1 : 0;
    
}

#define TEST_VALUE(expected,computed) test_value(expected,computed);
//#define TEST_VALUE(expected,computed) TEST_ASSERT_DOUBLE_WITHIN(1e-5, creal(expected), creal(computed)); TEST_ASSERT_DOUBLE_WITHIN(1e-5, cimag(expected), cimag(computed));
//#define TEST_VALUE(expected,computed) TEST_ASSERT_EQUAL_DOUBLE(creal(expected), creal(computed)); TEST_ASSERT_EQUAL_DOUBLE(cimag(expected), cimag(computed));
//#define TEST_VALUE(expected,computed) TEST_ASSERT_TRUE(test_value(expected, computed))

/* Wrapper around the library function */
complex_t
mittleff (const double alpha, const double beta, const complex_t z)
{
    double res[2];
    const int status = mittleff_cmplx(
        res,
        alpha, beta,
        creal(z), cimag(z),
        1.0e-15);
    assert(status == 0);
    return res[0] + res[1]*I;
}

bool
test_value(const complex_t expected, const complex_t computed)
{
    bool res;
    double tol, abs_err;

    tol = 1e-13;
    abs_err = cabs(expected - computed)/cabs(expected);
    
    printf("expected = %+.14e%+.14ej, computed = %+.14e%+.14ej, abs_err = %.14e\n",
           creal(expected), cimag(expected),
           creal(computed), cimag(computed),
           abs_err);

    res = (abs_err < tol) ? true : false;

    assert (isclose(expected, computed));
    
    return res;
    //TEST_ASSERT_DOUBLE_WITHIN(1e-10, creal(expected), creal(computed));
    //TEST_ASSERT_DOUBLE_WITHIN(1e-10, cimag(expected), cimag(computed));
}

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
test_z_zero (void)
{
    TEST_VALUE( +1.00000000000000e+00, mittleff(1.0, 1.0, 0.0) );
    TEST_VALUE( +1.08912442105834e+00, mittleff(1.0, 1.2, 0.0) );
    TEST_VALUE( -9.35778720912873e-02, mittleff(1.0, -0.1, 0.0) );
}

void
test_exp (void)
{
    TEST_VALUE( +2.22554092849247e+00+0.00000000000000e+00*I, mittleff(1.0, 1.0, 0.8) );    
    TEST_VALUE( +7.38905609893065e+00+0.00000000000000e+00*I, mittleff(1.0, 1.0, 2.0) );
    TEST_VALUE( -1.31287830814622e+01-1.52007844630680e+01*I, mittleff(1.0, 1.0, 3.0+4.0*I) );
}

void
test_cos_and_cosh (void)
{
    double complex z;

    /* cos */
    z = -1.42469226000000e+00-1.93094380000000e-01*I;
    TEST_VALUE( +1.48307362593645e-01-1.92226474311045e-01*I, mittleff(2.0, 1.0, -z*z) );
    TEST_VALUE( +2.15772016942189e+00+3.75733317193863e-01*I, mittleff(2.0, 1.0,  z*z) );

    z = +6.64679980000000e-01+1.46133978000000e+00*I;
    TEST_VALUE( +1.78818881219233e+00-1.25822734251227e+00*I, mittleff(2.0, 1.0, -z*z) );
    TEST_VALUE( +1.34270384074874e-01+7.10437973043830e-01*I, mittleff(2.0, 1.0,  z*z) );

    z = +8.94857660000000e-01+1.89827224000000e+00*I;
    TEST_VALUE( +2.13470542364750e+00-2.54495655003334e+00*I, mittleff(2.0, 1.0, -z*z) );
    TEST_VALUE( -4.59266292056142e-01+9.64999829685609e-01*I, mittleff(2.0, 1.0,  z*z) );

    z = -1.32190510000000e-01+2.22985315000000e+00*I;
    TEST_VALUE( +4.66199107535284e+00+6.05710930000310e-01*I, mittleff(2.0, 1.0, -z*z) );
    TEST_VALUE( -6.17729663704252e-01-1.04810500062338e-01*I, mittleff(2.0, 1.0,  z*z) );

    z = -1.60041154000000e+00+3.95133590000000e-01*I;
    TEST_VALUE( -3.19526988999422e-02+4.05318417916536e-01*I, mittleff(2.0, 1.0, -z*z) );
    TEST_VALUE( +2.37976084266865e+00-9.14839293944298e-01*I, mittleff(2.0, 1.0,  z*z) );

}

void
test_sin_and_sinh (void)
{
    double complex z;

    z = -1.42469226000000e+00-1.93094380000000e-01*I;
    TEST_VALUE( -1.00784724885197e+00-2.82866292071741e-02*I, z*mittleff(2.0, 2.0, -z*z) );
    TEST_VALUE( -1.92160887420070e+00-4.21900298087567e-01*I, z*mittleff(2.0, 2.0,  z*z) );

    z = +6.64679980000000e-01+1.46133978000000e+00*I;
    TEST_VALUE( +1.40128057728241e+00+1.60563708193143e+00*I, z*mittleff(2.0, 2.0, -z*z) );
    TEST_VALUE( +7.80741317390111e-02+1.22179750677017e+00*I, z*mittleff(2.0, 2.0,  z*z) );

    z = +8.94857660000000e-01+1.89827224000000e+00*I;
    TEST_VALUE( +2.66183979862789e+00+2.04096901440268e+00*I, z*mittleff(2.0, 2.0, -z*z) );
    TEST_VALUE( -3.27817271209846e-01+1.35194796777751e+00*I, z*mittleff(2.0, 2.0,  z*z) );

    z = -1.32190510000000e-01+2.22985315000000e+00*I;
    TEST_VALUE( -6.19885888230589e-01+4.55538511767942e+00*I, z*mittleff(2.0, 2.0, -z*z) );
    TEST_VALUE( +8.11856608972997e-02+7.97487564190490e-01*I, z*mittleff(2.0, 2.0,  z*z) );

    z = -1.60041154000000e+00+3.95133590000000e-01*I;
    TEST_VALUE( -1.07861309811772e+00-1.20071018874968e-02*I, z*mittleff(2.0, 2.0, -z*z) );
    TEST_VALUE( -2.19349810100817e+00+9.92523644338986e-01*I, z*mittleff(2.0, 2.0,  z*z) );

}

void
test_erfc (void)
{
    int i;
    
    double complex z[7] = {
        -7.33057219000000e-02-5.11934762000000e-01*I,
        +1.00809273000000e+01+2.22251668000000e+00*I,
        -8.81638303000000e+00+4.53794350000000e+00*I,
        -3.22342758000000e-01+8.45119872000000e+00*I,
        -3.75588680000000e-01-9.83203507000000e+00*I,
        +4.08373780000000e+00+2.53485316000000e+00*I,
        -5.00775165000000e+00+4.08876443000000e+00*I
    };

    double complex expected[7] = {
        +7.25165453880409e-01-4.32914368208589e-01*I,
        +1.32220943230009e+42+1.43926327412783e+42*I,
        +5.05454404812233e-02+2.57564201381802e-02*I,
        +2.59774904698260e-03+6.71347813921331e-02*I,
        +2.22360956770567e-03-5.75980077079615e-02*I,
        -1.58178109396067e+04+5.43930514682910e+04*I,
        +6.80477811127566e-02+5.42607316062483e-02*I
    };

    
    for (i = 0; i < 7; i++) {
        TEST_VALUE( expected[i], mittleff(0.5, 1.0, z[i]) );
    }
}

void
test_siam (void)
{
    int i;
    double alpha = 0.6;
    double beta[18] = {
        0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
        1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 
        -0.8, -0.8, -0.8, -0.8, -0.8, -0.8
    };
    double z[18] = {
        7.0, 20.0, -7.0, -50.0,
        -2.16311896062463+6.65739561406608*I,
        -6.18033988749895+19.0211303259031*I,
        7.0, 20.0, -7.0, -50.0,
        -2.16311896062463+6.65739561406608*I,
        -6.18033988749895+19.0211303259031*I,
        7.0, 20.0, -7.0, -50.0,
        -2.16311896062463+6.65739561406608*I,
        -6.18033988749895+19.0211303259031*I
    };
    double expected[18] = {
        +4.24680224735076e+11,
        +4.50513132816147e+64,        
        0.0364029650876388,
        0.0044638678216643,
        0.00509750816218177+0.0329981074690976*I,    
        0.00282134530403973+0.0107554765459201*I,
        98682128538.563,
        4.76359640442376e+63+9.21339224649432e-19*I,
        0.101261033685572,
        0.0144197663438114-7.6778268720786e-20*I,
        0.0333902562082633+0.0980431639835736*I,
        0.011289456355613+0.0342852434746551*I,
        76147703794042.9,
        1.32776365747668e+68,
        0.0501291913319244,
        0.00751163297262774,
        0.0193182614473201+0.0537209282676945*I,
        0.00592228306634142+0.0179734030934208*I
    };

    for (i = 0; i < 18; i++) {
        if (i != 2) TEST_VALUE( expected[i], mittleff(alpha, beta[i], z[i]) );
    }

}




/* #include "test_partition.h" */
/* #include "test_siam.h" */
/* #include "test_misc.h" */

int
main (void)
{
    UNITY_BEGIN();

    //RUN_TEST(test_siam_3);
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

    RUN_TEST(test_z_zero);
    RUN_TEST(test_exp);
    RUN_TEST(test_cos_and_cosh);
    RUN_TEST(test_sin_and_sinh);
    RUN_TEST(test_erfc);
    RUN_TEST(test_siam);
    
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

