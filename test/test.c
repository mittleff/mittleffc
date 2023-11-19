#include "unity.h"

#include "mittleff.h"
#include "new.h"
#include "num.h"
#include "partition.h"

#include <assert.h>
#include <complex.h>

//const double acc = 1.0e-15;

// #define TEST_VALUE(expected,computed) TEST_ASSERT_DOUBLE_WITHIN(1e-8, creal(expected), creal(computed)); TEST_ASSERT_DOUBLE_WITHIN(1e-8, cimag(expected), cimag(computed));
#define TEST_VALUE(expected,computed) TEST_ASSERT_TRUE(test_value(expected, computed))

/* Wrapper around the library function */
/* double complex */
/* mittleff (const double alpha, const double beta, const double complex z) */
/* { */
/*     double res[2]; */
/*     const int status = mittleff_cmplx( */
/*         res, */
/*         alpha, beta, */
/*         creal(z), cimag(z), */
/*         acc); */
/*     assert(status == 0); */
/*     return res[0] + res[1]*I; */
/* } */

/* bool */
/* test_value(const double complex expected, const double complex computed) */
/* { */
/*     bool res; */
/*     double tol, abs_err; */

/*     tol = 1e-14; */
/*     abs_err = fabs(expected - computed)/fabs(expected); */
    
/*     printf("expected = %+.8e%+.8ej, computed = %+.8e%+.8ej, abs_err = %.8e\n", */
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
}

void
tearDown (void)
{
    // clean stuff up here
}

void test_partition (void)
{
    int i, region;
    num_t z, alpha, acc;
        
    double complex zvalues[6][10] =
        {
            /* region G0 */
            {-7.33057219e-02-5.11934762e-01 * I, -2.71670473e-01-3.61159944e-01 * I,
             +2.65104392e-01+4.93166724e-01 * I, +3.83251931e-01-7.20912195e-01 * I,
             +6.15420692e-01+4.94737446e-01 * I, +5.77022431e-01+7.31978260e-01 * I,
             +5.50863360e-01+3.20026771e-01 * I, +3.19070375e-01-6.61764960e-01 * I,
             -4.45577540e-01-6.89742478e-01 * I, +1.59189786e-01+3.86077993e-02 * I},
            /* region G1 */
            {+1.00809273e+01+2.22251668e+00 * I, +4.22589759e+00+9.54335685e+00 * I,
             +7.07970270e+00-4.66671354e+00 * I, +1.03652185e+01+2.66482346e+00 * I,
             +7.36990152e+00+5.61334844e+00 * I, +1.02336488e+01-1.44223606e+00 * I,
             +8.17467292e+00-2.42910054e+00 * I, +4.34347484e+00+6.98890345e+00 * I,
             +9.93249059e+00+3.48082382e+00 * I, +5.72222055e+00+5.86522403e+00 * I}
        };

    z = new(num), alpha = new(num), acc = new(num);
    num_set_d(alpha, 0.5), num_set_d(acc, 1.0e-15);
    for (region = 0; region < 2; region++)
    {
        printf("region %d\n", region);
        for(i = 0; i < 10; i++)
        {
            num_set_d_d(z, creal(zvalues[region][i]), cimag(zvalues[region][i]));
            printf("z = %g %g\n", creal(zvalues[region][i]), cimag(zvalues[region][i]));
            switch(region)
            {
            case 0:
                TEST_ASSERT_TRUE(in_region_G0(z));
                break;
            case 1:
                TEST_ASSERT_TRUE(in_region_G1(z, alpha, acc));
                break;
            case 2:
                TEST_ASSERT_TRUE(in_region_G2(z, alpha, acc));
                break;
            case 3:
                TEST_ASSERT_TRUE(in_region_G3(z, alpha, acc));
                break;
            case 4:
                TEST_ASSERT_TRUE(in_region_G4(z, alpha, acc));
                break;
            case 5:
                TEST_ASSERT_TRUE(in_region_G5(z, alpha, acc));
                break;
            default:
                TEST_ASSERT_TRUE(in_region_G6(z, alpha, acc));
            }
        }
    }
    delete(z), delete(alpha), delete(acc);

    
    
    
    
    

    /* /\* in_region_G2 *\/ */
    /* TEST_ASSERT_TRUE(in_region_G2(-8.81638303e+00, +4.53794350e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G2(-9.83286816e+00, +3.89528981e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G2(-8.15754405e+00, +3.51012135e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G2(-9.75448297e+00, +4.41974518e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G2(-1.01655301e+01, -1.84058210e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G2(-7.74353934e+00, -7.18631674e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G2(-8.35629806e+00, +6.52829764e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G2(-5.54981503e+00, +6.60796993e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G2(-3.88573272e+00, -7.26142380e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G2(-7.05943454e+00, +7.46849452e+00, alpha, acc)); */

    /* /\* in_region_G3 *\/ */
    /* TEST_ASSERT_TRUE(in_region_G3(-3.22342758e-01, +8.45119872e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G3(+7.13347688e-01, +1.06864127e+01, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G3(+4.77462676e-01, +1.00665352e+01, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G3(+9.88018957e-01, +1.05615109e+01, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G3(-9.62473717e-01, +9.61880150e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G3(+1.39209663e+00, +8.13417968e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G3(-1.11915976e-01, +9.18185740e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G3(+1.93827983e+00, +1.05413022e+01, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G3(-9.93493529e-01, +1.02465564e+01, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G3(+1.13040766e+00, +9.12078090e+00, alpha, acc)); */

    /* /\* in_region_G4 *\/ */
    /* TEST_ASSERT_TRUE(in_region_G4(-3.75588680e-01, -9.83203507e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G4(-8.88774219e-01, -9.55277101e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G4(-6.41990869e-02, -8.27306571e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G4(-7.41899625e-02, -8.24863509e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G4(+7.94012305e-01, -9.21742744e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G4(+2.04734829e+00, -1.04412313e+01, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G4(+1.24541783e+00, -9.61437859e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G4(+1.57560906e+00, -9.41839253e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G4(+8.01100940e-01, -9.12867016e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G4(-1.52222011e+00, -9.69011955e+00, alpha, acc)); */

    /* /\* in_region_G5 *\/ */
    /* TEST_ASSERT_TRUE(in_region_G5(+4.08373780e+00, +2.53485316e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G5(+4.28734692e+00, +1.48338558e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G5(+5.58693005e+00, -5.46478678e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G5(+4.20513707e+00, -1.51605969e-01, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G5(+1.06269938e+00, -3.29532657e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G5(+5.31116784e+00, -5.90116020e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G5(+6.98695414e+00, +4.17832152e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G5(+6.38791794e+00, -2.61371742e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G5(+4.76970575e+00, -3.58214698e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G5(+2.03837467e+00, -3.27614785e+00, alpha, acc)); */

    /* /\* in_region_G6 *\/ */
    /* TEST_ASSERT_TRUE(in_region_G6(-5.00775165e+00, +4.08876443e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G6(+1.29206756e+00, +6.06330787e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G6(-3.61542114e+00, -6.31484347e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G6(-4.41578509e+00, +6.29748333e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G6(-1.98160394e+00, +2.44111893e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G6(-8.53012232e-02, +7.19321871e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G6(-7.28144099e+00, -2.17268099e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G6(-6.56039361e+00, +4.33423743e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G6(-3.85833401e+00, -4.25315083e+00, alpha, acc)); */
    /* TEST_ASSERT_TRUE(in_region_G6(-2.97631495e+00, +6.48320798e+00, alpha, acc)); */
}


/*******************************/
/* Tests for the main function */
/*******************************/

/************************************************************************/
/* sage: z = -7.33057219e-02-5.11934762e-01*I; exp(z**2) * (1 + erf(z)) */
/* # => 0.725165453880409 - 0.432914368208589*I                         */
/************************************************************************/
/* void */
/* test_mittleff0 (void)  */
/* { */
/*     const double complex x = -7.33057219e-02-5.11934762e-01*I; */
/*     const double complex expected = ref_value(x); */
/*     const double complex computed = mittleff(0.5, 1.0, x); */
/*     TEST_VALUE(expected, computed); */
/* } */

/* /\************************************************************************\/ */
/* /\* sage: z = +1.00809273e+01+2.22251668e+00*I; exp(z**2) * (1 + erf(z)) *\/ */
/* /\* # => 1.32220943230009e42 + 1.43926327412783e42*I                     *\/ */
/* /\************************************************************************\/ */
/* void */
/* test_mittleff1 (void) */
/* { */
/*     const double complex x = +1.00809273e+01+2.22251668e+00*I; */
/*     const double complex expected = ref_value(x); */
/*     const double complex computed = mittleff(0.5, 1.0, x); */
/*     TEST_VALUE(expected, computed); */
/* } */

/* /\************************************************************************\/ */
/* /\* sage: z = -8.81638303e+00+4.53794350e+00*I; exp(z**2) * (1 + erf(z)) *\/ */
/* /\* -113372.198954288 - 10701.8161144899*I                               *\/ */
/* /\************************************************************************\/ */
/* void */
/* test_mittleff2 (void) */
/* { */
/*     //TEST_IGNORE_MESSAGE("not working yet"); */
/*     const double complex x = -8.81638303e+00+4.53794350e+00*I; */
/*     const double complex expected = ref_value(x); */
/*     const double complex computed = mittleff(0.5, 1.0, x); */
/*     TEST_VALUE(expected, computed); */
/* } */

/* /\************************************************************************\/ */
/* /\* sage: z = -3.22342758e-01+8.45119872e+00*I; exp(z**2) * (1 + erf(z)) *\/ */
/* /\* #=> 0.00259774904698259 + 0.0671347813921331*I                       *\/ */
/* /\************************************************************************\/ */
/* void */
/* test_mittleff3 (void) */
/* { */
/*     //TEST_IGNORE_MESSAGE("not working yet"); */
/*     const double complex x = -3.22342758e-01+8.45119872e+00*I; */
/*     const double complex expected = ref_value(x); */
/*     const double complex computed = mittleff(0.5, 1.0, x); */
/*     TEST_VALUE(expected, computed); */
/* } */

/* /\************************************************************************\/ */
/* /\* sage: z = -3.75588680e-01-9.83203507e+00*I; exp(z**2) * (1 + erf(z)) *\/ */
/* /\* # => 0.00222360956770567 - 0.0575980077079615*I                      *\/ */
/* /\************************************************************************\/ */
/* void */
/* test_mittleff4 (void) */
/* { */
/*     TEST_IGNORE(); */
/*     const double complex x = -3.75588680e-01-9.83203507e+00*I; */
/*     const double complex expected = ref_value(x); */
/*     const double complex computed = mittleff(0.5, 1.0, x); */
/*     TEST_VALUE(expected, computed); */
/* } */

/* /\************************************************************************\/ */
/* /\* sage: z = +4.08373780e+00+2.53485316e+00*I; exp(z**2) * (1 + erf(z)) *\/ */
/* /\* # => -15817.8109396066 + 54393.0514682909*I                          *\/ */
/* /\************************************************************************\/ */
/* void */
/* test_mittleff5 (void) */
/* { */
/*     TEST_IGNORE(); */
/*     const double complex x = +4.08373780e+00+2.53485316e+00*I; */
/*     const double complex expected = ref_value(x); */
/*     const double complex computed = mittleff(0.5, 1.0, x); */
/*     TEST_VALUE(expected, computed); */
/* } */

/* /\************************************************************************\/ */
/* /\* sage: z = -5.00775165e+00+4.08876443e+00*I; exp(z**2) * (1 + erf(z)) *\/ */
/* /\* # => 0.0680477811130904 + 0.0542607316062114*I                       *\/ */
/* /\************************************************************************\/ */
/* void */
/* test_mittleff6 (void) */
/* { */
/*     TEST_IGNORE(); */
/*     const double complex x = -5.00775165e+00+4.08876443e+00*I; */
/*     const double complex expected = ref_value(x); */
/*     const double complex computed = mittleff(0.5, 1.0, x); */
/*     TEST_VALUE(expected, computed); */
/* } */

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

    RUN_TEST(test_partition);

    /* RUN_TEST(test_mittleff0); */
    /* RUN_TEST(test_mittleff1); */
    /* RUN_TEST(test_mittleff2); */
    /* RUN_TEST(test_mittleff3); */
    /* RUN_TEST(test_mittleff4); */
    /* RUN_TEST(test_mittleff5); */
    /* RUN_TEST(test_mittleff6); */
    /* RUN_TEST(test_G0_complex_parameter); */

    /* RUN_TEST(test_G1_complex_parameter); */

    /* RUN_TEST(test_G2_complex_parameter); */

    
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

