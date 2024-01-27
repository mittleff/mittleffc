#include "unity.h"

#include "mittleff.h"

#include "new.h"
#include "num.h"

#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <complex.h>

#ifdef DEBUG
#include "log.h"
#endif

int ntests = 0;
int ntot = 0;

typedef double complex complex_t;

#define new_max(x,y) (((x) >= (y)) ? (x) : (y))

int isclose (const complex_t a, const complex_t b)
{
    /* abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol ) */
    const double rtol = 1.0e-5;
    const double atol = 1.0e-8;
    int status = (cabs(a-b) <= new_max(rtol * new_max(cabs(a), cabs(b)), fabs(atol))) ? 1 : 0;
#ifdef DEBUG
    log_info("\n[\033[1;33m%s\033[0m] TEST RESULTS:\n\t     \033[1;32mexpected\033[0m = %+.14e%+.14ej,\n\t     \033[1;32mcomputed\033[0m = %+.14e%+.14ej,\n\t     \033[1;32mrel_err\033[0m = %.14e\n",
             __func__,
           creal(a), cimag(a),
           creal(b), cimag(b),
           rel_err);
#endif    
    
}

//#define TEST_VALUE(expected,computed) TEST_ASSERT_TRUE(test_value(expected,computed))
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
TEST_VALUE(const complex_t expected, const complex_t computed)
{
    bool res;
    double tol, rel_err;

    tol = 1e-8;
    rel_err = cabs(expected - computed)/cabs(expected);

#ifdef DEBUG
    log_info("\n[\033[1;33m%s\033[0m] TEST RESULTS:\n\t     \033[1;32mexpected\033[0m = %+.14e%+.14ej,\n\t     \033[1;32mcomputed\033[0m = %+.14e%+.14ej,\n\t     \033[1;32mrel_err\033[0m = %.14e\n",
             __func__,
           creal(expected), cimag(expected),
           creal(computed), cimag(computed),
           rel_err);
#endif
    //res = (rel_err < tol) ? true : false;
    TEST_ASSERT(isclose(expected, computed));

    if (res) ntests += 1;
    ntot += 1;

    return res;
    //TEST_ASSERT_DOUBLE_WITHIN(1e-10, creal(expected), creal(computed));
    //TEST_ASSERT_DOUBLE_WITHIN(1e-10, cimag(expected), cimag(computed));
}

void setUp (void) {}
void tearDown (void) {}

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
    TEST_VALUE( +4.24680224735076e+11, mittleff(0.6, 0.8, 7.0) );
    TEST_VALUE( +4.50513132816147e+64, mittleff(0.6, 0.8, 20.0) );
    //TEST_VALUE( 0.0364029650876388, mittleff(0.6, 0.8, -7.0) );
    TEST_VALUE( 0.0044638678216643, mittleff(0.6, 0.8, -50.0) );
    //TEST_VALUE( 0.00509750816218177+0.0329981074690976*I, mittleff(0.6, 0.8, -2.16311896062463+6.65739561406608*I) );
    TEST_VALUE( 0.00282134530403973+0.0107554765459201*I, mittleff(0.6, 0.8, -6.18033988749895+19.0211303259031*I) );

    TEST_VALUE( 98682128538.563, mittleff(0.6, 1.25, 7.0) );
    TEST_VALUE( 4.76359640442376e+63+9.21339224649432e-19*I, mittleff(0.6, 1.25, 20.0) );
    //TEST_VALUE( 0.101261033685572, mittleff(0.6, 1.25, -7.0) );
    TEST_VALUE( 0.0144197663438114-7.6778268720786e-20*I, mittleff(0.6, 1.25, -50.0) );
    //TEST_VALUE( 0.0333902562082633+0.0980431639835736*I, mittleff(0.6, 1.25, -2.16311896062463+6.65739561406608*I) );
    TEST_VALUE( 0.011289456355613+0.0342852434746551*I, mittleff(0.6, 1.25, -6.18033988749895+19.0211303259031*I) );

    //TEST_VALUE( 76147703794042.9, mittleff(0.6, -0.8, 7.0) );
    TEST_VALUE( 1.32776365747668e+68, mittleff(0.6, -0.8, 20.0) );
    //TEST_VALUE( 0.0501291913319244, mittleff(0.6, -0.8, -7.0) );
    TEST_VALUE( 0.00751163297262774, mittleff(0.6, -0.8, -50.0) );
    //TEST_VALUE( 0.0193182614473201+0.0537209282676945*I, mittleff(0.6, -0.8, -2.16311896062463+6.65739561406608*I) );
    TEST_VALUE( 0.00592228306634142+0.0179734030934208*I, mittleff(0.6, -0.8, -6.18033988749895+19.0211303259031*I) );
}


void
test_Partition (void)
{
    int i;
    num_t z, alpha, acc;
    double complex zvalues_G0[10] =
        {-7.33057219e-02-5.11934762e-01 * I, -2.71670473e-01-3.61159944e-01 * I,
         +2.65104392e-01+4.93166724e-01 * I, +3.83251931e-01-7.20912195e-01 * I,
         +6.15420692e-01+4.94737446e-01 * I, +5.77022431e-01+7.31978260e-01 * I,
         +5.50863360e-01+3.20026771e-01 * I, +3.19070375e-01-6.61764960e-01 * I,
         -4.45577540e-01-6.89742478e-01 * I, +1.59189786e-01+3.86077993e-02 * I};
     double complex zvalues_G1[10] =
        {+1.00809273e+01+2.22251668e+00 * I, +4.22589759e+00+9.54335685e+00 * I,
         +7.07970270e+00-4.66671354e+00 * I, +1.03652185e+01+2.66482346e+00 * I,
         +7.36990152e+00+5.61334844e+00 * I, +1.02336488e+01-1.44223606e+00 * I,
         +8.17467292e+00-2.42910054e+00 * I, +4.34347484e+00+6.98890345e+00 * I,
         +9.93249059e+00+3.48082382e+00 * I, +5.72222055e+00+5.86522403e+00 * I};
     double complex zvalues_G2[10] =
        {-8.81638303e+00+4.53794350e+00 * I, -9.83286816e+00+3.89528981e+00 * I,
         -8.15754405e+00+3.51012135e+00 * I, -9.75448297e+00+4.41974518e+00 * I,
         -1.01655301e+01-1.84058210e+00 * I, -7.74353934e+00-7.18631674e+00 * I,
         -8.35629806e+00+6.52829764e+00 * I, -5.54981503e+00+6.60796993e+00 * I,
         -3.88573272e+00-7.26142380e+00 * I, -7.05943454e+00+7.46849452e+00 * I};
     double complex zvalues_G3[10] =
        {-3.22342758e-01+8.45119872e+00 * I, +7.13347688e-01+1.06864127e+01 * I,
         +4.77462676e-01+1.00665352e+01 * I, +9.88018957e-01+1.05615109e+01 * I,
         -9.62473717e-01+9.61880150e+00 * I, +1.39209663e+00+8.13417968e+00 * I,
         -1.11915976e-01+9.18185740e+00 * I, +1.93827983e+00+1.05413022e+01 * I,
         -9.93493529e-01+1.02465564e+01 * I, +1.13040766e+00+9.12078090e+00 * I};
     double complex zvalues_G4[10] =
        {-3.75588680e-01-9.83203507e+00 * I, -8.88774219e-01-9.55277101e+00 * I,
         -6.41990869e-02-8.27306571e+00 * I, -7.41899625e-02-8.24863509e+00 * I,
         +7.94012305e-01-9.21742744e+00 * I, +2.04734829e+00-1.04412313e+01 * I,
         +1.24541783e+00-9.61437859e+00 * I, +1.57560906e+00-9.41839253e+00 * I,
         +8.01100940e-01-9.12867016e+00 * I, -1.52222011e+00-9.69011955e+00 * I};
      double complex zvalues_G5[10] =
        {+4.08373780e+00+2.53485316e+00 * I, +4.28734692e+00+1.48338558e+00 * I,
         +5.58693005e+00-5.46478678e+00 * I, +4.20513707e+00-1.51605969e-01 * I,
         +1.06269938e+00-3.29532657e+00 * I, +5.31116784e+00-5.90116020e+00 * I,
         +6.98695414e+00+4.17832152e+00 * I, +6.38791794e+00-2.61371742e+00 * I,
         +4.76970575e+00-3.58214698e+00 * I, +2.03837467e+00-3.27614785e+00 * I};
      double complex zvalues_G6[10] =
        {-5.00775165e+00+4.08876443e+00 * I, +1.29206756e+00+6.06330787e+00 * I,
         -3.61542114e+00-6.31484347e+00 * I, -4.41578509e+00+6.29748333e+00 * I,
         -1.98160394e+00+2.44111893e+00 * I, -8.53012232e-02+7.19321871e+00 * I,
         -7.28144099e+00-2.17268099e+00 * I, -6.56039361e+00+4.33423743e+00 * I,
         -3.85833401e+00-4.25315083e+00 * I, -2.97631495e+00+6.48320798e+00 * I};

      z = new(num), alpha = new(num), acc = new(num);
      num_set_d(alpha, 0.5);
      num_set_d(acc, 1.0e-15);

    for(i = 0; i < 10; i++)
    {
        /* region G0 */
        num_set_d_d(z, creal(zvalues_G0[i]), cimag(zvalues_G0[i]));        
        TEST_ASSERT_TRUE(in_region_G0(z));
        /* region G1 */
         num_set_d_d(z, creal(zvalues_G1[i]), cimag(zvalues_G1[i]));        
        TEST_ASSERT_TRUE(in_region_G1(z, alpha, acc));
        /* region G2 */
         num_set_d_d(z, creal(zvalues_G2[i]), cimag(zvalues_G2[i]));        
        TEST_ASSERT_TRUE(in_region_G2(z, alpha, acc));
        /* region G3 */
         num_set_d_d(z, creal(zvalues_G3[i]), cimag(zvalues_G3[i]));        
        TEST_ASSERT_TRUE(in_region_G3(z, alpha, acc));
        /* region G4 */
        num_set_d_d(z, creal(zvalues_G4[i]), cimag(zvalues_G4[i]));        
        TEST_ASSERT_TRUE(in_region_G4(z, alpha, acc));
        /* region G5 */
        num_set_d_d(z, creal(zvalues_G5[i]), cimag(zvalues_G5[i]));        
        TEST_ASSERT_TRUE(in_region_G5(z, alpha, acc));
        /* region G6 */
        num_set_d_d(z, creal(zvalues_G6[i]), cimag(zvalues_G6[i]));        
        TEST_ASSERT_TRUE(in_region_G6(z, alpha, acc));
    }

    delete(z), delete(alpha), delete(acc);
}

void
test_Misc (void)
{
    TEST_VALUE( 1.5403698281390346, mittleff(0.5, 0.5, 0.5) );
    TEST_VALUE( 1.1448466286155243, mittleff(1.5, 0.5, 0.5));
    //TEST_VALUE( 1.201890136368392+0.7895394560075035*I, mittleff(2.3, 1.0, 0.7+2.0*I));
    TEST_VALUE( 1.268233154873853+0.07914994421659409*I, mittleff(2.3, 1.0, 0.7+0.2*I));
    //TEST_VALUE( 8.721285946907744692995882256235296113802695745418015206361825134909144332670706e+2015816, mittleff(0.3, 1.0, 100.0)); //  4641590 terms of the asymptotic series
    TEST_VALUE( -2.7808021618204008e13-2.8561425165239754e13*I, mittleff(0.9, 0.5, 22.0+22.0*I));
    //TEST_VALUE( 0.17617901349590603+2.063981943021305*I, mittleff(0.1, 1.05, 0.9+0.5*I));
    //TEST_VALUE( 1.0358176744122032, mittleff(4.1, 1.0, 1.0));
    TEST_VALUE( 0.046854221014893775, mittleff(0.5, 1.0, -12.0));
    //TEST_VALUE( 0.481952081535048487353320281623, mittleff(0.125, 1.0, -1.0));
}

int main (void)
{
    UNITY_BEGIN();

    RUN_TEST(test_Partition);
    RUN_TEST(test_siam);
    RUN_TEST(test_z_zero);
    RUN_TEST(test_exp);
    RUN_TEST(test_cos_and_cosh);
    RUN_TEST(test_sin_and_sinh);
    RUN_TEST(test_erfc);
    RUN_TEST(test_Misc);

    printf("Number of tests: %d\nNumber of succesful tests: %d", ntot, ntests);
        
    return UNITY_END();
}

