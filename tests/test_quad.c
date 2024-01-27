#include "unity.h"

#include "new.h"
#include "num.h"
#include "quad.h"

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

void
fn_sin (num_t res,
        const num_t x,
        void * ctx)
{
    num_sin(res, x);
}

void
test_quad_Function_Call (void)
{
    num_t x, y;

    x = new(num), y = new(num);

    num_set_d(x, 1.57079632679490);
    fn_sin(y, x, NULL);

    TEST_ASSERT_TRUE( num_is_real(y) );
    TEST_ASSERT_EQUAL_DOUBLE( 1.0, num_to_d(y) );

    delete(x),delete(y);
}

/* void */
/* test_quad_Function_Call_f_integrand (void) */
/* { */
/*     num_function_t F; */
    
/*     F.function = &fn_sin; */
/*     F.params = NULL; */

/*     acb_t x, res; */
/*     acb_init(x), acb_init(res); */
/*     acb_set_d(x, 1.57079632679490); */
/*     acb_printn(x, 10, 10); puts("\n"); */
/*     int s = f_integrand(res, x, &F, 0, 100); */
/*     acb_printn(res, 10, 10); puts("\n"); */
/*     acb_clear(res), acb_clear(x);     */
/* } */

void
test_quad_Integrate_Sin_0_100 (void)
{
    num_t res, a, b;
    num_function_t F;

    res = new(num);
    a = new(num), b = new(num);
    num_set_d(a, 0.0);
    num_set_d(b, 100.0);
    
    F.function = &fn_sin;
    F.params = NULL;

    quad(res, &F, a, b);

    TEST_ASSERT_DOUBLE_WITHIN (1e-12, 0.137681127712316, num_real_d(res));

    delete(a), delete(b);
    delete(res);
}

int
main (void)
{
    UNITY_BEGIN();

    RUN_TEST(test_quad_Function_Call);
    /* RUN_TEST(test_quad_Function_Call_f_integrand); */
    RUN_TEST(test_quad_Integrate_Sin_0_100);

    return UNITY_END();
}

