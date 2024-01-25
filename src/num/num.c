#include <num.h>

#include <flint/flint.h>
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_gamma.h>

#include "log.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
// #include <complex.h>

/* Abstract Data Type */
struct num_instance_t
{
    //double* dat;
    acb_t z;
};

//#define NUM_REAL(x)((x)->dat[0])
//#define NUM_IMAG(x)((x)->dat[1])
#define new_max(x,y) (((x) >= (y)) ? (x) : (y))

static int
num_fcmp (const double x, const double y)
{
    return gsl_fcmp (x, y, 1.0e-15);
}

static bool
isclose (const double a, const double b)
{
    /* abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol ) */
    const double rtol = 1.0e-5;
    const double atol = 1.0e-8;

    const double x1 = rtol * new_max(fabs(a), fabs(b));

    //return (abs(a-b) < atol + rtol * abs(b)) ? true : false;
    return (fabs(a-b) <= new_max(x1, fabs(atol))) ? true: false;
    
}

static double
arbtod(const arb_t x)
{
    return arf_get_d(arb_midref(x), ARF_RND_NEAR);
}

static void
acbtonum (num_t res, acb_t x)
{
    /* arb_t re, im; */
    /* arb_init(re), arb_init(im); */
    /* acb_get_real(re, x), acb_get_imag(im, x); */
    /* NUM_REAL(res) = arbtod(re); */
    /* NUM_IMAG(res) = arbtod(im); */
    /* arb_clear(re), arb_clear(im); */
    acb_set(res->z, x);
}


static bool
iszero (const double x)
{
    //return (!isnan(x) && num_fcmp (x, 0.0) == 0) ? true : false;
    return isclose (0.0, x);
}

/* Memory management */
num_t
num_init (void)
{
    num_t self = malloc(sizeof(struct num_instance_t));
    assert (self != NULL);
    /* self -> dat = malloc(2 * sizeof(double)); */
    /* assert (self -> dat); */
    /* self -> dat[0] = 0.0; */
    /* self -> dat[1] = 0.0; */
    return self;
}

void
num_clear (num_t self)
{
    if (self)
    {
        acb_clear(self -> z);
        
        //if (self -> dat) free (self -> dat);
        free (self);
    }
}

/* Input and Output */
void
num_print (const num_t self)
{
    printf("(%+.8e, %+.8e)\n",
           num_real_d (self),
           num_imag_d (self));
}

/* Basic manipulation */

void
num_set_d_d (num_t self, const double x, const double y)
{
    /* self -> dat[0] = x; */
    /* self -> dat[1] = y; */
    acb_set_d_d (self -> z, x, y);
}

void
num_set_d (num_t self, const double x)
{
    num_set_d_d (self, x, 0.0);
}

void
num_set_num (num_t self, const num_t x)
{
    /* self -> dat[0] = num_real_d (x); */
    /* self -> dat[1] = num_imag_d (x); */
    acb_set(self -> z, x -> z);
}

/* void num_zero (num_t self) */
/* { */
    
/* } */

/* /\* Accessors *\/ */
/* void */
/* num_real (num_t res, num_t self) */
/* { */
/*     res -> r = self -> r; */
/*     res -> i = 0.0; */
/* } */

/* void */
/* num_real (num_t res, num_t self) */
/* { */
/*     res -> r = self -> r; */
/*     res -> i = 0.0; */
/* } */

/* double */
/* num_real_d (const num_t self) */
/* { */
/*     return self -> dat[0]; */
/*     acbtonum( */
/* } */

double
num_real_d (const num_t self)
{
    double ret;
    /* const struct num * _self = self; */    
    arb_t x;
    arb_init(x);
    acb_get_real(x, self -> z);
    ret = arbtod(x);
    arb_clear(x);
    return ret;
}


/* double */
/* num_imag_d (const num_t self) */
/* { */
/*     return self -> dat[1]; */
/* } */

double
num_imag_d (const num_t self)
{
    double ret;
    /* const struct num * _self = self; */    
    arb_t x;
    arb_init(x);
    acb_get_imag(x, self -> z);
    ret = arbtod(x);
    arb_clear(x);
    return ret;
}

/* /\* Initialization *\/ */





/* Precision and comparisons */

bool
num_is_zero (const num_t self)
{
    const bool p1 = isnan  (num_real_d (self));
    const bool p2 = isnan  (num_imag_d (self));
    const bool p3 = iszero (num_real_d (self));
    const bool p4 = iszero (num_imag_d (self));
    
    return (!p1 && !p2 && p3 && p4) ? true : false;
}

bool
num_eq (const num_t self, const num_t other)
{
    const double x = num_real_d(self) - num_real_d(other);
    const double y = num_imag_d(self) - num_imag_d(other);
    const double abs2 = x*x + y*y;
    return iszero(abs2);
}

bool
num_eq_d (const num_t self, const double other)
{
    num_t o;
    bool res;
    
    o = num_init();
    num_set_d (o, other);
    res = num_eq (self, o);
    num_clear (o);
    return res;
}

bool
num_ne  (const num_t self, const num_t other)
{
    return !num_eq(self, other);
}

bool num_gt (const num_t self, const num_t other)
{
    assert (num_is_real (self));
    assert (num_is_real (other));
    return (num_fcmp (num_real_d (self), num_real_d (other)) > 0) ? true : false;
}

bool num_gt_d (const num_t self, const double other)
{
    num_t o;
    bool res;
    
    assert (num_is_real (self));
    o = num_init ();
    num_set_d(o, other);
    res = num_gt (self, o);
    num_clear(o);
    return res;
}


bool num_ge (const num_t self, const num_t other)
{
    assert (num_is_real (self));
    assert (num_is_real (other));
    return (num_fcmp (num_real_d (self), num_real_d (other)) >= 0) ? true : false;
}

bool num_ge_d (const num_t self, const double other)
{
    num_t o;
    bool res;
    
    assert (num_is_real (self));
    o = num_init ();
    num_set_d(o, other);
    res = num_ge (self, o);
    num_clear(o);
    return res;
}


bool num_lt (const num_t self, const num_t other)
{
    assert (num_is_real (self));
    assert (num_is_real (other));
    return (num_fcmp (num_real_d (self), num_real_d (other)) < 0) ? true : false;
}

bool num_le (const num_t self, const num_t other)
{
    assert (num_is_real (self));
    assert (num_is_real (other));
    return (num_fcmp (num_real_d (self), num_real_d (other)) <= 0) ? true : false;
}
bool num_le_d (const num_t self, const double other)
{
    num_t o;
    bool res;
    
    assert (num_is_real (self));
    o = num_init ();
    num_set_d(o, other);
    res = num_le (self, o);
    num_clear(o);
    return res;
}

bool
num_is_real (const num_t self)
{
    //log_trace("[%s] %g%+g*I", __func__, NUM_REAL(self), NUM_IMAG(self));
    return (iszero(num_imag_d(self))) ? true : false;
}

/* bool */
/* num_isnan (num_t self) */
/* { */
/*     return (isnan(self -> r) || isnan(self -> i)) ? true : false; */
/* } */

/* Arithmetic */

void
num_neg (num_t res, const num_t self)
{
    /* NUM_REAL(res) = -1.0 * NUM_REAL(x); */
    /* NUM_IMAG(res) = -1.0 * NUM_IMAG(x); */
    acb_neg(res -> z, self -> z);
}

void
num_conj (num_t res, const num_t x)
{
    /* NUM_REAL(res) = +1.0 * NUM_REAL(x); */
    /* NUM_IMAG(res) = -1.0 * NUM_IMAG(x); */
     acb_conj(res -> z, x -> z);
}

void
num_add (num_t res, const num_t x, const num_t y)
{
    
    /* NUM_REAL(res) = NUM_REAL(x) + NUM_REAL(y); */
    /* NUM_IMAG(res) = NUM_IMAG(x) + NUM_IMAG(y); */

    slong prec = 63;
    acb_t _res, _x, _y;
    acb_init(_res), acb_init(_x), acb_init(_y);
    acb_set_d_d(_x, num_real_d(x), num_imag_d(x));
    acb_set_d_d(_y, num_real_d(y), num_imag_d(y));
    acb_add(_res, _x, _y, prec);
    acbtonum(res, _res);
    acb_clear(_res), acb_clear(_x),acb_clear(_y);
}

void
num_sub (num_t res, const num_t x, const num_t y)
{
    /* NUM_REAL(res) = NUM_REAL(x) - NUM_REAL(y); */
    /* NUM_IMAG(res) = NUM_IMAG(x) - NUM_IMAG(y); */

    slong prec = 63;
    acb_t _res, _x, _y;
    acb_init(_res), acb_init(_x), acb_init(_y);
    acb_set_d_d(_x, num_real_d(x), num_imag_d(x));
    acb_set_d_d(_y, num_real_d(y), num_imag_d(y));
    acb_sub(_res, _x, _y, prec);
    acbtonum(res, _res);
    acb_clear(_res), acb_clear(_x),acb_clear(_y);
}

void
num_mul (num_t res, const num_t x, const num_t y)
{
    /* gsl_complex a, b, c; */
    /* GSL_SET_COMPLEX(&a, NUM_REAL(x), NUM_IMAG(x)); */
    /* GSL_SET_COMPLEX(&b, NUM_REAL(y), NUM_IMAG(y)); */
    /* c = gsl_complex_mul (a, b); */
    /* NUM_REAL(res) = GSL_REAL(c); */
    /* NUM_IMAG(res) = GSL_IMAG(c); */

    /* if (num_is_zero (x) || num_is_zero (y)) */
    /*     num_set_d (res, 0.0); */
    /* else { */
    /*     slong prec = 63; */
    /*     acb_t _res, _x, _y; */
    /*     acb_init(_res), acb_init(_x), acb_init(_y); */
    /*     acb_set_d_d(_x, num_real_d(x), num_imag_d(x)); */
    /*     acb_set_d_d(_y, num_real_d(y), num_imag_d(y)); */
    /*     acb_mul(_res, _x, _y, prec); */
    /*     acbtonum(res, _res); */
    /*     acb_clear(_res), acb_clear(_x),acb_clear(_y); */
    /* } */

    slong prec = 53;
        acb_t _res, _x, _y;
        acb_init(_res), acb_init(_x), acb_init(_y);
        acb_set_d_d(_x, num_real_d(x), num_imag_d(x));
        acb_set_d_d(_y, num_real_d(y), num_imag_d(y));
        acb_mul(_res, _x, _y, prec);
        acbtonum(res, _res);
        acb_clear(_res), acb_clear(_x),acb_clear(_y);
}

void
num_mul_d (num_t res, const num_t x, const double y)
{
    /* gsl_complex a, b; */
    /* GSL_SET_COMPLEX(&a, NUM_REAL(x), NUM_IMAG(x)); */
    /* b = gsl_complex_mul_real (a, y); */
    /* NUM_REAL(res) = GSL_REAL(b); */
    /* NUM_IMAG(res) = GSL_IMAG(b); */

    num_t _y;
    _y = num_init();
    num_set_d (_y, y);
    num_mul (res, x, _y);
    num_clear (_y);
}

/* double num_to_d (const num_t self) */
/* { */
/*     //log_trace("%g%+g*I", NUM_REAL(self), NUM_IMAG(self)); */
/*     assert (num_is_real (self)); */
/*     return NUM_REAL(self); */
/* } */

double
num_to_d (const num_t self)
{
    assert(num_is_real(self));

    arb_t x;
    //const struct num * _self = self;
    double res;

    arb_init(x);
    acb_get_real(x, self -> z);

    res = arbtod(x);

    arb_clear (x);
    
    return res;
}


//complex double
/* num_to_complex (const num_t self) */
/* { */
/*     return NUM_REAL(self) + NUM_IMAG(self)*I; */
/* } */

complex double
num_to_complex (const num_t self)
{
    const double re = arbtod(acb_realref(self -> z));
    const double im = arbtod(acb_imagref(self -> z));

    return re + im * I;
}

/* void num_to_d_d (double* res, const num_t x) */
/* { */
/*     res[0] = NUM_REAL(x); */
/*     res[1] = NUM_IMAG(x); */
/* } */

void
num_to_d_d (double* res, const num_t self)
{
    arb_t x, y;
    //const struct num * _self = self;

    arb_init(x), arb_init(y);
    acb_get_real(x, self -> z);
    acb_get_imag(y, self -> z);

    res[0] = arbtod(x), res[1] = arbtod(y);
    arb_clear(x), arb_clear(y);
}

void
num_div (num_t res, const num_t x, const num_t y)
{
    /* gsl_complex a, b, c; */
    /* GSL_SET_COMPLEX(&a, NUM_REAL(x), NUM_IMAG(x)); */
    /* GSL_SET_COMPLEX(&b, NUM_REAL(y), NUM_IMAG(y)); */
    /* c = gsl_complex_div (a, b); */
    /* NUM_REAL(res) = GSL_REAL(c); */
    /* NUM_IMAG(res) = GSL_IMAG(c); */

    slong prec = 63;
    acb_t _res, _x, _y;
    acb_init(_res), acb_init(_x), acb_init(_y);
    acb_set_d_d(_x, num_real_d(x), num_imag_d(x));
    acb_set_d_d(_y, num_real_d(y), num_imag_d(y));
    acb_div(_res, _x, _y, prec);
    acbtonum(res, _res);
    acb_clear(_res), acb_clear(_x),acb_clear(_y);
}

void num_pow (num_t res, const num_t x, const num_t y)
{
    /* gsl_complex a, b, c; */
    /* GSL_SET_COMPLEX(&a, NUM_REAL(x), NUM_IMAG(x)); */
    /* GSL_SET_COMPLEX(&b, NUM_REAL(y), NUM_IMAG(y)); */
    /* c = gsl_complex_pow (a, b); */
    /* NUM_REAL(res) = GSL_REAL(c); */
    /* NUM_IMAG(res) = GSL_IMAG(c); */

    acb_t _res, _x, _y;
    slong prec;

    prec = 63;
    acb_init(_res), acb_init(_x), acb_init(_y);
    acb_set_d_d(_x, num_real_d(x), num_imag_d(x));
    acb_set_d_d(_y, num_real_d(y), num_imag_d(y));
    acb_pow(_res, _x, _y, prec);
    acbtonum(res, _res);
    acb_clear(_res), acb_clear(_x), acb_clear(_y);
}

void num_pow_d (num_t res, const num_t x, const double y)
{
    num_t b;
    b = num_init();
    num_set_d(b, y);
    num_pow(res, x, b);
    num_clear(b);
    /* gsl_complex a, c; */
    /* GSL_SET_COMPLEX(&a, NUM_REAL(x), NUM_IMAG(x)); */
    /* c = gsl_complex_pow_real (a, y); */
    /* NUM_REAL(res) = GSL_REAL(c); */
    /* NUM_IMAG(res) = GSL_IMAG(c); */
}

void
num_fmod (num_t res, const num_t self, const num_t other)
{
    assert(num_is_real (self) && num_is_real (other));

    const double _self = num_to_d(self);
    const double _other = num_to_d(other);
    
    num_set_d(res, fmod(_self, _other));
}

/* void */
/* num_ceil (num_t res, const num_t self) */
/* { */
/*     assert( num_is_real (self) ); */
/*     NUM_REAL(res) = ceil(NUM_REAL(self)); */
/* } */

void
num_ceil (num_t res, const num_t self)
{
    assert(num_is_real(self));
    //struct num * _res = res;
    const double x = num_to_d(self);
    acb_set_d(res -> z, ceil(x));
}


/* void */
/* num_inv (num_t res, const num_t x) */
/* { */
/*     gsl_complex a, c; */
/*     GSL_SET_COMPLEX(&a, NUM_REAL(x), NUM_IMAG(x)); */
/*     c = gsl_complex_inverse (a); */
/*     NUM_REAL(res) = GSL_REAL(c); */
/*     NUM_IMAG(res) = GSL_IMAG(c); */
/* } */

void
num_inv (num_t res, const num_t self)
{
    slong prec = 63;
    /* struct num * _res = res; */
    /* const struct num * _self = self; */
    /* acb_inv(_res -> dat, _self -> dat, PREC); */
    acb_inv(res -> z, self -> z, prec);
}

/* void */
/* num_abs (num_t res, const num_t x) */
/* { */
/*     gsl_complex a; */
/*     GSL_SET_COMPLEX(&a, NUM_REAL(x), NUM_IMAG(x)); */
/*     NUM_REAL(res) = gsl_complex_abs (a); */
/*     NUM_IMAG(res) = 0.0; */
/* } */

void
num_abs (num_t res, const num_t self)
{
    /* struct num * _res = res; */
    /* const struct num * _self = self; */

    slong prec = 63;
    arb_t x;
    arb_init(x);
    acb_abs(x, self -> z, prec);
    acb_set_arb(res -> z, x);
    arb_clear(x);
}


void num_exp (num_t res, const num_t x)
{
    /* gsl_complex a, b; */
    /* GSL_SET_COMPLEX(&a, NUM_REAL(x), NUM_IMAG(x)); */
    /* b = gsl_complex_exp (a); */
    /* NUM_REAL(res) = GSL_REAL(b); */
    /* NUM_IMAG(res) = GSL_IMAG(b); */

    acb_t _res, _x;
    slong prec;

    prec = 63;
    acb_init(_res), acb_init(_x);
    acb_set_d_d(_x, num_real_d(x), num_imag_d(x));
    acb_exp(_res, _x, prec);
    acbtonum(res, _res);
    acb_clear(_res), acb_clear(_x);
}

void num_log (num_t res, const num_t x)
{
    /* gsl_complex a, b; */
    /* GSL_SET_COMPLEX(&a, NUM_REAL(x), NUM_IMAG(x)); */
    /* b = gsl_complex_log (a); */
    /* NUM_REAL(res) = GSL_REAL(b); */
    /* NUM_IMAG(res) = GSL_IMAG(b); */

    acb_t _res, _x;
    slong prec;

    prec = 63;
    acb_init(_res), acb_init(_x);
    acb_set_d_d(_x, num_real_d(x), num_imag_d(x));
    acb_log(_res, _x, prec);
    acbtonum(res, _res);
    acb_clear(_res), acb_clear(_x);
}

void num_sqrt (num_t res, const num_t x)
{
    /* gsl_complex a, b; */
    /* GSL_SET_COMPLEX(&a, NUM_REAL(x), NUM_IMAG(x)); */
    /* b = gsl_complex_sqrt (a); */
    /* NUM_REAL(res) = GSL_REAL(b); */
    /* NUM_IMAG(res) = GSL_IMAG(b); */

    acb_t _res, _x;
    slong prec;

    prec = 63;
    acb_init(_res), acb_init(_x);
    acb_set_d_d(_x, num_real_d(x), num_imag_d(x));
    acb_sqrt(_res, _x, prec);
    acbtonum(res, _res);
    acb_clear(_res), acb_clear(_x);
}

void num_sin (num_t res, const num_t x)
{
    /* gsl_complex a, b; */
    /* GSL_SET_COMPLEX(&a, NUM_REAL(x), NUM_IMAG(x)); */
    /* b = gsl_complex_sin (a); */
    /* NUM_REAL(res) = GSL_REAL(b); */
    /* NUM_IMAG(res) = GSL_IMAG(b); */

    acb_t _res, _x;
    slong prec;

    prec = 63;
    acb_init(_res), acb_init(_x);
    acb_set_d_d(_x, num_real_d(x), num_imag_d(x));
    acb_sin(_res, _x, prec);
    acbtonum(res, _res);
    acb_clear(_res), acb_clear(_x);
}

void num_sinh (num_t res, const num_t x)
{
    /* gsl_complex a, b; */
    /* GSL_SET_COMPLEX(&a, NUM_REAL(x), NUM_IMAG(x)); */
    /* b = gsl_complex_sinh (a); */
    /* NUM_REAL(res) = GSL_REAL(b); */
    /* NUM_IMAG(res) = GSL_IMAG(b); */

    acb_t _res, _x;
    slong prec;

    prec = 63;
    acb_init(_res), acb_init(_x);
    acb_set_d_d(_x, num_real_d(x), num_imag_d(x));
    acb_sinh(_res, _x, prec);
    acbtonum(res, _res);
    acb_clear(_res), acb_clear(_x);
}

void num_cos (num_t res, const num_t x)
{
    /* gsl_complex a, b; */
    /* GSL_SET_COMPLEX(&a, NUM_REAL(x), NUM_IMAG(x)); */
    /* b = gsl_complex_cos (a); */
    /* NUM_REAL(res) = GSL_REAL(b); */
    /* NUM_IMAG(res) = GSL_IMAG(b); */

    acb_t _res, _x;
    slong prec;

    prec = 63;
    acb_init(_res), acb_init(_x);
    acb_set_d_d(_x, num_real_d(x), num_imag_d(x));
    acb_cos(_res, _x, prec);
    acbtonum(res, _res);
    acb_clear(_res), acb_clear(_x);
}

void num_cosh (num_t res, const num_t x)
{
    /* gsl_complex a, b; */
    /* GSL_SET_COMPLEX(&a, NUM_REAL(x), NUM_IMAG(x)); */
    /* b = gsl_complex_cosh (a); */
    /* NUM_REAL(res) = GSL_REAL(b); */
    /* NUM_IMAG(res) = GSL_IMAG(b); */

    acb_t _res, _x;
    slong prec;

    prec = 63;
    acb_init(_res), acb_init(_x);
    acb_set_d_d(_x, num_real_d(x), num_imag_d(x));
    acb_cosh(_res, _x, prec);
    acbtonum(res, _res);
    acb_clear(_res), acb_clear(_x);
}

void num_rgamma (num_t res, const num_t x)
{
    //log_trace ("computing Gamma(%g%+g)", num_real_d(x), num_imag_d(x));
    assert (num_is_real (x));

    /* if (fabs (NUM_REAL(x)) < 20.0) */
    /* { */
    /*     int status; */
    /*     gsl_sf_result result; */
        
    /*     status = gsl_sf_gammainv_e(NUM_REAL(x), &result); */
    /*     log_trace("%d gamma function error: %g\n", status, result.err); */
    /*     assert (status == 0); */
    /*     assert (result.err < 1.0e-12); */
    /*     num_set_d(res, result.val); */
    /* } */
    /* else */
    /* { */
    /*     slong prec = 63; */
    /*     acb_t _res, _x; */
    /*     acb_init(_res), acb_init(_x); */
    /*     acb_set_d_d(_x, NUM_REAL(x), NUM_IMAG(x)); */
    /*     acb_hypgeom_rgamma(_res, _x, prec); */
    /*     acbtonum(res, _res); */
    /*     acb_clear(_res), acb_clear(_x); */
    /* } */

    slong prec = 63;
    acb_t _res, _x;
    acb_init(_res), acb_init(_x);
    //acb_set_d_d(_x, NUM_REAL(x), NUM_IMAG(x));
    acb_hypgeom_rgamma(_res, x->z, prec);
    acbtonum(res, _res);
    acb_clear(_res), acb_clear(_x);
}

void num_erfc (num_t res, const num_t x)
{
    slong prec = 63;
    acb_t _res, _x;
    acb_init(_res), acb_init(_x);
    //acb_set_d_d(_x, NUM_REAL(x), NUM_IMAG(x));
    acb_hypgeom_erfc(_res, x->z, prec);
    acbtonum(res, _res);
    acb_clear(_res), acb_clear(_x);
}

void num_max2 (num_t res, num_t x, num_t y)
{
    if (num_ge (x, y))
        num_set_num (res, x);
    else
        num_set_num (res, y);
}
void num_max3 (num_t res, num_t x, num_t y, num_t z)
{
    num_max2 (res, x, y);
    num_max2 (res, res, z);
}

void num_arg (num_t res, num_t self)
{
    /* gsl_complex a; */
    /* GSL_SET_COMPLEX(&a, NUM_REAL(x), NUM_IMAG(x)); */
    /* NUM_REAL(res) = gsl_complex_arg (a); */
    /* NUM_IMAG(res) = 0.0; */

    slong prec = 63;
    arb_t x;
    arb_init(x);
    acb_arg(x, self -> z, prec);
    acb_set_arb(res -> z, x);
    arb_clear(x);
}

/* /\* Comparison *\/ */



