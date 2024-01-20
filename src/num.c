/*
 * This file is part of num.c (https://github.com/padawanphysicist/num.c).
 *
 * num.c is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * num.c is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * num.c. If not, see <https://www.gnu.org/licenses/>.
 */

/** 
 * @file num.c
 * @brief Implementation of the Abstract Data Type (ADT).
 */
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <complex.h>
#include <stdarg.h>

#include "abc.h"
#include "new.h"
#include "num.h"

#include <arb.h>
#include <acb.h>
#include <acb_hypgeom.h>

#define PREC 53
#define UNUSED(x) (void)(x)

struct num
{
    const void * class; /* must be first */
    acb_t dat;
};

static void *
num_ctor (void * self, va_list * app)
{
    UNUSED(app);
    struct num * _self = self;
    acb_init(_self -> dat);
    return _self;
}

static void *
num_dtor (void * self)
{
    struct num * _self = self;
    acb_clear(_self -> dat);
    return self;
}

static const struct ABC _num =
{
	sizeof(struct num),
	num_ctor, num_dtor
};

const void * num = & _num;

// Converts an arb_t number to double.
static double
arbtod (const arb_t x)
{
    return arf_get_d(arb_midref(x), ARF_RND_NEAR);
}

/****************************/
/* User interface functions */
/****************************/

/* Basic manipulation */

void
num_zero(num_t self)
{
    struct num * _self = self;
    acb_zero(_self -> dat);
}

void
num_one(num_t self)
{
    struct num * _self = self;
    acb_one(_self -> dat);
}

void
num_onei(num_t self)
{
    struct num * _self = self;
    acb_onei(_self -> dat);
}

void
num_set(num_t self, const num_t other)
{
    struct num * _self = self;
    const struct num * _other = other;
    acb_set(_self -> dat, _other -> dat);
}

static void
num_set_acb(num_t self, const acb_t other)
{
    struct num * _self = self;
    acb_set(_self -> dat, other);
}


void
num_set_d(num_t self, const double x)
{
    struct num * _self = self;
    acb_set_d(_self -> dat, x);
}

void
num_set_d_d(num_t self, const double x, const double y)
{
    struct num * _self = self;
    acb_set_d_d(_self -> dat, x, y);
}

/* Input and Output */
void
num_print (const num_t self, const bool endline)
{
    const struct num * _self = self;
    const slong digits = 8;
    const ulong flags = 0;

    acb_printn(_self -> dat, digits, flags);
    if (endline) puts("\n");
}

/* Accessors: Real and Imaginary parts */
void
num_real (num_t res, const num_t self)
{
    struct num * _res = res;
    const struct num * _self = self;
    
    arb_t x;
    arb_init(x);
    acb_get_real(x, _self -> dat);
    acb_set_arb(_res -> dat, x);
    arb_clear(x);
}
double
num_real_d (const num_t self)
{
    double ret;
    const struct num * _self = self;
    
    arb_t x;
    arb_init(x);
    acb_get_real(x, _self -> dat);
    ret = arbtod(x);
    arb_clear(x);

    return ret;
}

void
num_imag (num_t res, const num_t self)
{
    struct num * _res = res;
    const struct num * _self = self;
    
    arb_t x;
    arb_init(x);
    acb_get_imag(x, _self -> dat);
    acb_set_arb(_res -> dat, x);
    arb_clear(x);
}

double
num_imag_d (const num_t self)
{
    double ret;
    const struct num * _self = self;
    
    arb_t x;
    arb_init(x);
    acb_get_imag(x, _self -> dat);
    ret = arbtod(x);
    arb_clear(x);

    return ret;
}

/* Predicates */
bool
num_is_zero (const num_t self)
{
    const struct num * _self = self;

    const int res = acb_is_zero(_self -> dat);

    return (res != 0) ? true : false;
}

bool
num_is_real (const num_t self)
{
    const struct num * _self = self;

    const int res = acb_is_real(_self -> dat);

    return (res != 0) ? true : false;
}

/* /\* Type casting *\/ */

double
num_to_d (const num_t self)
{
    assert(num_is_real(self));

    arb_t x;
    const struct num * _self = self;
    double res;

    arb_init(x);
    acb_get_real(x, _self -> dat);

    res = arbtod(x);
    
    return res;
}

void
num_to_d_d (double* res, const num_t self)
{
    arb_t x, y;
    const struct num * _self = self;

    arb_init(x), arb_init(y);
    acb_get_real(x, _self -> dat);
    acb_get_imag(y, _self -> dat);

    res[0] = arbtod(x), res[1] = arbtod(y);
    arb_clear(x), arb_clear(y);
}


double complex
num_to_complex (const num_t self)
{
    const struct num * _self = self;
    const double re = arbtod(acb_realref(_self -> dat));
    const double im = arbtod(acb_imagref(_self -> dat));

    return re + im * I;
}

void
num_abs (num_t res, const num_t self)
{
    struct num * _res = res;
    const struct num * _self = self;

    arb_t x;
    arb_init(x);
    acb_abs(x, _self -> dat, PREC);
    acb_set_arb(_res -> dat, x);
    arb_clear(x);
}

void
num_neg (num_t res, const num_t self)
{
    struct num * _res = res;
    const struct num * _self = self;
    acb_neg(_res -> dat, _self -> dat);
}

void
num_inv (num_t res, const num_t self)
{
    struct num * _res = res;
    const struct num * _self = self;
    acb_inv(_res -> dat, _self -> dat, PREC);
}

void
num_conj (num_t res, const num_t self)
{
    struct num * _res = res;
    const struct num * _self = self;
    acb_conj(_res -> dat, _self -> dat);
}

void
num_ceil (num_t res, const num_t self)
{
    assert(num_is_real(self));
    struct num * _res = res;
    const double x = num_to_d(self);
    acb_set_d(_res -> dat, ceil(x));
}

void
num_arg (num_t res, const num_t self)
{
    struct num * _res = res;
    const struct num * _self = self;

    arb_t x;
    arb_init(x);
    acb_arg(x, _self -> dat, PREC);
    acb_set_arb(_res -> dat, x);
    arb_clear(x);
}

void
num_sqrt (num_t res, const num_t self)
{
    struct num * _res = res;
    const struct num * _self = self;
    acb_sqrt(_res -> dat, _self -> dat, PREC);
}

void
num_exp (num_t res, const num_t self)
{
    struct num * _res = res;
    const struct num * _self = self;
    acb_exp(_res -> dat, _self -> dat, PREC);
}

void
num_log (num_t res, const num_t self)
{
    struct num * _res = res;
    const struct num * _self = self;
    acb_log(_res -> dat, _self -> dat, PREC);
}

void
num_sin (num_t res, const num_t self)
{
    struct num * _res = res;
    const struct num * _self = self;
    acb_sin(_res -> dat, _self -> dat, PREC);
}

void
num_sinh (num_t res, const num_t self)
{
    struct num * _res = res;
    const struct num * _self = self;
    acb_sinh(_res -> dat, _self -> dat, PREC);
}

void
num_cos (num_t res, const num_t self)
{
    struct num * _res = res;
    const struct num * _self = self;
    acb_cos(_res -> dat, _self -> dat, PREC);
}

void
num_cosh (num_t res, const num_t self)
{
    struct num * _res = res;
    const struct num * _self = self;
    acb_cosh(_res -> dat, _self -> dat, PREC);
}

/* Binary operations */

/* Arithmetic */

void
num_add (num_t res, const num_t self, const num_t other)
{
    struct num * _res = res;
    const struct num * _self = self;
    const struct num * _other = other;
    acb_add(_res -> dat, _self -> dat, _other -> dat, PREC);
}

void
num_add_d (num_t res, const num_t self, const double other)
{
    num_t o;
    o = new(num);
    num_set_d(o, other);
    num_add(res, self, o);
    delete(o);
}


void
num_sub (num_t res, const num_t self, const num_t other)
{
    struct num * _res = res;
    const struct num * _self = self;
    const struct num * _other = other;
    acb_sub(_res -> dat, _self -> dat, _other -> dat, PREC);
}

void
num_mul (num_t res, const num_t self, const num_t other)
{
    struct num * _res = res;
    const struct num * _self = self;
    const struct num * _other = other;
    acb_mul(_res -> dat, _self -> dat, _other -> dat, PREC);    
}
void
num_mul_d (num_t res, const num_t self, const double other)
{
    num_t o;
    o = new(num);
    num_set_d(o, other);
    num_mul(res, self, o);
    delete(o);
}


void
num_div (num_t res, const num_t self, const num_t other)
{
    struct num * _res = res;
    const struct num * _self = self;
    const struct num * _other = other;
    acb_div(_res -> dat, _self -> dat, _other -> dat, PREC);
}

void
num_fmod (num_t res, const num_t self, const num_t other)
{
    assert(num_is_real(self) && num_is_real(other));

    const double _self = num_to_d(self);
    const double _other = num_to_d(other);
    
    num_set_d(res, fmod(_self, _other));
}

void
num_pow (num_t res, const num_t self, const num_t other)
{
    struct num * _res = res;
    const struct num * _self = self;
    const struct num * _other = other;
    acb_pow(_res -> dat, _self -> dat, _other -> dat, PREC);
}
void
num_pow_d (num_t res, const num_t self, const double other)
{
    num_t o;
    o = new(num);
    num_set_d(o, other);
    num_pow(res, self, o);
    delete(o);
}


/* Logical */

bool
num_eq (const num_t self, const num_t other)
{
    int re_eq, im_eq;
    
    const struct num * _self = self;
    const struct num * _other = other;

    arb_t _self_re, _self_im;
    arb_t _other_re, _other_im;
    arb_init(_self_re), arb_init(_self_im);
    arb_init(_other_re), arb_init(_other_im);
    acb_get_real(_self_re, _self -> dat);
    acb_get_real(_other_re, _other -> dat);
    acb_get_imag(_self_im, _self -> dat);
    acb_get_imag(_other_im, _other -> dat);
    re_eq = arb_eq(_self_re, _other_re);
    im_eq = arb_eq(_self_im, _other_im);
    arb_clear(_self_re), arb_clear(_self_im);
    arb_clear(_other_re), arb_clear(_other_im);

    return (re_eq != 0 && im_eq != 0) ? true : false;
}

bool
num_eq_d (const num_t self, const double other)
{
    num_t _other;
    bool ret;

    _other = new(num);
    num_set_d(_other, other);
    ret = num_eq(self, _other);
    delete(_other);

    return ret;
}

bool
num_lt (const num_t self, const num_t other)
{
    assert(num_is_real(self) && num_is_real(other));

    int is_lt;
    
    const struct num * _self = self;
    const struct num * _other = other;

    arb_t _self_re,_other_re;
    arb_init(_self_re), arb_init(_other_re);
    acb_get_real(_self_re, _self -> dat);
    acb_get_real(_other_re, _other -> dat);
    is_lt = arb_lt(_self_re, _other_re);
    arb_clear(_self_re), arb_clear(_other_re);

    return (is_lt != 0) ? true : false;
}

bool
num_gt (const num_t self, const num_t other)
{
    assert(num_is_real(self) && num_is_real(other));

    int is_gt;
    
    const struct num * _self = self;
    const struct num * _other = other;

    arb_t _self_re,_other_re;
    arb_init(_self_re), arb_init(_other_re);
    acb_get_real(_self_re, _self -> dat);
    acb_get_real(_other_re, _other -> dat);
    is_gt = arb_gt(_self_re, _other_re);
    arb_clear(_self_re), arb_clear(_other_re);

    return (is_gt != 0) ? true : false;
}

bool
num_gt_d (const num_t self, const double other)
{
    num_t _other;
    bool ret;

    _other = new(num);
    num_set_d(_other, other);
    ret = num_gt(self, _other);
    delete(_other);

    return ret;
}

bool
num_le (const num_t self, const num_t other)
{
    assert(num_is_real(self) && num_is_real(other));

    int is_le;
    
    const struct num * _self = self;
    const struct num * _other = other;

    arb_t _self_re,_other_re;
    arb_init(_self_re), arb_init(_other_re);
    acb_get_real(_self_re, _self -> dat);
    acb_get_real(_other_re, _other -> dat);
    is_le = arb_le(_self_re, _other_re);
    arb_clear(_self_re), arb_clear(_other_re);

    return (is_le != 0) ? true : false;
}

bool
num_le_d (const num_t self, const double other)
{
    bool ret;
    num_t _other;

    _other = new(num);
    num_set_d(_other, other);
    ret = num_le(self, _other);
    delete(_other);

    return ret;
}

bool
num_ge (const num_t self, const num_t other)
{
    assert(num_is_real(self) && num_is_real(other));

    int is_ge;
    
    const struct num * _self = self;
    const struct num * _other = other;

    arb_t _self_re,_other_re;
    arb_init(_self_re), arb_init(_other_re);
    acb_get_real(_self_re, _self -> dat);
    acb_get_real(_other_re, _other -> dat);
    is_ge = arb_ge(_self_re, _other_re);
    arb_clear(_self_re), arb_clear(_other_re);

    return (is_ge != 0) ? true : false;
}

bool
num_ge_d (const num_t self, const double other)
{
    bool ret;
    num_t _other;

    _other = new(num);
    num_set_d(_other, other);
    ret = num_ge(self, _other);
    delete(_other);

    return ret;
}

void
num_erf (num_t res, const num_t self)
{
    struct num * _res = res;
    const struct num * _self = self;
    acb_hypgeom_erf(_res -> dat, _self -> dat, PREC);    
}

void
num_erfc (num_t res, const num_t self)
{
    struct num * _res = res;
    const struct num * _self = self;
    acb_hypgeom_erfc(_res -> dat, _self -> dat, PREC);  
}

void
num_rgamma (num_t res, const num_t self)
{
    struct num * _res = res;
    const struct num * _self = self;
    acb_hypgeom_rgamma(_res -> dat, _self -> dat, PREC);    
}

void
num_max (num_t res, const num_t self, const num_t other)
{
    assert(num_is_real(self) && num_is_real(other));
    struct num * _res = res;
    const double _self = num_to_d(self);
    const double _other = num_to_d(other);
    acb_set_d(_res -> dat, (_self > _other) ? _self : _other );
}

void
num_max3 (num_t res, const num_t self, const num_t other, const num_t another)
{
    num_max(res, self, other);
    num_max(res, res, another);
}
