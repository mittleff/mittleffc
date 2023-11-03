/*
 * This file is part of mittleffc (https://github.com/mittleff/mittleffc).
 *
 * mittleffc is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * mittleffc is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * mittleffc. If not, see <https://www.gnu.org/licenses/>.
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

#include "log.h"

#include <arb.h>
#include <arb_hypgeom.h>

#include "arb2num.h"

#include "abc.h"
#include "new.h"
// #include "num.h"
//include "log.h"

/* #ifndef _TOLERANCE */
/* #define _TOLERANCE 1.0e-8 */
/* #endif */

const int prec = 53;

#define NCOMP 2

struct num {
	const void* class;	/* must be first */
	double* dat;        /* data */
};

static void*
num_ctor (void* _self, va_list* app)
{
    struct num* self = _self;
	
	self->dat = malloc(NCOMP*sizeof(double));
	assert(self->dat);

    for (int i = 0; i < NCOMP; i++)
    {
        const double res = va_arg(* app, const double);
        self->dat[i] = res;
    }

	return self;
}

static void*
num_dtor (void* _self)
{
    struct num* self = _self;

	free(self->dat), self->dat = NULL;
	return self;
}

static void*
num_clone (const void* _self)
{
    const struct num* self = _self;
	return new(num, self->dat);
}

static const struct ABC _num = {
	sizeof(struct num),
	num_ctor, num_dtor,	num_clone
};

const void * num = & _num;

#define scalar_zero 0.0

/****************************/
/* User interface functions */
/****************************/

/* Accessors: Real and Imaginary parts */
num_t
num_real (const num_t _self)
{
    const struct num* self = _self;
	return new(num, self->dat[0], scalar_zero);
}

num_t
num_imag (const num_t _self)
{
    const struct num* self = _self;
	return new(num, self->dat[1], scalar_zero);
}

double
num_real_d (const num_t _self)
{
    const struct num* self = _self;
	return self->dat[0];
}

double
num_imag_d (const num_t _self)
{
    const struct num* self = _self;
	return self->dat[1];
}

/* Predicates */
bool
num_is_zero (const num_t _self)
{
    const struct num* self = _self;

    return ((fabs(self->dat[0]) < _TOLERANCE) && (fabs(self->dat[1]) < _TOLERANCE));
}

bool
num_is_real (const num_t _self)
{
    return num_is_zero(num_imag(_self));
}

/* Type casting */

/* static double */
/* arbtod (const arb_t x); */

/* static num_t */
/* num_from_acb (const acb_t x); */

/* static num_t */
/* num_from_arb (const arb_t x); */

double
num_to_double (const num_t _self)
{
    assert(num_is_real(_self));
    const struct num* self = _self;
    return self->dat[0];
}

double complex
num_to_complex (const num_t _self)
{
    const struct num* self = _self;
    return self->dat[0] + self->dat[1] * I;
}


/* void */
/* num_to_pair (double* res, const num_t _self) */
/* { */
/*     const struct num* self = _self; */

/*     for (int i = 0; i < NCOMP; i++) */
/*         res[i] = self->dat[i]; */
/* } */

/* /\* Unary operations *\/ */
/* num_t */
/* num_abs2 (const num_t _self) */
/* { */
/*     const struct num* self = _self; */
/*     const double x = self->dat[0]; */
/*     const double y = self->dat[1]; */

/*     return new(num, x*x + y*y, 0.0); */
/* } */

num_t
num_abs (const num_t _self)
{
    num_t res;
    acb_t self;
    arb_t _res;
    
    acb_init(self); arb_init(_res);
    acb_set_d_d(self, num_real_d(_self), num_imag_d(_self));
    arb_set_d(_res, 0.0);
    acb_abs(_res, self, prec);
    res = num_from_arb(_res);
    acb_clear(self); arb_clear(_res);

    return res;
}

num_t
num_negative (const num_t _self)
{
    const struct num* self = _self;
    
	return new(num, -1 * self->dat[0], -1 * self->dat[1]);
}

num_t
num_inverse (const num_t _self)
{
    assert(!num_is_zero(_self));
    num_t res;
    acb_t _res, self;

    acb_init(self);
    acb_set_d_d(self, num_real_d(_self), num_imag_d(_self));
    
    acb_inv(_res, self, prec);
    
    res = num_from_acb(_res);
    
    acb_clear(self);

    return res;
}

/* num_t */
/* num_conjugate (const num_t _self) */
/* { */
/*     const struct num * self = _self; */
    
/* 	return new(num, self->dat[0], -1 * self->dat[1]); */
/* } */

num_t
num_arg (const num_t _self)
{
    num_t res;
    acb_t self;
    arb_t _res;

    acb_init(self); arb_init(_res);
    acb_set_d_d(self, num_real_d(_self), num_imag_d(_self));
    arb_set_d(_res, 0.0);       
    acb_arg(_res, self, prec);
    res = num_from_arb(_res);
    acb_clear(self); arb_clear(_res);

    return res;
}

num_t
num_sqrt (const num_t _self)
{
    num_t res;
    acb_t _res, self;

    acb_init(self);
    acb_set_d_d(self, num_real_d(_self), num_imag_d(_self));
    
    acb_sqrt(_res, self, prec);
    
    res = num_from_acb(_res);
    
    acb_clear(self);

    return res;
}

num_t
num_exp (const num_t _self)
{
    num_t res;
    acb_t _res, self;

    acb_init(self);
    acb_set_d_d(self, num_real_d(_self), num_imag_d(_self));
    
    acb_exp(_res, self, prec);
    
    res = num_from_acb(_res);
    
    acb_clear(self);

    return res;
}

num_t
num_log (const num_t _self)
{    
    num_t res;
    acb_t _res, self;

    acb_init(self); acb_init(_res);
    acb_set_d_d(self, num_real_d(_self), num_imag_d(_self));    
    acb_log_analytic(_res, self, 1, prec);
    res = num_from_acb(_res);
    acb_clear(self); acb_clear(_res);

    return res;
}

/* sin(x+iy) = sin(x) cosh(y) + i cos(x) sinh(y)) */
num_t
num_sin (const num_t _self)
{
      num_t res;
    acb_t _res, self;

    acb_init(self);
    acb_set_d_d(self, num_real_d(_self), num_imag_d(_self));
    
    acb_sin(_res, self, prec);

    
    res = num_from_acb(_res);
    
    acb_clear(self); //acb_clear(_res);
    //log_trace("%g %g", num_real_d(res), num_imag_d(res));

    return res;
}

/* cos(x+iy)=cos(x) cosh(y) − i sin(x) sinh(y) */
num_t
num_cos (const num_t _self)
{
    num_t res;
    acb_t _res, self;

    acb_init(self);
    acb_set_d_d(self, num_real_d(_self), num_imag_d(_self));
    
    acb_cos(_res, self, prec);

    
    res = num_from_acb(_res);
    
    acb_clear(self); //acb_clear(_res);
    //log_trace("%g %g", num_real_d(res), num_imag_d(res));

    return res;
}

/* Binary operations */

/* Arithmetic */

num_t
num_add (const num_t _self, const num_t _other)
{
    num_t res;
    acb_t self, other;
    acb_t _res;
    
    acb_init(self); acb_init(other); acb_init(_res);
    acb_set_d_d(self, num_real_d(_self), num_imag_d(_self));
    acb_set_d_d(other, num_real_d(_other), num_imag_d(_other));
    acb_add(_res, self, other, prec);
    res = num_from_acb(_res);
    acb_clear(self); acb_clear(other); acb_clear(_res);

    return res;
}

num_t
num_sub (const num_t _self, const num_t _other)
{
    num_t res;
    acb_t self, other;
    acb_t _res;
    
    acb_init(self); acb_init(other); acb_init(_res);
    acb_set_d_d(self, num_real_d(_self), num_imag_d(_self));
    acb_set_d_d(other, num_real_d(_other), num_imag_d(_other));
    acb_sub(_res, self, other, prec);
    res = num_from_acb(_res);
    acb_clear(self); acb_clear(other); acb_clear(_res);

    return res;
}

num_t
num_mul (const num_t _self, const num_t _other)
{
    num_t res;
    acb_t self, other, _res;
    
    acb_init(self); acb_init(other); acb_init(_res);
    acb_set_d_d(self, num_real_d(_self), num_imag_d(_self));
    acb_set_d_d(other, num_real_d(_other), num_imag_d(_other));
    acb_mul(_res, self, other, prec);
    res = num_from_acb(_res);
    acb_clear(self); acb_clear(other); acb_clear(_res);

    return res;
}

num_t
num_div (const num_t _self, const num_t _other)
{
    num_t res;
    acb_t self, other, _res;
    
    acb_init(self); acb_init(other); acb_init(_res);
    acb_set_d_d(self, num_real_d(_self), num_imag_d(_self));
    acb_set_d_d(other, num_real_d(_other), num_imag_d(_other));
    acb_div(_res, self, other, prec);
    res = num_from_acb(_res);
    acb_clear(self); acb_clear(other); acb_clear(_res);

    return res;
}

num_t
num_fmod (const num_t _self, const num_t _other)
{
    assert(num_is_real(_self) && num_is_real(_other));
    return new(num, fmod(num_to_double(_self), num_to_double(_other)), 0.0);
}

num_t
num_pow (const num_t _self, const num_t _other)
{
    num_t res;
    acb_t self, other, _res;
    
    acb_init(self); acb_init(other); acb_init(_res);
    acb_set_d_d(self, num_real_d(_self), num_imag_d(_self));
    acb_set_d_d(other, num_real_d(_other), num_imag_d(_other));
    acb_pow(_res, self, other, prec);
    res = num_from_acb(_res);
    acb_clear(self); acb_clear(other); acb_clear(_res);

    return res;
}

/* Logical */

bool
num_eq (const num_t _self, const num_t _other)
{
    bool res;
    acb_t self, other;

    acb_init(self); acb_init(other);
    acb_set_d_d(self, num_real_d(_self), num_imag_d(_self));
    acb_set_d_d(other, num_real_d(_other), num_imag_d(_other));
    res = acb_eq(self, other);
    acb_clear(self); acb_clear(other);

    return res;
}

bool
num_lt (const num_t _self, const num_t _other)
{
    assert(num_is_real(_self) && num_is_real(_other));
    
    bool res;
    arb_t self, other;

    arb_init(self); arb_init(other);
    arb_set_d(self, num_real_d(_self));
    arb_set_d(other, num_real_d(_other));
    res = arb_lt(self, other);
    arb_clear(self); arb_clear(other);

    return res;
}

bool
num_gt (const num_t _self, const num_t _other)
{
    assert(num_is_real(_self) && num_is_real(_other));
    
    bool res;
    arb_t self, other;

    arb_init(self); arb_init(other);
    arb_set_d(self, num_real_d(_self));
    arb_set_d(other, num_real_d(_other));
    res = arb_gt(self, other);
    arb_clear(self); arb_clear(other);

    return res;
}

bool
num_le (const num_t _self, const num_t _other)
{
    assert(num_is_real(_self) && num_is_real(_other));
   
    bool res;
    arb_t self, other;

    arb_init(self); arb_init(other);
    arb_set_d(self, num_real_d(_self));
    arb_set_d(other, num_real_d(_other));
    res = arb_le(self, other);
    arb_clear(self); arb_clear(other);

    return res;
}

bool
num_ge (const num_t _self, const num_t _other)
{
    assert(num_is_real(_self) && num_is_real(_other));

    bool res;
    arb_t self, other;

    arb_init(self); arb_init(other);
    arb_set_d(self, num_real_d(_self));
    arb_set_d(other, num_real_d(_other));
    res = arb_ge(self, other);
    arb_clear(self); arb_clear(other);

    return res;
}

num_t
num_max (const num_t _self, const num_t _other)
{
    assert(num_is_real(_self) && num_is_real(_other));
    
    num_t res;
    arb_t _res, self, other;

    arb_init(self); arb_init(other);
    arb_set_d(self, num_real_d(_self));
    arb_set_d(other, num_real_d(_other));
    arb_max(_res, self, other, prec);
    res = num_from_arb(_res);
    arb_clear(self); arb_clear(other); arb_clear(_res);

    return res;
}

num_t
num_ceil (const num_t _self)
{
    assert(num_is_real(_self));

    num_t res;
    arb_t _res, self, other;

    arb_init(self); 
    arb_set_d(self, num_real_d(_self));    
    arb_ceil(_res, self, prec);
    res = num_from_arb(_res);
    arb_clear(self); arb_clear(_res);

    return res;
}



num_t
num_rgamma (const num_t z)
{
    assert(num_is_real(z));

    arb_t res, x;    
    arb_init(res); arb_init(x);
    arb_set_d(x, num_to_double(z));
    arb_hypgeom_rgamma(res, x,  prec);
    const double ret = arbtod(res);
    arb_clear(res); arb_clear(x);

    return new(num, ret, 0.0);
}
