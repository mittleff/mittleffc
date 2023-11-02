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
 * @file num.h
 * @brief Interface of the Abstract Data Type (ADT).
 * @details Encapsulate arbitrary-precision library.
 */
#ifndef __NUM_H__
#define __NUM_H__

#include <stdbool.h>

/**
 * This should be used in the initialization of the variable
 */
extern const void* num;

/**
 * Type associated with the class
 *
 * Objects of this type may directly be created by the user.
 * 
 */
typedef void* num_t;


/* num_t */
/* num_zero (void); */

/***************************************/
/* Accessors: Real and Imaginary parts */
/***************************************/

/**
 * Returns the real part of \p _self (\f$x\f$), \f$\mathrm{Re}(x)\f$.
 */
num_t
num_real (const num_t _self);

/**
 * Returns the imaginary part of \p _self (\f$x\f$), \f$\mathrm{Im}(x)\f$.
 */
num_t
num_imag (const num_t _self);

/**
 * Returns the real part of \p _self (\f$x\f$), \f$\mathrm{Re}(x)\f$, as machine double.
 */
double
num_real_d (const num_t _self);

/**
 * Returns the imaginary part of \p _self (\f$x\f$), \f$\mathrm{Im}(x)\f$, as machine double.
 */
double
num_imag_d (const num_t _self);

/**************/
/* Predicates */
/**************/

/**
 * Returns nonzero iff \p _self is zero, within tolerance.
 */
bool
num_is_zero (const num_t _self);

/**
 * Returns nonzero iff \p _self has a zero imaginary part.
 */
bool
num_is_real (const num_t _self);

/****************/
/* Type casting */
/****************/

/**
 * Converts numeric types to machine's double.
 */
double
num_to_double (const num_t _self);

/**
 * Converts numeric types to machine's complex double.
 */
/* double complex */
/* num_to_complex (const num_t _self); */

/* void */
/* num_to_pair (double* res, const num_t _self); */

/* /\********************\/ */
/* /\* Unary operations *\/ */
/* /\********************\/ */

/* /\** */
/*  * Returns the square of absolute value of \p _self (\f$x\f$), \f$|x|\f$. */
/*  *\/ */
/* num_t */
/* num_abs2 (const num_t _self); */

/**
 * Returns the absolute value of \p _self (\f$x\f$), \f$|x|\f$.
 */
num_t
num_abs (const num_t _self);

/* /\** */
/*  * Returns the negative of \p _self (\f$x\f$), \f$-x\f$. */
/*  *\/ */
/* num_t */
/* num_negative (const num_t _self); */

/* /\** */
/*  * Returns the complex conjugate of \p _self (\f$x\f$), \f$x^\ast\f$. */
/*  *\/ */
/* num_t */
/* num_conjugate (const num_t _self); */

/**
 * Returns the argument of \p _self (\f$x\f$), \f$\mathrm{arg}(x)\f$.
 */
num_t
num_arg (const num_t _self);

/* /\** */
/*  * Returns the square root of the number */
/*  *\/ */
/* num_t */
/* num_sqrt (const num_t _self); */

/**
 * Returns the exponential of the number
 */
num_t
num_exp (const num_t _self);

/**
 * Returns the logarithm of the number
 */
num_t
num_log (const num_t _self);

/* /\** */
/*  * Returns the sine of the number */
/*  *\/ */
/* num_t */
/* num_sin (const num_t _self); */

/* /\** */
/*  * Returns the cosine of the number */
/*  *\/ */
/* num_t */
/* num_cos (const num_t _self); */

/* /\*********************\/ */
/* /\* Binary operations *\/ */
/* /\*********************\/ */

/* /\**************\/ */
/* /\* Arithmetic *\/ */
/* /\**************\/ */

/**
 * returns the addition of \p _self (\f$x\f$) and \p _other (\f$y\f$), \f$x + y\f$.
 */
num_t
num_add (const num_t _self, const num_t _other);

/**
 * returns the subtraction of \p _self (\f$x\f$) and \p _other (\f$y\f$), \f$x - y\f$.
 */
num_t
num_sub (const num_t _self, const num_t _other);

/**
 * returns the multiplication of \p _self (\f$x\f$) and \p _other (\f$y\f$), \f$x\cdot y\f$.
 */
num_t
num_mul (const num_t _self, const num_t _other);

/**
 * returns the division of \p _self (\f$x\f$) and \p _other (\f$y\f$), \f$x/y\f$.
 */
num_t
num_div (const num_t _self, const num_t _other);

/* /\** */
/*  *  Returns the remainder of the division of \p _self (\f$x\f$) and \p _other (\f$y\f$), \f$x/y\f$. */
/*  *\/ */
/* num_t */
/* num_fmod (const num_t _self, const num_t _other); */

/**
 * returns the exponentiation of \p _self (\f$x\f$) to \p _other (\f$y\f$), \f$x^y\f$.
 */
num_t
num_pow (const num_t _self, const num_t _other);

/***********/
/* Logical */
/***********/

/**
 * This function determines whether \p _self (\f$x\f$) and \p _other (\f$y\f$),
 * are approximately equal to a relative accuracy.
 */
bool
num_eq (const num_t _self, const num_t _other);

bool
num_gt (const num_t _self, const num_t _other);

bool
num_lt (const num_t _self, const num_t _other);

bool
num_ge (const num_t _self, const num_t _other);

bool
num_le (const num_t _self, const num_t _other);

/*********************/
/* Special functions */
/*********************/
num_t
num_rgamma(const num_t _self);

num_t
num_max(const num_t _self, const num_t _other);

num_t
num_ceil (const num_t _self);

#endif /* __NUM_H__ */
