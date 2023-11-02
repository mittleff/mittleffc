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
 * @file num.h
 * @brief Interface of the Abstract Data Type (ADT).
 * @details This should encapsulate any routine working with numbers.
 */
#ifndef __NUM_H__
#define __NUM_H__

#include <stdbool.h>
#include <complex.h>
#include <stdarg.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944
#endif


/**
 * This should be used in the initialization of the variable
 */
extern const void * num;

/**
 * Type associated with the class
 *
 * Objects of this type may directly be created by the user.
 * 
 */
typedef void * num_t;

void
num_print (const num_t self, const bool endline);

/**
 * Sets \p self to zero.
 */
void
num_zero (num_t self);

/**
 * Sets \p self to one.
 */
void
num_one (num_t self);

/**
 * Sets \p self to the imaginary unit.
 */
void
num_onei (num_t self);

void
num_set_num (num_t self, const num_t other);

void
num_set_d (num_t self, const double x);

void
num_set_d_d (num_t self, const double x, const double y);

/***************************************/
/* Accessors: Real and Imaginary parts */
/***************************************/

/**
 * Returns the real part of \p _self (\f$x\f$), \f$\mathrm{Re}(x)\f$.
 */
void
num_real (num_t res, const num_t self);
double
num_real_d (const num_t self);

/**
 * Returns the imaginary part of \p _self (\f$x\f$), \f$\mathrm{Im}(x)\f$.
 */
void
num_imag (num_t res, const num_t self);
double
num_imag_d (const num_t self);

/**************/
/* Predicates */
/**************/

/**
 * Returns nonzero iff \p _self is zero, within tolerance.
 */
bool
num_is_zero (const num_t self);

/**
 * Returns nonzero iff \p _self has a zero imaginary part.
 */
bool
num_is_real (const num_t self);

/****************/
/* Type casting */
/****************/

/**
 * Converts numeric types to machine's double.
 */
double
num_to_d (const num_t self);

void
num_to_d_d (double* res, const num_t self);

double complex
num_to_complex (const num_t self);

/********************/
/* Unary operations */
/********************/

/**
 * Returns the absolute value of \p _self (\f$x\f$), \f$|x|\f$.
 */
void
num_abs (num_t res, const num_t self);

/**
 * Returns the negative of \p _self (\f$x\f$), \f$-x\f$.
 */
void
num_neg (num_t res, const num_t self);

void
num_inv (num_t res, const num_t self);

/**
 * Returns the complex conjugate of \p _self (\f$x\f$), \f$x^\ast\f$.
 */
void
num_conj (num_t res, const num_t self);

void
num_ceil (num_t res, const num_t self);

/**
 * Returns the argument of \p _self (\f$x\f$), \f$\mathrm{arg}(x)\f$.
 */
void
num_arg (num_t res, const num_t self);

/**
 * Returns the square root of the number
 */
void
num_sqrt (num_t res, const num_t self);

/**
 * Returns the exponential of the number
 */
void
num_exp (num_t res, const num_t self);

/**
 * Returns the logarithm of the number
 */
void
num_log (num_t res, const num_t self);

/**
 * Returns the sine of the number
 */
void
num_sin (num_t res, const num_t self);

/**
 * Returns the hyperbolic sine of the number
 */
void
num_sinh (num_t res, const num_t self);

/**
 * Returns the cosine of the number
 */
void
num_cos (num_t res, const num_t self);

/**
 * Returns the hyperbolic cosine of the number
 */
void
num_cosh (num_t res, const num_t self);

/*********************/
/* Binary operations */
/*********************/

/**************/
/* Arithmetic */
/**************/

/**
 * returns the addition of \p _self (\f$x\f$) and \p _other (\f$y\f$), \f$x + y\f$.
 */
void
num_add (num_t res, const num_t self, const num_t other);
void
num_add_d (num_t res, const num_t self, const double other);

/**
 * returns the subtraction of \p _self (\f$x\f$) and \p _other (\f$y\f$), \f$x - y\f$.
 */
void
num_sub (num_t res, const num_t self, const num_t other);

/**
 * returns the multiplication of \p _self (\f$x\f$) and \p _other (\f$y\f$), \f$x\cdot y\f$.
 */
void
num_mul (num_t res, const num_t self, const num_t other);
void
num_mul_d (num_t res, const num_t self, const double other);

/**
 * returns the division of \p _self (\f$x\f$) and \p _other (\f$y\f$), \f$x/y\f$.
 */
void
num_div (num_t res, const num_t self, const num_t other);

/**
 *  Returns the remainder of the division of \p _self (\f$x\f$) and \p _other (\f$y\f$), \f$x/y\f$.
 */
void
num_fmod (num_t res, const num_t self, const num_t other);

/**
 * returns the exponentiation of \p _self (\f$x\f$) to \p _other (\f$y\f$), \f$x^y\f$.
 */
void
num_pow (num_t res, const num_t self, const num_t other);
void
num_pow_d (num_t res, const num_t self, const double other);

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
num_eq_d (const num_t _self, const double _other);

bool
num_gt (const num_t _self, const num_t _other);

bool
num_gt_d (const num_t _self, const double _other);

bool
num_lt (const num_t _self, const num_t _other);

bool
num_ge (const num_t _self, const num_t _other);
bool
num_ge_d (const num_t _self, const double _other);

bool
num_le (const num_t _self, const num_t _other);
bool
num_le_d (const num_t _self, const double _other);

/*********************/
/* Special functions */
/*********************/

void
num_erf (num_t res, const num_t self);

void
num_erfc (num_t res, const num_t self);

void
num_rgamma (num_t res, const num_t self);

void
num_max (num_t res, const num_t self, const num_t other);

void
num_max3 (num_t res, const num_t self, const num_t other, const num_t another);
#endif /* __NUM_H__ */
