#ifndef _NUM_H_
#define _NUM_H_

#include <complex.h>
#include <stdbool.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944
#endif

/**
 * Opaque pointer to ADT
 */
typedef struct num_instance_t * num_t;

/* Memory Management */

/**
 * Initializes the variable x for use, and sets its value to zero.
 */
num_t num_init (void);

/**
 * Clears the variable x, freeing or recycling its allocated memory.
 */
void num_clear (num_t self);

/* Input and Output */
void num_print (const num_t self);

/* Basic manipulation */

/**
 * Sets the real and imaginary parts of self to the values x and y respectively
 */
void num_set_d_d (num_t self, const double x, const double y);

/**
 * Sets self to the value of x.
 */
void num_set_d   (num_t self, const double x);

/**
 * Sets self to the value of x.
 */
void num_set_num (num_t self, const num_t x);
/* void num_zero (num_t self); */
/* void num_zero (num_t self); */

/* Accessors */
double num_real_d (const num_t self);
double num_imag_d (const num_t self);

/* /\* Predicates *\/ */
/* bool num_isnan  (num_t self); */

/* bool num_iszero (num_t self); */

/* Precision and comparisons */

/**
 * Returns true iff self is zero.
 */
bool num_is_zero (const num_t self);

/**
 * Returns true iff self equals other
 */
bool num_eq (const num_t self, const num_t other);
bool num_eq_d (const num_t self, const double other);

/**
 * Returns true iff self is not equal to other
 */
bool num_ne (const num_t self, const num_t other);

bool num_gt (const num_t self, const num_t other);
bool num_gt_d (const num_t self, const double other);
bool num_ge (const num_t self, const num_t other);
bool num_ge_d (const num_t self, const double other);

bool num_lt (const num_t self, const num_t other);
bool num_le (const num_t self, const num_t other);
bool num_le_d (const num_t self, const double other);

bool num_is_real (const num_t self);

void num_max2 (num_t res, num_t x, num_t y);
void num_max3 (num_t res, num_t x, num_t y, num_t z);
/* Arithmetic */

/**
 * Sets res to the negation of x.
 */
void num_neg (num_t res, const num_t x);

/**
 *     Sets z to the complex conjugate of x.
 */
void num_conj (num_t res, const num_t x);

/**
 * Sets z to the sum of x and y.
 */
void num_add (num_t res, const num_t x, const num_t y);

/**
 * Sets z to the difference of x and y.
 */
void num_sub (num_t res, const num_t x, const num_t y);

void num_mul (num_t res, const num_t x, const num_t y);
void num_mul_d (num_t res, const num_t x, const double y);

void num_div (num_t res, const num_t x, const num_t y);

void num_pow (num_t res, const num_t x, const num_t y);
void num_pow_d (num_t res, const num_t x, const double y);

double num_to_d (const num_t self);
void num_to_d_d (double* res, const num_t x);
complex double num_to_complex (const num_t self);

void
num_inv (num_t res, const num_t x);
void num_abs (num_t res, const num_t x);
void num_ceil (num_t res, const num_t x);
void num_sqrt (num_t res, const num_t x);
void num_arg (num_t res, const num_t x);
void num_fmod (num_t res, const num_t x, const num_t y);
void num_exp (num_t res, const num_t x);
void num_log (num_t res, const num_t x);
void num_sin (num_t res, const num_t x);
void num_sinh (num_t res, const num_t x);
void num_cos (num_t res, const num_t x);
void num_cosh (num_t res, const num_t x);
void num_rgamma (num_t res, const num_t x);
void num_erfc (num_t res, const num_t x);
void num_arg (num_t res, num_t x);
#endif /* _NUM_H_ */
