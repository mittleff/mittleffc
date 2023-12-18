#ifndef __INTEGRATE_H__
#define __INTEGRATE_H__

#include "num.h"

void
A (num_t res,
   const num_t z,
   const num_t alpha,
   const num_t beta,
   const num_t x);

void
B (num_t res,
   const num_t r,
   const num_t alpha,
   const num_t beta,
   const num_t z,
   const num_t phi);

void
integrate_B (num_t res,
             const num_t alpha,
             const num_t beta,
             const num_t z,
             const num_t phi,
             const num_t from,
             const num_t to);

void
integrate_C (num_t res,
             const num_t alpha,
             const num_t beta,
             const num_t z,
             const num_t rho,
             const num_t from,
             const num_t to);

#endif  /* __INTEGRATE_H__ */
