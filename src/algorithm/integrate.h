#ifndef __INTEGRATE_H__
#define __INTEGRATE_H__

#include "num.h"

num_t
A (const num_t z,
   const num_t alpha,
   const num_t beta,
   const num_t x);

num_t
B (const num_t r,
   const num_t alpha,
   const num_t beta,
   const num_t z,
   const num_t phi);

num_t
integrate_B (const num_t alpha,
             const num_t beta,
             const num_t z,
             const num_t phi,
             const num_t from,
             const num_t to);

num_t
integrate_C (const num_t alpha,
             const num_t beta,
             const num_t z,
             const num_t rho,
             const num_t from,
             const num_t to);

#endif  /* __INTEGRATE_H__ */
