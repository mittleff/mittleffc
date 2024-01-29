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
 * @file integrate.h
 * @brief Interface for the routines integrating the functions on complex plane.
 */

#ifndef __INTEGRATE_H__
#define __INTEGRATE_H__

#include <flint/arb.h>
#include <flint/acb.h>
//#include "quad.h

void
A (acb_t res,
   const acb_t z,
   const arb_t alpha,
   const arb_t beta,
   const arb_t x);

/* void */
/* A (num_t res, */
/*    const num_t z, */
/*    const num_t alpha, */
/*    const num_t beta, */
/*    const num_t x); */

/* void */
/* B (num_t res, */
/*    const num_t r, */
/*    const num_t alpha, */
/*    const num_t beta, */
/*    const num_t z, */
/*    const num_t phi); */

/* void */
/* integrate_B (acb_t res, */
/*              const num_t alpha, */
/*              const num_t beta, */
/*              const num_t z, */
/*              const num_t phi, */
/*              const num_t from, */
/*              const num_t to, */
/*              const num_t acc); */

/* void */
/* integrate_C (num_t res, */
/*              const num_t alpha, */
/*              const num_t beta, */
/*              const num_t z, */
/*              const num_t rho, */
/*              const num_t from, */
/*              const num_t to); */

#endif  /* __INTEGRATE_H__ */
