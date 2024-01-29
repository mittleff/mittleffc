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
 * @file algorithm.h
 * @brief Interface for the routines regarding the algorithm for each region of the complex plane.
 */
#ifndef __ALGORITHM_H__
#define __ALGORITHM_H__

#include <flint/acb.h>
//#include "num.h"

/**
 * Compute ML function in region G_0
 */
void
mittleff0 (acb_t res,
           const arb_t alpha, const arb_t beta,
           const acb_t z,
           const acb_t acc);

/**
 * Compute ML function in region G_1
 */
void
mittleff1 (acb_t res,
           const arb_t alpha, const arb_t beta,
           const acb_t z,
           const acb_t acc);
/**
 * Compute ML function in region G_2
 */
void
mittleff2 (acb_t res,
           const arb_t alpha, const arb_t beta,
           const acb_t z,
           const acb_t acc);

/**
 * Compute ML function in region G_3
 */
void
mittleff3 (acb_t res,
           const arb_t alpha, const arb_t beta,
           const acb_t z,
           const acb_t acc);

/**
 * Compute ML function in region G_4
 */
void
mittleff4 (acb_t res,
           const arb_t alpha, const arb_t beta,
           const acb_t z,
           const acb_t acc);

/**
 * Compute ML function in region G_5
 */
void
mittleff5 (acb_t res,
           const arb_t alpha, const arb_t beta,
           const acb_t z,
           const acb_t acc);

/**
 * Compute ML function in region G_6
 */
void
mittleff6 (acb_t res,
           const arb_t alpha,
           const arb_t beta,
           const acb_t z,
           const acb_t acc);

#endif /* __ALGORITHM_H__ */
