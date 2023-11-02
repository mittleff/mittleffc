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

#include "num.h"

/**
 * Compute ML function in region G_0
 */
num_t
mittleff0 (const num_t alpha,
           const num_t beta,
           const num_t z,
           const num_t acc);

#endif /* __ALGORITHM_H__ */
