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
 * @file partition.h
 * @brief Interface for the routines regarding the partition of the complex plane.
 */
#ifndef __PARTITION_H__
#define __PARTITION_H__

#include "num.h"
#include <stdbool.h>

/**
 * Size in bytes of the object
 */
bool
in_region_G0 (const num_t z);

/**
 * Size in bytes of the object
 */
bool
in_region_G1 (const num_t z, const num_t alpha, const num_t acc);

/**
 * Size in bytes of the object
 */
bool
in_region_G2 (const num_t z, const num_t alpha, const num_t acc);

/**
 * Size in bytes of the object
 */
bool
in_region_G3 (const num_t z, const num_t alpha, const num_t acc);

/**
 * Size in bytes of the object
 */
bool
in_region_G4 (const num_t z, const num_t alpha, const num_t acc);

/**
 * Size in bytes of the object
 */
bool
in_region_G5 (const num_t z, const num_t alpha, const num_t acc);

/**
 * Size in bytes of the object
 */
bool
in_region_G6 (const num_t z, const num_t alpha, const num_t acc);

#endif /* __PARTITION_H__ */
