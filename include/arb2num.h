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
 * @brief Interface for the routines to interact with arblib.
 */
#ifndef __ARB2NUM_H__
#define __ARB2NUM_H__

#include <arb.h>
#include <acb.h>
#include "num.h"

double
arbtod (const arb_t x);

num_t
num_from_acb (const acb_t x);

num_t
num_from_arb (const arb_t x);


#endif /* __ARB2NUM_H__ */
