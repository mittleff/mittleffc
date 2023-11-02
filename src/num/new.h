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
 * @file new.h
 * @brief Interface for the initialization routines of the classes.
 */
#ifndef __NEW_H__
#define __NEW_H__

#include <stddef.h>

/**
 * Attempts to allocate and initialize an object
 */
void *
new (const void* class, ...);

/**
 * Deallocate the memory allocated by the new() function
 */
void
delete (void* item);

/**
 * Copy objects
 */
/* void * */
/* clone (const void* self); */

/**
 * Size in bytes of the object
 */
size_t
size_of (const void* self);

#endif /* __NEW_H__ */
