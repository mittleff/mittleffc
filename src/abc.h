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
 * @file abc.h
 * @brief Interface for the Abstract Base Class (ABC).
 */
#ifndef __ABC_H__
#define __ABC_H__

#include <stdarg.h>
#include <stdio.h>

struct ABC
{
    size_t size;
    /**
     * Constructor
     */
    void* (* ctor) (void* self, va_list* app);
    /**
     * Destructor
     */
    void* (* dtor) (void* self);
    /**
     * Clone object
     */
    /* void* (* clone) (const void* self); */
};

#endif /* __ABC_H__ */
