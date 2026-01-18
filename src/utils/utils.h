/**
 * @file utils.h
 * @brief Private utility functions for internal use
 *
 * This is a PRIVATE header file and is NOT installed with the library.
 * These utility functions are for internal use within libmittleff only.
 *
 * @version 0.1
 * @author Victor P. Santos
 * @copyright MIT License
 */

#ifndef ML_UTILS_PRIVATE_H
#define ML_UTILS_PRIVATE_H

#include <stddef.h>

/**
 * @brief Example utility function - adds two integers
 *
 * This is a simple utility function that demonstrates the utils module.
 *
 * @param a First integer
 * @param b Second integer
 * @return The sum of a and b
 */
int
utils_add(int a, int b);

/**
 * @brief Example utility function - multiplies two integers
 *
 * This is a simple utility function that demonstrates the utils module.
 *
 * @param a First integer
 * @param b Second integer
 * @return The product of a and b
 */
int
utils_multiply(int a, int b);

#endif /* ML_UTILS_PRIVATE_H */
