/**
 * @file mittleff.h
 * @brief Mittleff C Library - Utility Functions
 *
 * This header provides utility functions for the Mittleff C library.
 * All public API functions are prefixed with the library name to avoid
 * namespace conflicts.
 *
 * @version 0.1
 * @author Victor P. Santos
 * @copyright MIT License
 */

#ifndef ML_MITTLEFF_H
#define ML_MITTLEFF_H

#include <stdio.h>

/* Symbol visibility macros */
#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_MITTLEFF
    #define MITTLEFF_EXPORT __declspec(dllexport)
  #else
    #define MITTLEFF_EXPORT __declspec(dllimport)
  #endif
#else
  #if __GNUC__ >= 4
    #define MITTLEFF_EXPORT __attribute__((visibility("default")))
  #else
    #define MITTLEFF_EXPORT
  #endif
#endif

/**
 * @brief Prints a hello message with library version information
 *
 * This function prints two lines to standard output:
 * 1. "Hello World!"
 * 2. The library package string with version information
 *
 * @return Always returns 5
 *
 * @example
 * @code
 * #include <mittleff/mittleff.h>
 *
 * int main(void) {
 *     int result = hello();
 *     printf("Function returned: %d\n", result);
 *     return 0;
 * }
 * @endcode
 *
 * @note This is a demonstration function for the library.
 * @since 0.1
 */
MITTLEFF_EXPORT int
mittleff (void);

#endif /* ML_MITTLEFF_H */
