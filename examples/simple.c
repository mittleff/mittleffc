/**
 * @file simple.c
 * @brief Simple example demonstrating mittleff library usage
 *
 * This example shows how to link against and use the mittleff library.
 * Compile with: gcc -o simple simple.c $(pkg-config --cflags --libs mittleff)
 */

#include <stdio.h>
#include <mittleff/mittleff.h>

int main(void)
{
    printf("=== mittleff Library Example ===\n\n");

    printf("Calling hello() function:\n");
    int result = mittleff();

    printf("\nFunction returned: %d\n", result);
    printf("\nExample completed successfully!\n");

    return 0;
}
