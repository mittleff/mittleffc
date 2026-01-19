# mittleff Examples

This directory contains example programs demonstrating how to use the mittleff library.

## Building Examples

Examples are built automatically when you run:

```bash
./configure
make
```

The compiled examples will be in this directory but are not installed.

## Running Examples

After building, you can run the examples directly:

```bash
./simple
./region5_integration
./flint_integration_issue
```

## Examples

### simple.c

A basic example showing how to:
- Include mittleff headers
- Link against the library
- Call library functions

### region5_integration.c

An advanced example demonstrating:
- Computing the Mittag-Leffler function using integral representation (region 5)
- Handling edge cases with beta=1.0 that involve 0^0 expressions
- Using the custom adaptive Simpson integration to avoid numerical issues
- Working with FLINT's arbitrary-precision complex numbers (acb_t)

This example computes E_{0.5,1.0}(4.08373780+2.53485316i) and shows:
- Direct use of the algorithm module's `ml_integral_rep` function
- Proper initialization and cleanup of FLINT types
- Formatting and displaying results with error bounds

**Test case details:**
- Parameters: α=0.5, β=1.0, z=4.08373780+2.53485316i
- Expected result: -15817.811 + 54393.051i
- This case was challenging because β=1.0 creates 0^0 expressions in the integrand
- The implementation uses adaptive Simpson integration instead of FLINT's acb_calc_integrate

### flint_integration_issue.c

A diagnostic example that **demonstrates why FLINT's integration fails** for this package:
- Shows the "ulp error not defined for special value!" error in action
- Compares FLINT's acb_calc_integrate vs. custom adaptive Simpson integration
- Tests both simple (x^0) and Mittag-Leffler-like integrands
- Uses signal handling to catch FLINT's SIGABRT and continue execution
- Provides detailed explanation of the root cause

**What you'll see when running this:**
1. FLINT attempting to integrate x^0 → **FAILS** with special value error
2. Adaptive Simpson integrating x^0 → **SUCCEEDS** with correct result (π)
3. FLINT attempting ML-like integrand → **FAILS** again
4. Adaptive Simpson integrating ML-like integrand → **SUCCEEDS** with result 2.0

**Key insight:** FLINT correctly computes x^0 = 1, but marks it as a "special value"
(like 0^0). During adaptive subdivision, FLINT cannot compute ULP (unit in last place)
error bounds for special values, causing it to abort. Our adaptive Simpson integration
doesn't have this limitation.

## Using mittleff in Your Own Projects

### Using pkg-config (Recommended)

After installing mittleff, compile your program with:

```bash
gcc -o myprogram myprogram.c $(pkg-config --cflags --libs mittleff)
```

### Manual Compilation

Alternatively, specify flags manually:

```bash
gcc -o myprogram myprogram.c -I/usr/local/include -L/usr/local/lib -lutils
```

### In Your Makefile

```makefile
CFLAGS += $(shell pkg-config --cflags mittleff)
LDLIBS += $(shell pkg-config --libs mittleff)
```

### In Your CMakeLists.txt

```cmake
find_package(PkgConfig REQUIRED)
pkg_check_modules(MITTLEFFC REQUIRED mittleff)

target_include_directories(your_target PRIVATE ${MITTLEFFC_INCLUDE_DIRS})
target_link_libraries(your_target ${MITTLEFFC_LIBRARIES})
```
