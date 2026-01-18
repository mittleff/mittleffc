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
```

## Example: simple.c

A basic example showing how to:
- Include mittleff headers
- Link against the library
- Call library functions

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
