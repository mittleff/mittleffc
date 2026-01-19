AC_INIT([mittleff], [0.1], [your-email@example.com])
AM_INIT_AUTOMAKE([subdir-objects -Wall -Werror foreign])
AC_PROG_CC
AM_PROG_AR
LT_INIT

# Check for packages
PKG_CHECK_MODULES([CHECK], [check >= 0.9.4])
PKG_CHECK_MODULES([FLINT], [flint >= 3.0.1])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
 src/mittleff/Makefile
 src/utils/Makefile
 src/partition/Makefile
 src/quad/Makefile
 src/algorithm/Makefile
 tests/Makefile
 examples/Makefile
 mittleff.pc
])
AC_OUTPUT
