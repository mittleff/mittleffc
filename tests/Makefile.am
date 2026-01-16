TESTS = check_partition.x check_algorithm.x check_z_zero.x check_mittleff.x check_recursion.x
AM_TESTS_ENVIRONMENT=CK_VERBOSITY=verbose

check_PROGRAMS = check_partition.x check_mittleff.x check_recursion.x check_z_zero.x check_algorithm.x

check_partition_x_SOURCES = check_partition.c
check_partition_x_CFLAGS = @CHECK_CFLAGS@ @FLINT_CFLAGS@ -I$(top_srcdir)/src
check_partition_x_LDADD = $(top_builddir)/src/partition/libpartition.la @CHECK_LIBS@ @FLINT_LIBS@ -lflint

check_algorithm_x_SOURCES = check_algorithm.c
check_algorithm_x_CFLAGS = @CHECK_CFLAGS@ @FLINT_CFLAGS@ -I$(top_srcdir)/src
check_algorithm_x_LDADD = $(top_builddir)/src/algorithm/libalgorithm.la $(top_builddir)/src/partition/libpartition.la $(top_builddir)/src/quad/libquad.la @CHECK_LIBS@ @FLINT_LIBS@ -lflint

check_mittleff_x_SOURCES = check_mittleff.c
check_mittleff_x_CFLAGS = @CHECK_CFLAGS@ @FLINT_CFLAGS@ -I$(top_srcdir)/src
check_mittleff_x_LDADD = $(top_builddir)/src/mittleff/libmittleff.la @CHECK_LIBS@ @FLINT_LIBS@ -lflint

check_recursion_x_SOURCES = check_recursion.c
check_recursion_x_CFLAGS = @CHECK_CFLAGS@ @FLINT_CFLAGS@ -I$(top_srcdir)/src
check_recursion_x_LDADD = $(top_builddir)/src/mittleff/libmittleff.la @CHECK_LIBS@ @FLINT_LIBS@ -lflint

# check_utils_SOURCES = check_utils.c
# check_utils_CFLAGS = @CHECK_CFLAGS@ -I$(top_srcdir)/src
# check_utils_LDADD = $(top_builddir)/src/utils/libutils.la @CHECK_LIBS@

check_z_zero_x_SOURCES = check_z_zero.c
check_z_zero_x_CFLAGS = @CHECK_CFLAGS@ @FLINT_CFLAGS@ -I$(top_srcdir)/src
check_z_zero_x_LDADD = $(top_builddir)/src/mittleff/libmittleff.la @CHECK_LIBS@ @FLINT_LIBS@ -lflint
