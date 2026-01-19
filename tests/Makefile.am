TESTS = check_utils check_partition check_algorithm check_mittleff check_recursion
AM_TESTS_ENVIRONMENT=CK_VERBOSITY=verbose

check_PROGRAMS = check_mittleff check_recursion check_utils check_partition check_algorithm

check_mittleff_SOURCES = check_mittleff.c
check_mittleff_CFLAGS = @CHECK_CFLAGS@ @FLINT_CFLAGS@ -I$(top_srcdir)/src
check_mittleff_LDADD = $(top_builddir)/src/mittleff/libmittleff.la @CHECK_LIBS@ @FLINT_LIBS@ -lflint

check_recursion_SOURCES = check_recursion.c
check_recursion_CFLAGS = @CHECK_CFLAGS@ @FLINT_CFLAGS@ -I$(top_srcdir)/src
check_recursion_LDADD = $(top_builddir)/src/mittleff/libmittleff.la @CHECK_LIBS@ @FLINT_LIBS@ -lflint

check_utils_SOURCES = check_utils.c
check_utils_CFLAGS = @CHECK_CFLAGS@ -I$(top_srcdir)/src
check_utils_LDADD = $(top_builddir)/src/utils/libutils.la @CHECK_LIBS@

check_partition_SOURCES = check_partition.c
check_partition_CFLAGS = @CHECK_CFLAGS@ @FLINT_CFLAGS@ -I$(top_srcdir)/src
check_partition_LDADD = $(top_builddir)/src/partition/libpartition.la $(top_builddir)/src/utils/libutils.la @CHECK_LIBS@ @FLINT_LIBS@ -lflint

check_algorithm_SOURCES = check_algorithm.c
check_algorithm_CFLAGS = @CHECK_CFLAGS@ @FLINT_CFLAGS@ -I$(top_srcdir)/src
check_algorithm_LDADD = $(top_builddir)/src/algorithm/libalgorithm.la $(top_builddir)/src/partition/libpartition.la $(top_builddir)/src/quad/libquad.la $(top_builddir)/src/utils/libutils.la @CHECK_LIBS@ @FLINT_LIBS@ -lflint
