# Examples are not installed, just built for demonstration
noinst_PROGRAMS = simple.x region5_integration.x flint_integration_issue.x diagnose_z_zero.x

simple_x_SOURCES = simple.c
simple_x_CFLAGS = -I$(top_srcdir)/src @FLINT_CFLAGS@
simple_x_LDADD = $(top_builddir)/src/mittleff/libmittleff.la @FLINT_LIBS@ -lflint

region5_integration_x_SOURCES = region5_integration.c
region5_integration_x_CFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/algorithm @FLINT_CFLAGS@
region5_integration_x_LDADD = $(top_builddir)/src/algorithm/libalgorithm.la \
                            $(top_builddir)/src/partition/libpartition.la \
                            $(top_builddir)/src/quad/libquad.la \
                            @FLINT_LIBS@ -lflint

flint_integration_issue_x_SOURCES = flint_integration_issue.c
flint_integration_issue_x_CFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/quad @FLINT_CFLAGS@
flint_integration_issue_x_LDADD = $(top_builddir)/src/quad/libquad.la \
                                 @FLINT_LIBS@ -lflint

diagnose_z_zero_x_SOURCES = diagnose_z_zero.c
diagnose_z_zero_x_CFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/mittleff @FLINT_CFLAGS@
diagnose_z_zero_x_LDADD = $(top_builddir)/src/mittleff/libmittleff.la \
                            @FLINT_LIBS@ -lflint

# Include examples in distribution
EXTRA_DIST = README.md
