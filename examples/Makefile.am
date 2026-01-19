# Examples are not installed, just built for demonstration
noinst_PROGRAMS = simple region5_integration flint_integration_issue

simple_SOURCES = simple.c
simple_CFLAGS = -I$(top_srcdir)/src
simple_LDADD = $(top_builddir)/src/mittleff/libmittleff.la

region5_integration_SOURCES = region5_integration.c
region5_integration_CFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/algorithm @FLINT_CFLAGS@
region5_integration_LDADD = $(top_builddir)/src/algorithm/libalgorithm.la \
                            $(top_builddir)/src/partition/libpartition.la \
                            $(top_builddir)/src/quad/libquad.la \
                            $(top_builddir)/src/utils/libutils.la \
                            @FLINT_LIBS@ -lflint

flint_integration_issue_SOURCES = flint_integration_issue.c
flint_integration_issue_CFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/quad @FLINT_CFLAGS@
flint_integration_issue_LDADD = $(top_builddir)/src/quad/libquad.la \
                                 @FLINT_LIBS@ -lflint

# Include examples in distribution
EXTRA_DIST = README.md
